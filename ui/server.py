"""Dependency-free local server for configuring and visualizing NegPar runs."""

from __future__ import annotations

import argparse
import collections
import dataclasses
import json
import math
import mimetypes
import os
from pathlib import Path
import re
import subprocess
import threading
import time
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any
from urllib.parse import parse_qs, urlparse
import webbrowser


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
STATIC_ROOT = Path(__file__).resolve().parent / "static"
STEP_PATTERN = re.compile(r"^step\s+(\d+)/(\d+)\s*$")
SEED_PATTERN = re.compile(r"^seed\s*=\s*(\d+)\s*$")
SNAPSHOT_PATTERN = re.compile(r"^(?P<name>.+)_(?P<index>\d{3})\.txt$")


class RequestError(ValueError):
    """A client-visible request validation error."""


def _inside(parent: Path, child: Path) -> bool:
    return child == parent or parent in child.parents


def resolve_output_directory(root: Path, value: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise RequestError("Output directory is required.")
    raw = Path(value.strip())
    resolved = (raw if raw.is_absolute() else root / raw).resolve()
    if not _inside(root.resolve(), resolved) or resolved == root.resolve():
        raise RequestError("Output directory must be inside the repository.")
    return resolved


def find_executable(root: Path, override: str | None = None) -> Path:
    if override:
        candidate = Path(override)
        return (candidate if candidate.is_absolute() else root / candidate).resolve()

    candidates = (
        root / "build/release/Release/negpar_inhomo.exe",
        root / "build/debug/Debug/negpar_inhomo.exe",
        root / "x64/Release/negpar_inhomo.exe",
        root / "x64/Debug/negpar_inhomo.exe",
        root / "build/release/negpar_inhomo",
        root / "build/debug/negpar_inhomo",
        root / "out",
    )
    return next((path.resolve() for path in candidates if path.is_file()),
                candidates[0].resolve())


def _integer(payload: dict[str, Any], name: str, minimum: int,
             maximum: int) -> int:
    value = payload.get(name)
    if isinstance(value, bool):
        raise RequestError(f"{name.capitalize()} must be an integer.")
    try:
        parsed = int(value)
    except (TypeError, ValueError) as error:
        raise RequestError(f"{name.capitalize()} must be an integer.") from error
    if str(value).strip() != str(parsed) and not isinstance(value, int):
        raise RequestError(f"{name.capitalize()} must be an integer.")
    if not minimum <= parsed <= maximum:
        raise RequestError(
            f"{name.capitalize()} must be between {minimum} and {maximum}.")
    return parsed


@dataclasses.dataclass(frozen=True)
class RunRequest:
    seed: int | None
    threads: int
    steps: int
    output_directory: Path
    output_label: str
    command: tuple[str, ...]


def validate_run_request(payload: dict[str, Any], root: Path,
                         executable: Path) -> RunRequest:
    if not executable.is_file():
        raise RequestError(
            "Simulation executable was not found. Build the Release or Debug preset first.")

    seed_mode = payload.get("seedMode", "fixed")
    if seed_mode not in {"fixed", "random"}:
        raise RequestError("Seed mode must be fixed or random.")
    seed = (_integer(payload, "seed", 0, 2**32 - 1)
            if seed_mode == "fixed" else None)
    threads = _integer(payload, "threads", 1, 2**31 - 1)
    steps = _integer(payload, "steps", 1, 2**31 - 1)
    output_label = str(payload.get("outputDirectory", "")).strip()
    output_directory = resolve_output_directory(root, output_label)
    if output_directory.exists() and not output_directory.is_dir():
        raise RequestError("Output directory points to an existing file.")
    if output_directory.exists() and any(output_directory.iterdir()):
        raise RequestError(
            "Output directory is not empty. Choose a new directory to preserve run isolation.")

    command = [str(executable)]
    if seed is not None:
        command.extend(("--seed", str(seed)))
    command.extend(("--threads", str(threads), "--steps", str(steps),
                    "--output-dir", str(output_directory)))
    return RunRequest(seed, threads, steps, output_directory, output_label,
                      tuple(command))


def read_numeric_file(path: Path) -> list[float]:
    before = path.stat()
    text = path.read_text(encoding="utf-8")
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise OSError("Output changed while it was being read.")
    values = [float(token) for token in text.split()]
    if any(not math.isfinite(value) for value in values):
        raise ValueError(f"{path.name} contains NaN or Inf.")
    return values


def read_metadata(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    metadata: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        key, separator, value = line.partition(" ")
        if separator:
            metadata[key] = value
    return metadata


class ResultReader:
    SPATIAL_METRICS = {
        "density": (("rho", "rho"), ("rhoF", "rhoF"), ("rhoM", "rhoM")),
        "velocity": (("u1", "u1"), ("u1F", "u1F"), ("u1M", "u1M")),
        "temperature": (("Tprt", "Tprt"), ("TprtF", "TprtF"),
                        ("TprtM", "TprtM")),
        "electric": (("elecfield", "elecfield"),
                     ("elecfield_F", "elecfield_F")),
    }
    DIAGNOSTICS = {
        "conservation": (("total_energy", "total_energy"),
                         ("total_energy_F", "total_energy_F"),
                         ("elec_energy", "elec_energy")),
        "particles": (("Np_rec", "Np"), ("Nn_rec", "Nn"),
                      ("Nf_rec", "Nf")),
    }

    def __init__(self, output_directory: Path | None):
        self.output_directory = output_directory
        self.warnings: list[str] = []

    def _read(self, filename: str) -> list[float]:
        if self.output_directory is None:
            return []
        path = self.output_directory / filename
        if not path.is_file():
            return []
        try:
            return read_numeric_file(path)
        except (OSError, UnicodeError, ValueError) as error:
            self.warnings.append(f"{filename}: {error}")
            return []

    def snapshots(self, metric: str) -> list[dict[str, Any]]:
        if self.output_directory is None or not self.output_directory.is_dir():
            return []
        names = self.SPATIAL_METRICS.get(metric, self.SPATIAL_METRICS["density"])
        indices: set[int] = set()
        has_final = False
        for base, _ in names:
            if (self.output_directory / f"{base}.txt").is_file():
                has_final = True
            for path in self.output_directory.glob(f"{base}_???.txt"):
                match = SNAPSHOT_PATTERN.match(path.name)
                if match and match.group("name") == base:
                    indices.add(int(match.group("index")))

        times = self._read("time_dist.txt")
        result = [{
            "id": index,
            "label": (f"t = {times[index]:.3f}" if index < len(times)
                      else f"snapshot {index:03d}"),
        } for index in sorted(indices)]
        if has_final:
            result.append({"id": "final", "label": "final"})
        return result

    def spatial(self, metric: str, snapshot: str | int | None) -> dict[str, Any]:
        metric = metric if metric in self.SPATIAL_METRICS else "density"
        snapshots = self.snapshots(metric)
        valid_ids = [item["id"] for item in snapshots]
        if snapshot is None or snapshot not in valid_ids:
            snapshot = valid_ids[-1] if valid_ids else None

        suffix = "" if snapshot == "final" else (
            f"_{int(snapshot):03d}" if snapshot is not None else "")
        series = []
        for base, label in self.SPATIAL_METRICS[metric]:
            values = self._read(f"{base}{suffix}.txt") if snapshot is not None else []
            if values:
                series.append({"name": label, "values": values})

        longest = max((len(item["values"]) for item in series), default=0)
        x_values = self._read("x.txt")
        if len(x_values) != longest:
            x_values = [float(index) for index in range(longest)]
        for item in series:
            item["values"] = item["values"][:longest]

        return {
            "metric": metric,
            "snapshot": snapshot,
            "snapshots": snapshots,
            "x": x_values,
            "series": series,
        }

    def diagnostics(self) -> dict[str, list[dict[str, Any]]]:
        result: dict[str, list[dict[str, Any]]] = {}
        for group, definitions in self.DIAGNOSTICS.items():
            result[group] = []
            for filename, label in definitions:
                values = self._read(f"{filename}.txt")
                if values:
                    result[group].append({"name": label, "values": values})
        return result

    def payload(self, metric: str = "density",
                snapshot: str | int | None = None) -> dict[str, Any]:
        spatial = self.spatial(metric, snapshot)
        return {
            **spatial,
            "diagnostics": self.diagnostics(),
            "metadata": (read_metadata(self.output_directory / "run_metadata.txt")
                         if self.output_directory else {}),
            "warnings": self.warnings,
        }


class SimulationManager:
    def __init__(self, root: Path, executable: Path):
        self.root = root.resolve()
        self.executable = executable.resolve()
        self._lock = threading.RLock()
        self._process: subprocess.Popen[str] | None = None
        self._consumer: threading.Thread | None = None
        self._logs: collections.deque[dict[str, str]] = collections.deque(maxlen=600)
        self._status = "ready"
        self._step = 0
        self._steps = 0
        self._started_at: float | None = None
        self._finished_at: float | None = None
        self._return_code: int | None = None
        self._output_directory: Path | None = None
        self._output_label = ""
        self._seed: int | None = None
        self._threads = 0
        self._command: tuple[str, ...] = ()
        self._stop_requested = False
        self._error = ""

    @property
    def output_directory(self) -> Path | None:
        with self._lock:
            return self._output_directory

    def _log(self, text: str) -> None:
        self._logs.append({"time": time.strftime("%H:%M:%S"), "text": text})

    def start(self, payload: dict[str, Any]) -> dict[str, Any]:
        request = validate_run_request(payload, self.root, self.executable)
        with self._lock:
            if self._process is not None and self._process.poll() is None:
                raise RequestError("A simulation is already running.")
            if self._consumer is not None and self._consumer.is_alive():
                raise RequestError("The previous simulation is still finalizing. Try again shortly.")
            self._logs.clear()
            self._status = "starting"
            self._step = 0
            self._steps = request.steps
            self._started_at = time.time()
            self._finished_at = None
            self._return_code = None
            self._output_directory = request.output_directory
            self._output_label = request.output_label
            self._seed = request.seed
            self._threads = request.threads
            self._command = request.command
            self._stop_requested = False
            self._error = ""
            self._log("Launching simulation...")
            try:
                self._process = subprocess.Popen(
                    request.command,
                    cwd=self.root,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    bufsize=1,
                )
            except OSError as error:
                self._status = "failed"
                self._error = str(error)
                self._finished_at = time.time()
                self._log(f"Unable to launch simulation: {error}")
                raise RequestError(f"Unable to launch simulation: {error}") from error
            self._status = "running"
            self._consumer = threading.Thread(target=self._consume,
                                              name="negpar-output",
                                              daemon=True)
            self._consumer.start()
        return self.snapshot()

    def _consume(self) -> None:
        process = self._process
        if process is None or process.stdout is None:
            return
        for raw_line in process.stdout:
            line = raw_line.rstrip("\r\n")
            with self._lock:
                self._log(line)
                step_match = STEP_PATTERN.match(line)
                if step_match:
                    self._step = int(step_match.group(1))
                    self._steps = int(step_match.group(2))
                seed_match = SEED_PATTERN.match(line)
                if seed_match:
                    self._seed = int(seed_match.group(1))
        return_code = process.wait()
        with self._lock:
            self._return_code = return_code
            self._finished_at = time.time()
            if self._stop_requested:
                self._status = "stopped"
                self._log("Simulation stopped by user.")
            elif return_code == 0:
                self._status = "completed"
                self._step = self._steps
            else:
                self._status = "failed"
                self._error = f"Simulation exited with code {return_code}."

    def stop(self) -> dict[str, Any]:
        with self._lock:
            if self._process is None or self._process.poll() is not None:
                raise RequestError("No simulation is currently running.")
            self._stop_requested = True
            self._status = "stopping"
            self._log("Stop requested...")
            self._process.terminate()
        return self.snapshot()

    def snapshot(self) -> dict[str, Any]:
        with self._lock:
            now = self._finished_at or time.time()
            elapsed = max(0.0, now - self._started_at) if self._started_at else 0.0
            output_directory = self._output_directory
            file_count = (sum(1 for path in output_directory.glob("*.txt")
                              if path.is_file())
                          if output_directory and output_directory.is_dir() else 0)
            return {
                "status": self._status,
                "step": self._step,
                "steps": self._steps,
                "elapsedSeconds": elapsed,
                "returnCode": self._return_code,
                "outputDirectory": self._output_label,
                "outputDirectoryAbsolute": (str(output_directory)
                                            if output_directory else ""),
                "outputFiles": file_count,
                "seed": self._seed,
                "threads": self._threads,
                "command": list(self._command),
                "logs": list(self._logs),
                "error": self._error,
            }

    def open_output_directory(self) -> None:
        with self._lock:
            directory = self._output_directory
        if directory is None or not directory.is_dir():
            raise RequestError("The output directory does not exist yet.")
        if os.name == "nt":
            os.startfile(directory)  # type: ignore[attr-defined]
        elif os.uname().sysname == "Darwin":
            subprocess.Popen(("open", str(directory)))
        else:
            subprocess.Popen(("xdg-open", str(directory)))


def default_output_label() -> str:
    return time.strftime("result/ui-run-%Y%m%d-%H%M%S")


class StudioHandler(BaseHTTPRequestHandler):
    manager: SimulationManager

    def log_message(self, format: str, *args: Any) -> None:
        return

    def _json(self, payload: Any, status: HTTPStatus = HTTPStatus.OK) -> None:
        body = json.dumps(payload, separators=(",", ":"), ensure_ascii=False).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def _body(self) -> dict[str, Any]:
        try:
            length = int(self.headers.get("Content-Length", "0"))
        except ValueError as error:
            raise RequestError("Invalid content length.") from error
        if length <= 0 or length > 64 * 1024:
            raise RequestError("Request body must contain at most 64 KiB of JSON.")
        try:
            payload = json.loads(self.rfile.read(length))
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise RequestError("Request body is not valid JSON.") from error
        if not isinstance(payload, dict):
            raise RequestError("Request body must be a JSON object.")
        return payload

    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        if parsed.path == "/api/config":
            executable = self.manager.executable
            self._json({
                "executable": str(executable),
                "executableAvailable": executable.is_file(),
                "build": ("release" if "release" in str(executable).lower()
                          else "debug"),
                "repository": str(self.manager.root),
                "defaults": {
                    "seedMode": "fixed",
                    "seed": 20260809,
                    "threads": max(1, (os.cpu_count() or 2) - 2),
                    "steps": 100,
                    "outputDirectory": default_output_label(),
                },
            })
            return
        if parsed.path == "/api/status":
            self._json(self.manager.snapshot())
            return
        if parsed.path == "/api/results":
            query = parse_qs(parsed.query)
            metric = query.get("metric", ["density"])[0]
            raw_snapshot = query.get("snapshot", [None])[0]
            snapshot: str | int | None = raw_snapshot
            if raw_snapshot not in {None, "final"}:
                try:
                    snapshot = int(raw_snapshot)
                except ValueError:
                    snapshot = None
            reader = ResultReader(self.manager.output_directory)
            self._json(reader.payload(metric, snapshot))
            return
        self._serve_static(parsed.path)

    def do_POST(self) -> None:
        try:
            if self.path == "/api/runs":
                self._json(self.manager.start(self._body()), HTTPStatus.ACCEPTED)
                return
            if self.path == "/api/runs/stop":
                self._json(self.manager.stop(), HTTPStatus.ACCEPTED)
                return
            if self.path == "/api/open-output":
                self.manager.open_output_directory()
                self._json({"opened": True})
                return
            self._json({"error": "Not found."}, HTTPStatus.NOT_FOUND)
        except RequestError as error:
            self._json({"error": str(error)}, HTTPStatus.BAD_REQUEST)

    def _serve_static(self, request_path: str) -> None:
        relative = "index.html" if request_path == "/" else request_path.lstrip("/")
        target = (STATIC_ROOT / relative).resolve()
        if not _inside(STATIC_ROOT.resolve(), target) or not target.is_file():
            self.send_error(HTTPStatus.NOT_FOUND)
            return
        content = target.read_bytes()
        mime_type = mimetypes.guess_type(target.name)[0] or "application/octet-stream"
        self.send_response(HTTPStatus.OK)
        self.send_header("Content-Type", f"{mime_type}; charset=utf-8")
        self.send_header("Content-Length", str(len(content)))
        self.send_header("Cache-Control", "no-cache")
        self.end_headers()
        self.wfile.write(content)


def create_server(host: str, port: int, executable: str | None = None,
                  root: Path = REPOSITORY_ROOT) -> ThreadingHTTPServer:
    manager = SimulationManager(root, find_executable(root, executable))
    handler = type("ConfiguredStudioHandler", (StudioHandler,), {"manager": manager})
    return ThreadingHTTPServer((host, port), handler)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--executable", help="Path to negpar_inhomo executable")
    parser.add_argument("--no-open", action="store_true",
                        help="Do not open the dashboard in the default browser")
    arguments = parser.parse_args()

    server = create_server(arguments.host, arguments.port, arguments.executable)
    url = f"http://{arguments.host}:{server.server_port}"
    executable = server.RequestHandlerClass.manager.executable
    print(f"NegPar Simulation Studio: {url}")
    print(f"Executable: {executable}")
    if not executable.is_file():
        print("Warning: executable not found; build a CMake preset before running.")
    if not arguments.no_open:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping Simulation Studio.")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
