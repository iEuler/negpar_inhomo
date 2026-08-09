from __future__ import annotations

import json
from pathlib import Path
import tempfile
import threading
import unittest
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from ui.server import (
    RequestError,
    ResultReader,
    create_server,
    read_numeric_file,
    resolve_output_directory,
    validate_run_request,
)


class RunRequestTests(unittest.TestCase):
    def test_run_request_maps_supported_cli_options(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = root / "negpar_inhomo.exe"
            executable.touch()

            request = validate_run_request({
                "seedMode": "fixed",
                "seed": 20260809,
                "threads": 4,
                "steps": 100,
                "outputDirectory": "result/run-a",
            }, root, executable)

            self.assertEqual(request.seed, 20260809)
            self.assertEqual(request.threads, 4)
            self.assertEqual(request.steps, 100)
            self.assertEqual(request.output_directory, root / "result/run-a")
            self.assertEqual(request.command[1:7], (
                "--seed", "20260809", "--threads", "4", "--steps", "100"))
            self.assertEqual(request.command[-2], "--output-dir")

    def test_random_seed_is_omitted_from_command(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = root / "negpar_inhomo"
            executable.touch()
            request = validate_run_request({
                "seedMode": "random",
                "threads": 1,
                "steps": 1,
                "outputDirectory": "result/random",
            }, root, executable)
            self.assertIsNone(request.seed)
            self.assertNotIn("--seed", request.command)

    def test_output_directory_cannot_escape_repository(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaises(RequestError):
                resolve_output_directory(root, "../outside")

    def test_nonempty_output_directory_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = root / "negpar_inhomo"
            executable.touch()
            output = root / "result/existing"
            output.mkdir(parents=True)
            (output / "run_metadata.txt").write_text("seed 1\n", encoding="utf-8")
            with self.assertRaisesRegex(RequestError, "not empty"):
                validate_run_request({
                    "seedMode": "fixed", "seed": 1, "threads": 1, "steps": 1,
                    "outputDirectory": "result/existing",
                }, root, executable)

    def test_output_directory_cannot_be_an_existing_file(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = root / "negpar_inhomo"
            executable.touch()
            (root / "result").mkdir()
            (root / "result/not-a-directory").write_text("occupied", encoding="utf-8")
            with self.assertRaisesRegex(RequestError, "existing file"):
                validate_run_request({
                    "seedMode": "fixed", "seed": 1, "threads": 1, "steps": 1,
                    "outputDirectory": "result/not-a-directory",
                }, root, executable)

    def test_invalid_numeric_inputs_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            executable = root / "negpar_inhomo"
            executable.touch()
            with self.assertRaisesRegex(RequestError, "Threads"):
                validate_run_request({
                    "seedMode": "fixed", "seed": 1, "threads": 0, "steps": 1,
                    "outputDirectory": "result/invalid",
                }, root, executable)


class ResultReaderTests(unittest.TestCase):
    def _write(self, directory: Path, name: str, values: list[float]) -> None:
        (directory / name).write_text(
            "\n".join(str(value) for value in values) + "\n", encoding="utf-8")

    def test_numeric_reader_rejects_nonfinite_values(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "bad.txt"
            path.write_text("1\nnan\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "NaN or Inf"):
                read_numeric_file(path)

    def test_reader_discovers_snapshots_and_diagnostics(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            self._write(output, "x.txt", [0.0, 0.5, 1.0])
            for name, values in {
                "rho_000.txt": [1.0, 1.1, 1.0],
                "rhoF_000.txt": [0.9, 1.0, 0.9],
                "rhoM_000.txt": [1.0, 1.0, 1.0],
                "rho.txt": [1.0, 1.2, 1.0],
                "rhoF.txt": [0.9, 1.1, 0.9],
                "rhoM.txt": [1.0, 1.0, 1.0],
                "time_dist.txt": [0.0],
                "total_energy.txt": [5.0, 5.001],
                "Np_rec.txt": [100.0, 101.0],
                "Nn_rec.txt": [20.0, 21.0],
                "Nf_rec.txt": [1000.0, 1000.0],
            }.items():
                self._write(output, name, values)
            (output / "run_metadata.txt").write_text(
                "seed 123\nthreads 1\nbuild_type release\n", encoding="utf-8")

            payload = ResultReader(output).payload("density", 0)

            self.assertEqual(payload["snapshot"], 0)
            self.assertEqual(payload["snapshots"], [
                {"id": 0, "label": "t = 0.000"},
                {"id": "final", "label": "final"},
            ])
            self.assertEqual(payload["x"], [0.0, 0.5, 1.0])
            self.assertEqual([series["name"] for series in payload["series"]],
                             ["rho", "rhoF", "rhoM"])
            self.assertEqual(len(payload["diagnostics"]["conservation"]), 1)
            self.assertEqual(len(payload["diagnostics"]["particles"]), 3)
            self.assertEqual(payload["metadata"]["seed"], "123")

    def test_bad_series_is_reported_without_breaking_other_results(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            self._write(output, "x.txt", [0.0, 1.0])
            self._write(output, "rho.txt", [1.0, 2.0])
            (output / "rhoF.txt").write_text("nan\n", encoding="utf-8")

            payload = ResultReader(output).payload("density", "final")

            self.assertEqual([series["name"] for series in payload["series"]], ["rho"])
            self.assertTrue(any("rhoF.txt" in warning for warning in payload["warnings"]))


class HttpServerTests(unittest.TestCase):
    def test_static_dashboard_and_api_are_served_locally(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            missing_executable = root / "missing-negpar"
            server = create_server("127.0.0.1", 0, str(missing_executable), root)
            thread = threading.Thread(target=server.serve_forever, daemon=True)
            thread.start()
            base_url = f"http://127.0.0.1:{server.server_port}"
            try:
                with urlopen(f"{base_url}/", timeout=5) as response:
                    page = response.read().decode("utf-8")
                self.assertIn("NegPar Simulation Studio", page)

                with urlopen(f"{base_url}/api/config", timeout=5) as response:
                    config = json.load(response)
                self.assertFalse(config["executableAvailable"])
                self.assertEqual(config["repository"], str(root))

                request = Request(
                    f"{base_url}/api/runs",
                    data=json.dumps({
                        "seedMode": "fixed", "seed": 1, "threads": 1,
                        "steps": 1, "outputDirectory": "result/http-test",
                    }).encode("utf-8"),
                    headers={"Content-Type": "application/json"},
                    method="POST",
                )
                with self.assertRaises(HTTPError) as context:
                    urlopen(request, timeout=5)
                self.assertEqual(context.exception.code, 400)
                error = json.loads(context.exception.read().decode("utf-8"))
                self.assertIn("executable", error["error"].lower())
            finally:
                server.shutdown()
                server.server_close()
                thread.join(timeout=5)


if __name__ == "__main__":
    unittest.main()
