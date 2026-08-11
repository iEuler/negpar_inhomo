"use strict";

const state = {
  config: null,
  status: null,
  metric: "density",
  results: null,
  followLatestSnapshot: true,
  selectedSnapshot: null,
  seedMode: "fixed",
  autoScroll: true,
  refreshing: false,
  mode: "run",
  savedRuns: [],
  comparison: null,
  comparisonMetric: "density",
  comparisonSnapshot: null,
};

const colors = ["#1a3c2b", "#ff8c69", "#f4d35e", "#2b7a78", "#64748b"];
const svgNamespace = "http://www.w3.org/2000/svg";

const elements = {};

function element(id) {
  return document.getElementById(id);
}

async function api(path, options = {}) {
  const response = await fetch(path, {
    headers: { "Content-Type": "application/json", ...(options.headers || {}) },
    ...options,
  });
  let payload = {};
  try {
    payload = await response.json();
  } catch (_) {
    payload = { error: `Server returned ${response.status}.` };
  }
  if (!response.ok) {
    throw new Error(payload.error || `Request failed with ${response.status}.`);
  }
  return payload;
}

function cacheElements() {
  [
    "run-form", "preset-input", "seed-input", "threads-input", "steps-input",
    "output-input", "seed-help", "fixed-seed-button", "random-seed-button",
    "command-preview", "form-error", "top-run-button", "sidebar-run-button",
    "output-folder-button", "engine-dot", "engine-label", "build-label",
    "top-status-badge", "top-status-label", "status-value", "step-value",
    "steps-value", "elapsed-value", "file-count-value", "snapshot-slider",
    "snapshot-label", "macro-plot", "macro-summary", "conservation-plot",
    "particle-plot", "console-output", "autoscroll-button", "footer-seed",
    "footer-threads", "footer-build", "footer-path",
    "run-mode-button", "compare-mode-button", "run-command-actions",
    "run-workspace", "run-footer", "comparison-workspace",
    "refresh-runs-button", "run-a-select", "run-b-select",
    "run-a-metadata", "run-b-metadata", "swap-runs-button",
    "comparison-error", "comparison-snapshot-slider",
    "comparison-snapshot-label", "max-delta-value", "relative-l2-value",
    "drift-delta-value", "runtime-ratio-value", "run-a-plot-title",
    "run-b-plot-title", "run-a-plot", "run-b-plot", "delta-plot",
    "comparison-conservation-plot", "comparison-particle-plot",
    "provenance-body",
  ].forEach((id) => { elements[id] = element(id); });
}

function bindEvents() {
  elements["run-form"].addEventListener("submit", (event) => {
    event.preventDefault();
    runOrStop();
  });
  elements["top-run-button"].addEventListener("click", runOrStop);
  elements["fixed-seed-button"].addEventListener("click", () => setSeedMode("fixed"));
  elements["random-seed-button"].addEventListener("click", () => setSeedMode("random"));
  elements["preset-input"].addEventListener("change", applyPreset);
  ["seed-input", "threads-input", "steps-input", "output-input"].forEach((id) => {
    elements[id].addEventListener("input", () => {
      clearError();
      validateForm(false);
      updateCommandPreview();
    });
  });
  document.querySelectorAll(".metric-tab").forEach((button) => {
    button.addEventListener("click", () => {
      state.metric = button.dataset.metric;
      state.followLatestSnapshot = true;
      state.selectedSnapshot = null;
      document.querySelectorAll(".metric-tab").forEach((candidate) => {
        candidate.classList.toggle("is-active", candidate === button);
      });
      refreshResults();
    });
  });
  elements["snapshot-slider"].addEventListener("input", () => {
    const snapshots = state.results?.snapshots || [];
    const selected = snapshots[Number(elements["snapshot-slider"].value)];
    if (!selected) return;
    state.followLatestSnapshot = false;
    state.selectedSnapshot = selected.id;
    elements["snapshot-label"].textContent = selected.label;
    refreshResults();
  });
  elements["autoscroll-button"].addEventListener("click", () => {
    state.autoScroll = !state.autoScroll;
    elements["autoscroll-button"].setAttribute("aria-pressed", String(state.autoScroll));
    elements["autoscroll-button"].lastChild.textContent = state.autoScroll
      ? " Auto-scroll on" : " Auto-scroll off";
    if (state.autoScroll) scrollConsole();
  });
  elements["output-folder-button"].addEventListener("click", async () => {
    try {
      await api("/api/open-output", { method: "POST", body: "{}" });
    } catch (error) {
      showError(error.message);
    }
  });
  elements["run-mode-button"].addEventListener("click", () => setMode("run"));
  elements["compare-mode-button"].addEventListener("click", () => setMode("compare"));
  elements["refresh-runs-button"].addEventListener("click", loadSavedRuns);
  ["run-a-select", "run-b-select"].forEach((id) => {
    elements[id].addEventListener("change", refreshComparison);
  });
  elements["swap-runs-button"].addEventListener("click", () => {
    const runA = elements["run-a-select"].value;
    elements["run-a-select"].value = elements["run-b-select"].value;
    elements["run-b-select"].value = runA;
    refreshComparison();
  });
  document.querySelectorAll("[data-compare-metric]").forEach((button) => {
    button.addEventListener("click", () => {
      state.comparisonMetric = button.dataset.compareMetric;
      state.comparisonSnapshot = null;
      document.querySelectorAll("[data-compare-metric]").forEach((candidate) => {
        candidate.classList.toggle("is-active", candidate === button);
      });
      refreshComparison();
    });
  });
  elements["comparison-snapshot-slider"].addEventListener("input", () => {
    const snapshots = state.comparison?.snapshots || [];
    const selected = snapshots[Number(elements["comparison-snapshot-slider"].value)];
    if (!selected) return;
    state.comparisonSnapshot = selected.id;
    elements["comparison-snapshot-label"].textContent = selected.label;
    refreshComparison();
  });
}

function setSeedMode(mode) {
  state.seedMode = mode;
  const fixed = mode === "fixed";
  elements["fixed-seed-button"].classList.toggle("is-active", fixed);
  elements["random-seed-button"].classList.toggle("is-active", !fixed);
  elements["fixed-seed-button"].setAttribute("aria-pressed", String(fixed));
  elements["random-seed-button"].setAttribute("aria-pressed", String(!fixed));
  elements["seed-input"].disabled = !fixed || isRunActive();
  elements["seed-help"].textContent = fixed ? "Valid uint32" : "Generated by runtime";
  validateForm(false);
  updateCommandPreview();
}

function applyPreset() {
  const preset = elements["preset-input"].value;
  const defaults = state.config?.defaults || {};
  if (preset === "reference") {
    setSeedMode("fixed");
    elements["seed-input"].value = "123";
    elements["threads-input"].value = "1";
    elements["steps-input"].value = "1";
  } else if (preset === "standard") {
    setSeedMode(defaults.seedMode || "fixed");
    elements["seed-input"].value = String(defaults.seed ?? 20260809);
    elements["threads-input"].value = String(defaults.threads ?? 1);
    elements["steps-input"].value = String(defaults.steps ?? 100);
  }
  if (preset !== "custom") {
    elements["output-input"].value = freshOutputDirectory();
  }
  validateForm(false);
  updateCommandPreview();
}

function freshOutputDirectory() {
  const now = new Date();
  const pad = (value, width = 2) => String(value).padStart(width, "0");
  return `result/ui-run-${now.getFullYear()}${pad(now.getMonth() + 1)}${pad(now.getDate())}-${pad(now.getHours())}${pad(now.getMinutes())}${pad(now.getSeconds())}-${pad(now.getMilliseconds(), 3)}`;
}

function formPayload() {
  return {
    seedMode: state.seedMode,
    seed: Number(elements["seed-input"].value),
    threads: Number(elements["threads-input"].value),
    steps: Number(elements["steps-input"].value),
    outputDirectory: elements["output-input"].value.trim(),
  };
}

function validateForm(showMessages = true) {
  const payload = formPayload();
  const errors = [];
  const seedValid = state.seedMode === "random" || (Number.isInteger(payload.seed)
    && payload.seed >= 0 && payload.seed <= 4294967295);
  elements["seed-input"].classList.toggle("is-invalid", !seedValid);
  elements["seed-help"].classList.toggle("is-invalid", !seedValid);
  if (!seedValid) errors.push("Seed must be an integer from 0 through 4294967295.");
  if (!Number.isInteger(payload.threads) || payload.threads < 1) {
    errors.push("Threads must be a positive integer.");
  }
  if (!Number.isInteger(payload.steps) || payload.steps < 1) {
    errors.push("Steps must be a positive integer.");
  }
  if (!payload.outputDirectory) errors.push("Output directory is required.");
  if (showMessages && errors.length) showError(errors[0]);
  return errors.length === 0;
}

function quoteArgument(value) {
  return /\s/.test(value) ? `"${value.replaceAll('"', '\\"')}"` : value;
}

function updateCommandPreview() {
  const payload = formPayload();
  const pieces = ["negpar_inhomo"];
  if (state.seedMode === "fixed") pieces.push("--seed", String(payload.seed));
  pieces.push("--threads", String(payload.threads), "--steps", String(payload.steps),
    "--output-dir", quoteArgument(payload.outputDirectory || "<path>"));
  elements["command-preview"].textContent = pieces.join(" ");
}

function isRunActive() {
  return ["starting", "running", "stopping"].includes(state.status?.status);
}

async function runOrStop() {
  if (isRunActive()) {
    try {
      clearError();
      state.status = await api("/api/runs/stop", { method: "POST", body: "{}" });
      renderStatus();
    } catch (error) {
      showError(error.message);
    }
    return;
  }
  if (!validateForm(true)) return;
  if (!state.config?.executableAvailable) {
    showError("Build the Debug or Release CMake preset before starting a simulation.");
    return;
  }
  if (["completed", "failed", "stopped"].includes(state.status?.status)
      && elements["output-input"].value.trim() === state.status.outputDirectory) {
    elements["output-input"].value = freshOutputDirectory();
    updateCommandPreview();
  }
  try {
    clearError();
    state.followLatestSnapshot = true;
    state.selectedSnapshot = null;
    state.results = null;
    renderAllPlots();
    state.status = await api("/api/runs", {
      method: "POST",
      body: JSON.stringify(formPayload()),
    });
    renderStatus();
  } catch (error) {
    showError(error.message);
  }
}

function showError(message) {
  elements["form-error"].textContent = message;
  elements["form-error"].hidden = false;
}

function clearError() {
  elements["form-error"].hidden = true;
  elements["form-error"].textContent = "";
}

function setMode(mode) {
  state.mode = mode;
  const comparing = mode === "compare";
  elements["run-mode-button"].classList.toggle("is-active", !comparing);
  elements["compare-mode-button"].classList.toggle("is-active", comparing);
  elements["run-mode-button"].setAttribute("aria-pressed", String(!comparing));
  elements["compare-mode-button"].setAttribute("aria-pressed", String(comparing));
  elements["run-command-actions"].hidden = comparing;
  elements["run-workspace"].hidden = comparing;
  elements["run-footer"].hidden = comparing;
  elements["comparison-workspace"].hidden = !comparing;
  document.querySelector(".studio-shell").classList.toggle("is-comparing", comparing);
  if (comparing) loadSavedRuns();
}

function runByPath(path) {
  return state.savedRuns.find((run) => run.path === path);
}

function populateRunSelect(select, selectedPath) {
  select.replaceChildren();
  state.savedRuns.forEach((run) => {
    const option = document.createElement("option");
    option.value = run.path;
    option.textContent = run.name + " · " + run.status;
    select.append(option);
  });
  if (selectedPath && state.savedRuns.some((run) => run.path === selectedPath)) {
    select.value = selectedPath;
  }
}

function renderRunMetadata(container, run) {
  container.replaceChildren();
  if (!run) return;
  [
    ["Status", run.status],
    ["Seed", run.seed],
    ["Threads / steps", run.threads + " / " + run.steps],
    ["Build", run.build],
    ["Files", run.fileCount],
    ["Runtime", run.elapsedSeconds == null ? "-" : formatElapsed(run.elapsedSeconds)],
  ].forEach(([label, value]) => {
    const row = document.createElement("div");
    const term = document.createElement("dt");
    const detail = document.createElement("dd");
    term.textContent = label;
    detail.textContent = String(value);
    row.append(term, detail);
    container.append(row);
  });
}

function showComparisonError(message) {
  elements["comparison-error"].textContent = message;
  elements["comparison-error"].hidden = !message;
}

async function loadSavedRuns() {
  const selectedA = elements["run-a-select"].value;
  const selectedB = elements["run-b-select"].value;
  try {
    const payload = await api("/api/saved-runs");
    state.savedRuns = payload.runs || [];
    populateRunSelect(elements["run-a-select"], selectedA || state.savedRuns[0]?.path);
    populateRunSelect(elements["run-b-select"], selectedB || state.savedRuns[1]?.path
      || state.savedRuns[0]?.path);
    if (state.savedRuns.length) {
      await refreshComparison();
    } else {
      showComparisonError("No saved runs found under result/. Complete a run to compare it.");
      renderComparison();
    }
  } catch (error) {
    showComparisonError("Unable to discover saved runs: " + error.message);
  }
}

async function refreshComparison() {
  const runA = elements["run-a-select"].value;
  const runB = elements["run-b-select"].value;
  renderRunMetadata(elements["run-a-metadata"], runByPath(runA));
  renderRunMetadata(elements["run-b-metadata"], runByPath(runB));
  if (!runA || !runB) return;
  const parameters = new URLSearchParams({
    runA, runB, metric: state.comparisonMetric,
  });
  if (state.comparisonSnapshot !== null) {
    parameters.set("snapshot", String(state.comparisonSnapshot));
  }
  try {
    showComparisonError("");
    state.comparison = await api("/api/compare?" + parameters.toString());
    state.comparisonSnapshot = state.comparison.snapshot;
    renderComparisonSnapshotControl();
    renderComparison();
  } catch (error) {
    showComparisonError("Unable to compare runs: " + error.message);
  }
}

function renderComparisonSnapshotControl() {
  const snapshots = state.comparison?.snapshots || [];
  const slider = elements["comparison-snapshot-slider"];
  slider.disabled = snapshots.length < 2;
  slider.min = "0";
  slider.max = String(Math.max(0, snapshots.length - 1));
  let index = snapshots.findIndex((item) => item.id === state.comparison?.snapshot);
  if (index < 0) index = Math.max(0, snapshots.length - 1);
  slider.value = String(index);
  elements["comparison-snapshot-label"].textContent =
    snapshots[index]?.label || "No common snapshots";
}

function formatComparisonValue(value, suffix = "") {
  if (value == null || !Number.isFinite(value)) return "—";
  const magnitude = Math.abs(value);
  const formatted = (magnitude > 0 && magnitude < 0.001) || magnitude >= 1000
    ? value.toExponential(2) : value.toPrecision(3);
  return formatted + suffix;
}

function sharedDomain(...groups) {
  const values = groups.flatMap((series) =>
    (series || []).flatMap((item) => item.values || [])).filter(Number.isFinite);
  return values.length ? [Math.min(...values), Math.max(...values)] : null;
}

function pairedDiagnostics(runA, runB) {
  const result = [];
  const byNameB = new Map((runB || []).map((item) => [item.name, item.values]));
  (runA || []).forEach((item) => {
    result.push({ name: "A " + item.name, values: item.values });
    if (byNameB.has(item.name)) {
      result.push({ name: "B " + item.name, values: byNameB.get(item.name) });
    }
  });
  return result;
}

function renderProvenance(runA, runB) {
  const body = elements["provenance-body"];
  body.replaceChildren();
  [
    ["Output path", runA?.path, runB?.path],
    ["Seed", runA?.seed, runB?.seed],
    ["Threads", runA?.threads, runB?.threads],
    ["Steps", runA?.steps, runB?.steps],
    ["Build", runA?.build, runB?.build],
    ["Compiler", runA?.compiler, runB?.compiler],
    ["Status", runA?.status, runB?.status],
  ].forEach((values) => {
    const row = document.createElement("tr");
    values.forEach((value) => {
      const cell = document.createElement("td");
      cell.textContent = value == null ? "—" : String(value);
      row.append(cell);
    });
    body.append(row);
  });
}

function renderComparison() {
  const comparison = state.comparison;
  const summary = comparison?.summary || {};
  elements["max-delta-value"].textContent = formatComparisonValue(summary.maxAbsDelta);
  elements["relative-l2-value"].textContent = formatComparisonValue(summary.relativeL2);
  elements["drift-delta-value"].textContent =
    formatComparisonValue(summary.energyDriftDelta);
  elements["runtime-ratio-value"].textContent =
    formatComparisonValue(summary.runtimeRatio, summary.runtimeRatio == null ? "" : "×");
  const runA = comparison?.runs?.a;
  const runB = comparison?.runs?.b;
  renderRunMetadata(elements["run-a-metadata"], runA || runByPath(elements["run-a-select"].value));
  renderRunMetadata(elements["run-b-metadata"], runB || runByPath(elements["run-b-select"].value));
  elements["run-a-plot-title"].textContent = runA?.name || "Run A";
  elements["run-b-plot-title"].textContent = runB?.name || "Run B";

  const seriesA = comparison?.a?.series || [];
  const seriesB = comparison?.b?.series || [];
  const domain = sharedDomain(seriesA, seriesB);
  renderLineChart(elements["run-a-plot"], comparison?.a?.x || [], seriesA, false,
    { yDomain: domain });
  renderLineChart(elements["run-b-plot"], comparison?.b?.x || [], seriesB, false,
    { yDomain: domain });
  renderLineChart(elements["delta-plot"], comparison?.delta?.x || [],
    comparison?.delta?.series || [], true, { zeroBaseline: true });

  const conservation = pairedDiagnostics(
    comparison?.a?.diagnostics?.conservation,
    comparison?.b?.diagnostics?.conservation);
  const particles = pairedDiagnostics(
    comparison?.a?.diagnostics?.particles,
    comparison?.b?.diagnostics?.particles);
  renderLineChart(elements["comparison-conservation-plot"], [], conservation, true,
    { paired: true });
  renderLineChart(elements["comparison-particle-plot"], [], particles, true,
    { paired: true });
  renderProvenance(runA, runB);
}

function statusLabel(status) {
  return ({
    ready: "Ready", starting: "Starting", running: "Running", stopping: "Stopping",
    completed: "Completed", failed: "Failed", stopped: "Stopped",
  })[status] || "Ready";
}

function formatElapsed(seconds) {
  const total = Math.max(0, Math.floor(seconds || 0));
  const hours = Math.floor(total / 3600);
  const minutes = Math.floor((total % 3600) / 60);
  const remaining = total % 60;
  const pair = (value) => String(value).padStart(2, "0");
  return hours ? `${pair(hours)}:${pair(minutes)}:${pair(remaining)}`
    : `${pair(minutes)}:${pair(remaining)}`;
}

function renderStatus() {
  const current = state.status || {
    status: "ready", step: 0, steps: Number(elements["steps-input"].value) || 0,
    elapsedSeconds: 0, outputFiles: 0, logs: [],
  };
  const label = statusLabel(current.status);
  const active = ["starting", "running", "stopping"].includes(current.status);
  const failed = current.status === "failed";
  const warning = ["starting", "stopping", "stopped"].includes(current.status);

  elements["status-value"].textContent = label;
  elements["top-status-label"].textContent = label;
  elements["step-value"].textContent = String(current.step || 0);
  elements["steps-value"].textContent = String(current.steps || Number(elements["steps-input"].value) || 0);
  elements["elapsed-value"].textContent = formatElapsed(current.elapsedSeconds);
  elements["file-count-value"].textContent = String(current.outputFiles || 0);
  elements["output-folder-button"].disabled = !current.outputDirectoryAbsolute;

  const badgeDot = elements["top-status-badge"].querySelector(".status-square");
  [badgeDot, elements["engine-dot"]].forEach((dot) => {
    if (!dot) return;
    dot.classList.toggle("is-error", failed);
    dot.classList.toggle("is-warning", warning);
  });

  document.querySelectorAll(".run-button").forEach((button) => {
    button.classList.toggle("button-stop", active);
    button.disabled = current.status === "stopping" || !state.config?.executableAvailable;
    const text = active ? (current.status === "stopping" ? "Stopping" : "Stop simulation")
      : "Run simulation";
    const target = button.querySelector("span:last-child") || button;
    target.textContent = text;
  });

  ["preset-input", "threads-input", "steps-input", "output-input",
    "fixed-seed-button", "random-seed-button"].forEach((id) => {
    elements[id].disabled = active;
  });
  elements["seed-input"].disabled = active || state.seedMode === "random";

  elements["footer-seed"].textContent = current.seed ?? (state.seedMode === "random" ? "random" : elements["seed-input"].value);
  elements["footer-threads"].textContent = current.threads || elements["threads-input"].value || "—";
  elements["footer-build"].textContent = state.config?.build || "—";
  elements["footer-path"].textContent = current.outputDirectoryAbsolute || elements["output-input"].value || "—";
  renderConsole(current.logs || [], current.error);
}

function renderConsole(logs, error) {
  const container = elements["console-output"];
  const previousCount = container.childElementCount;
  container.replaceChildren();
  const visibleLogs = logs.length ? logs : [{ time: "--:--:--", text: "Ready to launch simulation." }];
  visibleLogs.forEach((entry) => {
    const row = document.createElement("div");
    row.className = "console-line";
    const timestamp = document.createElement("time");
    timestamp.textContent = entry.time;
    const text = document.createElement("span");
    text.textContent = entry.text;
    row.append(timestamp, text);
    container.append(row);
  });
  if (error && !visibleLogs.some((entry) => entry.text.includes(error))) {
    const row = document.createElement("div");
    row.className = "console-line is-error";
    row.innerHTML = `<time>${new Date().toLocaleTimeString([], { hour12: false })}</time>`;
    const text = document.createElement("span");
    text.textContent = error;
    row.append(text);
    container.append(row);
  }
  if (state.autoScroll && (previousCount !== container.childElementCount || logs.length)) {
    scrollConsole();
  }
}

function scrollConsole() {
  const container = elements["console-output"];
  container.scrollTop = container.scrollHeight;
}

async function refreshStatus() {
  if (state.refreshing) return;
  state.refreshing = true;
  try {
    state.status = await api("/api/status");
    renderStatus();
    if (state.status.outputFiles > 0) await refreshResults();
  } catch (error) {
    showError(`Dashboard connection failed: ${error.message}`);
  } finally {
    state.refreshing = false;
  }
}

async function refreshResults() {
  const parameters = new URLSearchParams({ metric: state.metric });
  if (!state.followLatestSnapshot && state.selectedSnapshot !== null) {
    parameters.set("snapshot", String(state.selectedSnapshot));
  }
  try {
    state.results = await api(`/api/results?${parameters}`);
    state.selectedSnapshot = state.results.snapshot;
    renderSnapshotControl();
    renderAllPlots();
  } catch (error) {
    showError(`Unable to read results: ${error.message}`);
  }
}

function renderSnapshotControl() {
  const snapshots = state.results?.snapshots || [];
  const slider = elements["snapshot-slider"];
  slider.disabled = snapshots.length < 2;
  slider.min = "0";
  slider.max = String(Math.max(0, snapshots.length - 1));
  let index = snapshots.findIndex((item) => item.id === state.results?.snapshot);
  if (index < 0) index = Math.max(0, snapshots.length - 1);
  slider.value = String(index);
  elements["snapshot-label"].textContent = snapshots[index]?.label || "No snapshots";
}

function renderAllPlots() {
  const results = state.results;
  renderLineChart(elements["macro-plot"], results?.x || [], results?.series || [], false);
  renderLineChart(elements["conservation-plot"], [], results?.diagnostics?.conservation || [], true);
  renderLineChart(elements["particle-plot"], [], results?.diagnostics?.particles || [], true);
  const series = results?.series || [];
  elements["macro-summary"].textContent = series.length
    ? `${series.map((item) => item.name).join(" · ")} · ${results.x.length} spatial samples`
    : "Run a simulation to populate spatial results.";
}

function svgElement(name, attributes = {}) {
  const node = document.createElementNS(svgNamespace, name);
  Object.entries(attributes).forEach(([key, value]) => node.setAttribute(key, String(value)));
  return node;
}

function renderLineChart(container, xValues, rawSeries, compact, options = {}) {
  container.replaceChildren();
  const series = rawSeries.filter((item) => Array.isArray(item.values) && item.values.length);
  if (!series.length) {
    const empty = document.createElement("div");
    empty.className = "plot-empty";
    empty.textContent = "Awaiting stable simulation output";
    container.append(empty);
    return;
  }

  const width = 1000;
  const height = compact ? 210 : 340;
  const margin = compact
    ? { top: 30, right: 16, bottom: 24, left: 54 }
    : { top: 38, right: 24, bottom: 38, left: 68 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const longest = Math.max(...series.map((item) => item.values.length));
  const xs = xValues.length === longest ? xValues : Array.from({ length: longest }, (_, index) => index);
  const finiteValues = series.flatMap((item) => item.values).filter(Number.isFinite);
  let yMin = options.yDomain ? options.yDomain[0] : Math.min(...finiteValues);
  let yMax = options.yDomain ? options.yDomain[1] : Math.max(...finiteValues);
  if (options.zeroBaseline) {
    yMin = Math.min(0, yMin);
    yMax = Math.max(0, yMax);
  }
  if (yMax === yMin) {
    const padding = Math.max(1, Math.abs(yMax) * 0.05);
    yMin -= padding;
    yMax += padding;
  } else {
    const padding = (yMax - yMin) * 0.08;
    yMin -= padding;
    yMax += padding;
  }
  let xMin = Math.min(...xs);
  let xMax = Math.max(...xs);
  if (xMax === xMin) xMax = xMin + 1;

  const xScale = (value) => margin.left + ((value - xMin) / (xMax - xMin)) * plotWidth;
  const yScale = (value) => margin.top + ((yMax - value) / (yMax - yMin)) * plotHeight;
  const svg = svgElement("svg", { viewBox: `0 0 ${width} ${height}`, preserveAspectRatio: "none" });
  svg.setAttribute("aria-hidden", "true");

  const horizontalLines = compact ? 3 : 5;
  for (let index = 0; index <= horizontalLines; index += 1) {
    const ratio = index / horizontalLines;
    const y = margin.top + ratio * plotHeight;
    svg.append(svgElement("line", {
      x1: margin.left, x2: width - margin.right, y1: y, y2: y,
      stroke: "rgba(58,58,56,0.18)", "stroke-width": 1,
    }));
    const label = svgElement("text", {
      x: margin.left - 9, y: y + 4, "text-anchor": "end",
      fill: "#626862", "font-family": "Consolas, monospace", "font-size": compact ? 17 : 13,
    });
    label.textContent = formatAxis(yMax - ratio * (yMax - yMin));
    svg.append(label);
  }
  if (options.zeroBaseline && yMin <= 0 && yMax >= 0) {
    const zeroY = yScale(0);
    svg.append(svgElement("line", {
      x1: margin.left, x2: width - margin.right, y1: zeroY, y2: zeroY,
      stroke: "#202321", "stroke-width": 1.5, "stroke-dasharray": "5 4",
      "vector-effect": "non-scaling-stroke",
    }));
  }
  const verticalLines = compact ? 4 : 6;
  for (let index = 0; index <= verticalLines; index += 1) {
    const ratio = index / verticalLines;
    const x = margin.left + ratio * plotWidth;
    svg.append(svgElement("line", {
      x1: x, x2: x, y1: margin.top, y2: height - margin.bottom,
      stroke: "rgba(58,58,56,0.12)", "stroke-width": 1,
    }));
  }

  series.forEach((item, seriesIndex) => {
    const colorIndex = options.paired ? Math.floor(seriesIndex / 2) : seriesIndex;
    const path = item.values.map((value, index) => {
      const x = xs[Math.min(index, xs.length - 1)] ?? index;
      return `${index ? "L" : "M"}${xScale(x).toFixed(2)},${yScale(value).toFixed(2)}`;
    }).join(" ");
    const attributes = {
      d: path, fill: "none", stroke: colors[colorIndex % colors.length],
      "stroke-width": compact ? 2 : 2.25, "vector-effect": "non-scaling-stroke",
    };
    if ((options.paired && seriesIndex % 2 === 1) || (!options.paired && seriesIndex === 1)) {
      attributes["stroke-dasharray"] = "7 4";
    }
    if (!options.paired && seriesIndex === 2) attributes["stroke-dasharray"] = "2 3";
    svg.append(svgElement("path", attributes));

    const legendX = margin.left + seriesIndex * (compact ? 180 : 150);
    svg.append(svgElement("line", {
      x1: legendX, x2: legendX + 24, y1: 15, y2: 15,
      stroke: colors[colorIndex % colors.length], "stroke-width": 3,
    }));
    const legend = svgElement("text", {
      x: legendX + 31, y: 20, fill: "#202321",
      "font-family": "Consolas, monospace", "font-size": compact ? 17 : 13,
    });
    legend.textContent = item.name;
    svg.append(legend);
  });
  container.append(svg);
}

function formatAxis(value) {
  const magnitude = Math.abs(value);
  if ((magnitude > 0 && magnitude < 0.001) || magnitude >= 10000) return value.toExponential(1);
  if (magnitude >= 100) return value.toFixed(0);
  if (magnitude >= 1) return value.toFixed(2);
  return value.toFixed(3);
}

async function initialize() {
  cacheElements();
  bindEvents();
  try {
    state.config = await api("/api/config");
    const defaults = state.config.defaults;
    elements["seed-input"].value = String(defaults.seed);
    elements["threads-input"].value = String(defaults.threads);
    elements["steps-input"].value = String(defaults.steps);
    elements["output-input"].value = defaults.outputDirectory;
    setSeedMode(defaults.seedMode);
    elements["engine-label"].textContent = state.config.executableAvailable
      ? "NEGPAR_INHOMO: ACTIVE" : "EXECUTABLE MISSING";
    elements["engine-dot"].classList.toggle("is-error", !state.config.executableAvailable);
    elements["build-label"].textContent = `${state.config.build} · local`;
    elements["footer-build"].textContent = state.config.build;
    elements["top-run-button"].disabled = !state.config.executableAvailable;
    elements["sidebar-run-button"].disabled = !state.config.executableAvailable;
    updateCommandPreview();
    validateForm(false);
    await refreshStatus();
    await loadSavedRuns();
  } catch (error) {
    showError(`Unable to initialize dashboard: ${error.message}`);
  }
  window.setInterval(refreshStatus, 1000);
}

window.addEventListener("DOMContentLoaded", initialize);
