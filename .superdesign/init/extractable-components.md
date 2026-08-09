# Extractable component catalog

This application is one static document, so none of these regions is currently a standalone source component. Keep them inline for this design round.

## StudioShell
- Source: ui/static/index.html
- Category: layout
- Description: Command bar, Run/Compare mode switch, active workspace, and run footer.
- Extractable props: activeMode, status
- Hardcoded: NP mark, product title, CSS classes

## CommandBar
- Source: ui/static/index.html
- Category: layout
- Description: Identity, engine state, output action, status, and primary action.
- Extractable props: activeMode, engineStatus, runStatus
- Hardcoded: labels, inline symbols, CSS

## ConfigurationPanel
- Source: ui/static/index.html
- Category: layout
- Description: Run inputs, compiled configuration, and launch action.
- Extractable props: preset, seedMode, seed, threads, steps, outputDirectory
- Hardcoded: labels and CSS

## MetricTile
- Source: ui/static/index.html
- Category: basic
- Description: Run-summary card with label and monospaced value.
- Extractable props: label, value, suffix
- Hardcoded: article markup and CSS

## PlotPanel
- Source: ui/static/index.html
- Category: basic
- Description: Chart surface with header, controls, grid plot, and caption.
- Extractable props: title, activeMetric, snapshotLabel, summary
- Hardcoded: metric labels and CSS

## ConsolePanel
- Source: ui/static/index.html
- Category: basic
- Description: Forest console with timestamped lines.
- Extractable props: autoScroll, lines
- Hardcoded: title and CSS

## ComparisonPanel
- Source: ui/static/index.html
- Category: layout
- Description: Saved Run A/B selectors, metadata summaries, refresh, and swap controls.
- Extractable props: runA, runB, savedRuns
- Hardcoded: comparison labels and CSS

## ComparisonPlotGrid
- Source: ui/static/index.html
- Category: basic
- Description: Synchronized A/B charts with shared scale and A-minus-B delta.
- Extractable props: activeMetric, snapshot, runA, runB
- Hardcoded: metric labels and CSS
