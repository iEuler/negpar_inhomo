# NegPar Simulation Studio Design System

## Product context

NegPar Simulation Studio is a local, browser-based control surface for the
`negpar_inhomo` kinetic plasma simulation. It is used by computational physics
researchers who need to configure reproducible runs, watch execution progress,
inspect numerical outputs, and compare conservation and distribution behavior
without assembling command lines or plotting text files manually.

The first release wraps the existing executable and its supported run options:
seed, OpenMP thread count, step count, and isolated output directory. Advanced
physical parameters may be exposed later, but the UI must clearly distinguish
currently editable run controls from read-only compiled configuration.

## Jobs to be done

1. Configure a deterministic or randomly seeded simulation run.
2. Validate inputs before launching an expensive computation.
3. Start, monitor, and safely stop one local run.
4. See command, seed, status, elapsed time, step progress, and console output.
5. Inspect density, velocity, temperature, electric field, energy, particle
   count, resampling, and runtime series produced in the output directory.
6. Switch between spatial snapshots and time-series diagnostics.
7. Locate the output directory and retain enough metadata to reproduce a run.

## Information architecture

Use one desktop-first workspace with three persistent regions:

- Top command bar: product identity, executable/build status, run status,
  output-folder action, and primary Run/Stop action.
- Left configuration rail (320-360px): run preset, seed policy/value, threads,
  steps, output directory, effective CLI preview, validation, and launch action.
- Main analysis canvas: run overview metrics followed by a large primary plot,
  plot selector, series legend, snapshot/time controls, and secondary panels for
  conservation, particle counts, and console output.

The empty state should remain useful: show a restrained example plot scaffold
and explain which artifacts will populate each panel. During a run, pin status
and progress without covering results. After completion, foreground plots and
reproducibility metadata.

## Real data contract

Run options map exactly to:

`negpar_inhomo --seed <uint32> --threads <positive int> --steps <positive int> --output-dir <path>`

The result viewer should recognize these artifact families:

- Spatial macro snapshots: `rho`, `u1`, `u2`, `u3`, `Tprt`, `rhoF`, `u1F`,
  `TprtF`, `rhoM`, `u1M`, `TprtM`, `elecfield`, and `elecfield_F`.
- Conservation/time series: `elec_energy`, `elec_energy_F`, `total_energy`,
  `total_energy_F`, `m0`, `m1`, and `m2`.
- Particle/time series: `Np_rec`, `Nn_rec`, `Nf_rec`, `Neff_F_rec`, and
  `num_resample`.
- Distribution outputs: `dist`, particle coordinates, radial distributions,
  and velocity histograms.
- Provenance: `run_metadata.txt`, `parameter.txt`, and `parameter2.txt`.

Do not imply live streaming of numerical files unless the implementation can
read a fully written snapshot. Console output may stream live; plots should
update only from stable files.

## Visual direction

Style source: Mosaic Grid Architecture, adapted from a marketing-page system
into a dense technical application. The mood is technical minimalist: precise,
calm, credible, and instrument-like rather than futuristic or decorative.

### Palette

- Paper background: `#F7F7F5`
- Surface: `#FFFFFF`
- Primary forest: `#1A3C2B`
- Ink: `#202321`
- Secondary ink: `#626862`
- Hairline grid: `rgba(58, 58, 56, 0.20)`
- Selected/healthy mint: `#9EFFBF`
- Warning gold: `#F4D35E`
- Error/stop coral: `#FF8C69`
- Plot series should use forest, coral, gold, teal, and slate with line-style
  differences so meaning does not depend on color alone.

### Typography

- Space Grotesk for product title and section headings; tight tracking.
- General Sans or system sans-serif for body copy and form values.
- JetBrains Mono for labels, numeric metrics, paths, commands, logs, axes, and
  provenance metadata.
- Application headings are compact (20-32px), not landing-page hero scale.
- Tabular numeric alignment is mandatory for metrics and axis values.

### Geometry and spacing

- 1px dividers and explicit grid structure; no shadows and no gradients.
- Square corners by default; 2px radius maximum.
- 8px base spacing unit, with dense 8/12/16px control spacing and 24/32px panel
  padding.
- Use flat white panels on Paper; selected panels may use a very pale mint fill.
- Charts use a clean plotting field, subtle gridlines, 1.5-2px series strokes,
  visible focus points, and compact monospaced legends.

## Components

- Status badge: square 8px indicator, hairline border, uppercase mono label.
- Run button: solid forest; Stop button: coral outline or fill with explicit
  confirmation only when a run is active.
- Inputs: label above value, unit/help beneath when necessary, visible invalid
  border and inline reason, no placeholder-only labels.
- Metric tile: mono eyebrow, large tabular value, unit, optional delta/sparkline.
- Plot toolbar: segmented metric selector, snapshot slider/dropdown, scale and
  download actions.
- Console: dark forest/ink background is allowed as the single inverted panel,
  with compact mono output and clear auto-scroll state.
- Corner markers may be used on the configuration rail and primary plot to
  reinforce technical precision, but never as ornamental clutter.

## Interaction and motion

- Snappy 120-180ms ease-out transitions for hover, selection, and panel state.
- Linear progress motion; never use indefinite decorative animation while idle.
- Respect reduced-motion preferences.
- Running state disables configuration fields that affect the launched process.
- Destructive Stop is visually distinct and preserves completed output.
- Keyboard focus must be obvious; every chart control and form field is usable
  without a pointer.

## Accessibility and responsive behavior

- Meet WCAG AA contrast and never encode run status or series identity by color
  alone.
- Desktop is primary at 1440px. At tablet widths, configuration becomes a
  collapsible drawer. At narrow widths, stack configuration, metrics, plots,
  and console; retain the top Run/Stop action.
- Provide textual summaries for plots and expose raw values/downloads.

## Design constraints

- Use ONLY the fonts, colors, spacing, and component styles defined here.
- Do not introduce shadows, gradients, glass effects, neon colors, decorative
  serif fonts, excessive rounding, or generic analytics-dashboard chrome.
- Prefer exact simulation vocabulary over business terms such as revenue,
  conversion, users, or KPI.
- Show realistic values and filenames from the documented simulation contract.
