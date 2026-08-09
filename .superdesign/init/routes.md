# Routes

The app uses Python standard-library HTTP routing in ui/server.py and has no client router.

| URL | Method | Source | Layout |
| --- | --- | --- | --- |
| / | GET | ui/static/index.html | StudioShell |
| /styles.css | GET | ui/static/styles.css | n/a |
| /app.js | GET | ui/static/app.js | n/a |
| /api/config | GET | StudioHandler.do_GET | JSON |
| /api/status | GET | StudioHandler.do_GET | JSON |
| /api/results | GET | StudioHandler.do_GET | JSON |
| /api/saved-runs | GET | StudioHandler.do_GET | JSON |
| /api/compare | GET | StudioHandler.do_GET | JSON |
| /api/runs | POST | StudioHandler.do_POST | JSON |
| /api/runs/stop | POST | StudioHandler.do_POST | JSON |
| /api/open-output | POST | StudioHandler.do_POST | JSON |

The root page renders Run and Compare modes in one shell. Run provides configuration, progress, plots, console, and reproducibility output. Compare provides saved-run selectors, synchronized spatial plots, A-minus-B deltas, diagnostics, and provenance.
