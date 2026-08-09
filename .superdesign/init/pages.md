# Page dependency trees

## / (Simulation Studio)

Entry: ui/static/index.html

Dependencies:
- ui/static/styles.css
- ui/static/app.js
  - GET /api/config
  - GET /api/status
  - GET /api/results
  - GET /api/saved-runs
  - GET /api/compare
  - POST /api/runs
  - POST /api/runs/stop
  - POST /api/open-output
- ui/server.py
  - SimulationManager
  - ResultReader
  - compiled negpar_inhomo executable

There are no browser-side local imports or external framework, font, icon, or image dependencies.
