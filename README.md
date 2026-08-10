# Galaxy Clustering Correlation Lab

Two-point correlation functions, bias and BAO-scale structure.

Created and maintained by Biswajit Jana.

## Scientific Purpose

This zero-build browser laboratory puts a compact reference-data bundle in front of the simulation. The app loads `data/reference.json`, renders those published anchors first, then sends the adjustable model to `physicsWorker.js` so numerical work stays off the UI thread.

## Architecture

- `index.html`: mission-control interface.
- `styles.css`: dense dark scientific dashboard.
- `app.js`: UI state, Canvas rendering and worker orchestration.
- `physicsWorker.js`: numerical model and heatmap generation.
- `data/reference.json`: small auditable reference-data bundle.
- `scripts/validate.js`: no-dependency repository validation.

## Run

```bash
python -m http.server 8080
```

Open `http://localhost:8080`.

## Validate

```bash
npm run check
```

The validation script checks required files, JSON reference data, worker syntax, citations and absence of unfinished scaffold tokens.

## Reference Data

Representative xi(r) anchors showing power-law clustering and the baryon-acoustic feature.

## References

- Peebles, P.J.E., 1980. The large-scale structure of the universe. Princeton University Press.
- Eisenstein, D.J. et al., 2005. Detection of the baryon acoustic peak in the large-scale correlation function of SDSS luminous red galaxies. The Astrophysical Journal, 633(2), pp.560-574.
