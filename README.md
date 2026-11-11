# Synechocystis MEP model

Kinetic model of terpenoid (squalene) production via the MEP pathway in
*Synechocystis* sp. PCC 6803, coupled to a surrogate of the pfennig2024
photosynthetic electron transport model and a growth model. Goals: maximal squalene
productivity of an MCA-optimised strain in a batch culture, and the best LED light
composition (ADR 0006).

## Installation

Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then

```bash
uv sync --all-groups            # model + notebook
uv sync --all-groups --extra surrogate   # + mxlmodels, to retrain the PETC surrogate
```

## Layout

- `src/synmep/` — data loading, surrogates and the model families `single`, `double`,
  `minimal` (`models/`; families and names: ADR 0006)
- `main.ipynb` — figures
- `data/` — raw data; `data/reference/` holds outputs of the legacy modelbase model
- `docs/` — ADRs (`docs/adrs/CONTEXT.md` is the entry point), units, assumptions, reports
- `scripts/legacy/` — scripts that only ran in the former pixi/modelbase environment
