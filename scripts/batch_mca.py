"""MCA under batch conditions for one variant (ADR 0006 §7, goal 1 follow-up).

    uv run --extra surrogate python scripts/batch_mca.py <variant>

Writes docs/reports/batch-mca/<variant>.json: the sequence from
`synmep.strains.optimise_batch` (warm white, Fig 9 protocol) and, per strain, the
batch summary next to the paper-protocol sequence from the goal runs.
"""

from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path

from synmep import evaluate, strains, variants

warnings.filterwarnings("ignore")
OUT = Path("docs/reports/batch-mca")


def main() -> None:
    name = sys.argv[1]
    a = variants.adapter(name)
    t0 = time.perf_counter()
    res = strains.optimise_batch(a)
    print(name, "batch sequence", res["sequence"], flush=True)
    per = []
    for s in strains.strains(res):
        b = evaluate.batch(a, s)
        summary = None if b is None else {
            k: b[k] for k in ("best_harvest_h", "best_productivity_mg_L_h",
                              "titer_end_mg_L")
        }  # fmt: skip
        per.append({"strain": list(s), "batch": summary})
        print(" ", list(s), summary, flush=True)
    report = {"variant": name, "batch_mca": res, "strains": per,
              "wall_s": time.perf_counter() - t0}  # fmt: skip
    (OUT / f"{name}.json").write_text(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
