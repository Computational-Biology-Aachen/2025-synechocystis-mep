"""Fit the PETC forms m0, m1, m2 and the degree-5 polynomial to pfennig2024 steady states
under white light and LED mixtures (ADR 0006 §2.4, checkpoint A).

    uv run python scripts/fit_petc_balance.py

Data: data/surrogates/petc_training.csv (white light) and petc_training_led.csv,
each split 80/20 (seed 0). All forms are fitted to the union of the training parts.
Hold-out R² is reported on: white, LED (all), and *imbalanced* = LED samples with one
LED colour >= 80 % of the light. Criterion for m1: R² >= ~0.97 on imbalanced.
Writes data/surrogates/petc_balance.json and docs/reports/petc-balance.json.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from synmep.petc import balance as b
from synmep.surrogates.poly import fit_polynomial, predict

DATA = Path("data/surrogates")
OUT = DATA / "petc_balance.json"
REPORT = Path("docs/reports/petc-balance.json")
LEDS = ("w_led_blue", "w_led_orange", "w_led_red")
OUTPUTS = ("J_3PGA", "J_ATP")


def load() -> tuple[pd.DataFrame, pd.DataFrame]:
    white = pd.read_csv(DATA / "petc_training.csv").dropna(subset=list(OUTPUTS))
    led = pd.read_csv(DATA / "petc_training_led.csv").dropna(subset=list(OUTPUTS))
    white["set"], led["set"] = "white", "led"
    led["imbalanced"] = led[list(LEDS)].max(axis=1) >= 0.8
    white["imbalanced"] = False
    rng = np.random.default_rng(0)
    parts = []
    for d in (white, led):
        d = d.reset_index(drop=True)
        test = np.zeros(len(d), bool)
        test[rng.permutation(len(d))[: len(d) // 5]] = True
        parts.append(d.assign(test=test))
    data = pd.concat(parts, ignore_index=True)
    return data[~data["test"]], data[data["test"]]


def r2(pred: np.ndarray, y: np.ndarray) -> float:
    return float(1 - np.sum((pred - y) ** 2) / np.sum((y - y.mean()) ** 2))


# parameters: name -> (start, transform); "log" positive, "logit" in (0, 1)
SPECS = {
    "m0": {"J_max": 20.0, "I_k": 2e4, "w_ps2": 10.0, "w_pbs": 15.0, "r": 0.3},
    "m1": {"J_max": 20.0, "I_k": 2e4, "w2": 3.0, "b": 15.0, "phi": 0.7,
           "k_ocp": 1e-3, "p": 4.0, "r": 0.3},
    "m2": {"w1": 1e-3, "w2": 3e-3, "b": 15.0, "phi": 0.7, "k_ocp": 1e-3,
           "k_b": 30.0, "V": 25.0, "K": 0.5, "r": 0.3},
}  # fmt: skip
LOGIT = {"phi"}


def unpack(names: list[str], z: np.ndarray) -> dict[str, float]:
    return {
        n: float(1 / (1 + np.exp(-v)) if n in LOGIT else np.exp(v))
        for n, v in zip(names, z, strict=True)
    }


def pack(p: dict[str, float]) -> np.ndarray:
    return np.array(
        [np.log(v / (1 - v)) if n in LOGIT else np.log(v) for n, v in p.items()]
    )


def j3_and_excess(form: str, p: dict[str, float], x: np.ndarray):
    if form == "m0":
        return b.m0(p, x), np.zeros(len(x))
    return getattr(b, form)(p, x)


def fit_form(form: str, train: pd.DataFrame) -> dict[str, float]:
    x, y = train[list(b.INPUTS)].to_numpy(), train["J_3PGA"].to_numpy()
    names = list(SPECS[form])
    starts = [pack(SPECS[form])]
    rng = np.random.default_rng(1)
    starts += [starts[0] + rng.normal(0, 0.7, len(names)) for _ in range(5)]
    best = min(
        (
            least_squares(
                lambda z: j3_and_excess(form, unpack(names, z), x)[0] - y,
                z0,
                max_nfev=400,
            )
            for z0 in starts
        ),
        key=lambda f: f.cost,
    )
    p = unpack(names, best.x)
    j3, ex = j3_and_excess(form, p, x)
    cols = [np.ones_like(j3), j3, j3**2] + ([ex] if form != "m0" else [])
    coef, *_ = np.linalg.lstsq(np.column_stack(cols), train["J_ATP"].to_numpy())
    p |= dict(zip(b.ATP, [*map(float, coef), 0.0][:4], strict=True))
    return p


def evaluate(pred: pd.DataFrame, test: pd.DataFrame) -> dict:
    out = {}
    for subset, mask in (
        ("white", test["set"] == "white"),
        ("led", test["set"] == "led"),
        ("imbalanced", test["imbalanced"]),
    ):
        out[subset] = {
            o: r2(pred.loc[mask, o].to_numpy(), test.loc[mask, o].to_numpy())
            for o in OUTPUTS
        }
    return out


def main() -> None:
    train, test = load()
    xt = test[list(b.INPUTS)].to_numpy()
    params, report = {}, {"n_train": len(train), "n_test": len(test),
                          "n_test_imbalanced": int(test["imbalanced"].sum())}  # fmt: skip
    for form in ("m0", "m1", "m2"):
        p = fit_form(form, train)
        j3, ex = j3_and_excess(form, p, xt)
        pred = pd.DataFrame(
            {"J_3PGA": j3, "J_ATP": b.j_atp(p, j3, ex)}, index=test.index
        )
        params[form] = p
        report[form] = evaluate(pred, test) | {"parameters": p}
        print(form, json.dumps(report[form]["imbalanced"]), flush=True)
    poly = fit_polynomial(train[list(b.INPUTS)], train[list(OUTPUTS)], 5)
    report["poly5"] = evaluate(
        predict(poly, test[list(b.INPUTS)]).set_axis(test.index), test
    )
    print("poly5", json.dumps(report["poly5"]["imbalanced"]), flush=True)
    OUT.write_text(json.dumps(params, indent=1))
    (DATA / "petc_poly_led.json").write_text(json.dumps(poly, indent=1))
    REPORT.write_text(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
