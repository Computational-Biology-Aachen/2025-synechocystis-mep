"""Pigment contents from the proteomics spectra (fit targets, ADR 0006 §2.8)."""

from synmep import pigments


def test_extract_matches_the_models_chl_per_od():
    ext = pigments.methanol_extract()
    # the models assume 7.5 mg Chl OD^-1 L^-1
    assert ext["chl_ug_per_ml_od"].between(6.0, 8.5).all()


def test_contents_are_physiological():
    t = pigments.pigment_targets()
    assert t["chl_mM"].between(12, 25).all()
    car = t.loc["sqs", "car_mM"]
    assert car.between(2, 10).all()
    assert car.is_monotonic_decreasing  # carotenoids per cell fall over the batch


def test_in_vivo_chl_is_calibrated_on_the_extracts():
    c = pigments.pigment_contents()
    sqs = c[c["strain"] == "sqs"]
    ratio = (sqs["chl_mM"] / sqs["chl_extract_mM"]).median()
    assert abs(ratio - 1) < 0.02
