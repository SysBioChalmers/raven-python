"""Tests for compare_fluxes (analysis/compare_fluxes.py, compareFluxes port)."""
import cobra
import numpy as np
import pandas as pd
import pytest

from raven_toolbox.analysis import compare_fluxes


@pytest.fixture
def model():
    """A <-> B <-> C, plus a side branch D touched only by r3.

    r1: A -> B, r2: B -> C, r3: B -> D. metabolite_list filters exercise
    "B" (touches all three) vs "D" (touches only r3).
    """
    m = cobra.Model("t")
    a, b, c, d = (cobra.Metabolite(x, name=x, compartment="c") for x in "ABCD")
    m.add_metabolites([a, b, c, d])
    r1 = cobra.Reaction("r1", name="first", lower_bound=-1000, upper_bound=1000)
    r1.add_metabolites({a: -1, b: 1})
    r2 = cobra.Reaction("r2", name="second", lower_bound=-1000, upper_bound=1000)
    r2.add_metabolites({b: -1, c: 1})
    r3 = cobra.Reaction("r3", name="third", lower_bound=-1000, upper_bound=1000)
    r3.add_metabolites({b: -1, d: 1})
    m.add_reactions([r1, r2, r3])
    return m


def test_identical_fluxes_change_nothing(model):
    fluxes = pd.Series({"r1": 10.0, "r2": 10.0, "r3": 10.0})
    result = compare_fluxes(model, fluxes, fluxes)
    assert result.changed.empty
    assert list(result.changed.columns) == [
        "reaction", "name", "flux1", "flux2", "abs_delta", "rel_change", "type",
    ]
    assert result.turned_on == []
    assert result.turned_off == []
    assert result.flipped == []


def test_changes_are_sorted_by_absolute_difference(model):
    a = pd.Series({"r1": 10.0, "r2": 1.0, "r3": 100.0})
    b = pd.Series({"r1": 11.0, "r2": 6.0, "r3": 97.0})
    result = compare_fluxes(model, a, b)
    assert list(result.changed["reaction"]) == ["r2", "r3", "r1"]
    assert list(result.changed["abs_delta"]) == [5.0, 3.0, 1.0]
    assert list(result.changed["name"]) == ["second", "third", "first"]


def test_ties_keep_model_order(model):
    a = pd.Series({"r1": 0.0, "r2": 0.0, "r3": 0.0})
    b = pd.Series({"r1": 5.0, "r2": 5.0, "r3": 5.0})
    result = compare_fluxes(model, a, b)
    assert list(result.changed["reaction"]) == ["r1", "r2", "r3"]


def test_turned_on_off_and_flipped_are_classified(model):
    a = pd.Series({"r1": 0.0, "r2": 4.0, "r3": 3.0})
    b = pd.Series({"r1": 7.0, "r2": 0.0, "r3": -3.0})
    result = compare_fluxes(model, a, b)
    assert result.turned_on == ["r1"]
    assert result.turned_off == ["r2"]
    assert result.flipped == ["r3"]
    types = dict(zip(result.changed["reaction"], result.changed["type"], strict=True))
    assert types == {"r1": "on", "r2": "off", "r3": "flip"}


def test_relative_change_is_nan_without_a_baseline(model):
    a = pd.Series({"r1": 0.0, "r2": 4.0, "r3": 0.0})
    b = pd.Series({"r1": 7.0, "r2": 6.0, "r3": 0.0})
    result = compare_fluxes(model, a, b)
    rel = dict(zip(result.changed["reaction"], result.changed["rel_change"], strict=True))
    assert np.isnan(rel["r1"])            # emerged from nothing
    assert rel["r2"] == pytest.approx(0.5)


def test_relative_change_is_negative_for_a_drop(model):
    a = pd.Series({"r1": 10.0, "r2": 0.0, "r3": 0.0})
    b = pd.Series({"r1": 4.0, "r2": 0.0, "r3": 0.0})
    result = compare_fluxes(model, a, b)
    assert result.changed.loc[0, "rel_change"] == pytest.approx(-0.6)


def test_cutoff_excludes_small_changes(model):
    a = pd.Series({"r1": 10.0, "r2": 1.0, "r3": 100.0})
    b = pd.Series({"r1": 10.4, "r2": 1.0, "r3": 102.0})
    result = compare_fluxes(model, a, b, cutoff=0.5)
    assert list(result.changed["reaction"]) == ["r3"]


def test_cutoff_governs_activity_not_just_difference(model):
    # r1 moves 0.2 -> 0.9 with cutoff 1.0: neither flux counts as active, so it
    # is neither turned on nor changed.
    a = pd.Series({"r1": 0.2, "r2": 0.0, "r3": 0.0})
    b = pd.Series({"r1": 0.9, "r2": 0.0, "r3": 0.0})
    result = compare_fluxes(model, a, b, cutoff=1.0)
    assert result.changed.empty
    assert result.turned_on == []


def test_metabolite_list_restricts_the_comparison(model):
    a = pd.Series({"r1": 10.0, "r2": 10.0, "r3": 10.0})
    b = pd.Series({"r1": 20.0, "r2": 20.0, "r3": 20.0})
    result = compare_fluxes(model, a, b, metabolite_list=["D"])
    # Only r3 touches D; the others are absent from the result, not reported
    # as unchanged.
    assert list(result.changed["reaction"]) == ["r3"]
    assert result.n_considered == 1


def test_metabolite_list_matches_names_case_insensitively(model):
    a = pd.Series({"r1": 10.0, "r2": 10.0, "r3": 10.0})
    b = pd.Series({"r1": 20.0, "r2": 10.0, "r3": 10.0})
    result = compare_fluxes(model, a, b, metabolite_list=["a"])
    assert list(result.changed["reaction"]) == ["r1"]


def test_metabolite_list_suppresses_classification_outside_the_filter(model):
    a = pd.Series({"r1": 0.0, "r2": 5.0, "r3": 0.0})
    b = pd.Series({"r1": 5.0, "r2": 0.0, "r3": 0.0})
    result = compare_fluxes(model, a, b, metabolite_list=["D"])
    assert result.turned_on == []
    assert result.turned_off == []
    assert result.changed.empty


def test_unknown_metabolite_is_reported_not_raised(model):
    fluxes = pd.Series({"r1": 1.0, "r2": 1.0, "r3": 1.0})
    result = compare_fluxes(model, fluxes, fluxes, metabolite_list=["D", "nope"])
    assert result.missing_metabolites == ["nope"]


def test_empty_metabolite_list_means_no_filter(model):
    a = pd.Series({"r1": 10.0, "r2": 10.0, "r3": 10.0})
    b = pd.Series({"r1": 20.0, "r2": 20.0, "r3": 20.0})
    result = compare_fluxes(model, a, b, metabolite_list=[])
    assert len(result.changed) == 3
    assert result.n_considered == 3


def test_missing_reaction_in_a_flux_vector_raises(model):
    full = pd.Series({"r1": 1.0, "r2": 1.0, "r3": 1.0})
    partial = pd.Series({"r1": 1.0, "r2": 1.0})
    with pytest.raises(ValueError, match="comparison has no entry for r3"):
        compare_fluxes(model, full, partial)


def test_accepts_a_plain_mapping(model):
    result = compare_fluxes(
        model,
        {"r1": 1.0, "r2": 0.0, "r3": 0.0},
        {"r1": 3.0, "r2": 0.0, "r3": 0.0},
    )
    assert list(result.changed["reaction"]) == ["r1"]
