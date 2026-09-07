"""Tests for close_model / close_model_in_place (RAVEN closeModel port)."""
import cobra
import pytest

from raven_toolbox.manipulation import close_model
from raven_toolbox.manipulation.boundary import close_model_in_place


@pytest.fixture
def model():
    m = cobra.Model("t")
    m.compartments = {"c": "cytoplasm"}
    a = cobra.Metabolite("a_c", name="a", compartment="c")
    b = cobra.Metabolite("b_c", name="b", compartment="c")
    m.add_metabolites([a, b])

    real = cobra.Reaction("R1")  # a -> b, a genuine reaction, not a boundary one
    real.add_metabolites({a: -1, b: 1})
    real.bounds = (0, 1000)

    unit_sink = cobra.Reaction("R2")  # a =>
    unit_sink.add_metabolites({a: -1})
    unit_sink.bounds = (0, 1000)

    scaled_sink = cobra.Reaction("R3")  # 2 a =>
    scaled_sink.add_metabolites({a: -2})
    scaled_sink.bounds = (0, 1000)

    multi_met_sink = cobra.Reaction("R4")  # 0.5 a + 0.5 b =>
    multi_met_sink.add_metabolites({a: -0.5, b: -0.5})
    multi_met_sink.bounds = (0, 1000)

    m.add_reactions([real, unit_sink, scaled_sink, multi_met_sink])
    return m


def test_genuine_reaction_left_open(model):
    closed = close_model_in_place(model)
    assert "R1" not in closed
    assert model.reactions.R1.bounds == (0, 1000)


def test_unit_single_met_sink_closed(model):
    closed = close_model_in_place(model)
    assert "R2" in closed
    assert model.reactions.R2.bounds == (0.0, 0.0)


def test_scaled_single_met_sink_closed(model):
    # A non-unit coefficient must not exempt a single-metabolite sink.
    closed = close_model_in_place(model)
    assert "R3" in closed
    assert model.reactions.R3.bounds == (0.0, 0.0)


def test_multi_met_sink_closed(model):
    # Several metabolites on the one populated side is still a boundary
    # reaction, not a genuine two-sided one.
    closed = close_model_in_place(model)
    assert "R4" in closed
    assert model.reactions.R4.bounds == (0.0, 0.0)


def test_close_model_copies_and_leaves_original_untouched(model):
    out = close_model(model)
    assert out is not model
    assert out.reactions.R2.bounds == (0.0, 0.0)
    assert model.reactions.R2.bounds == (0, 1000)
