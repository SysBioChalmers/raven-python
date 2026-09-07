"""Tests for the highest-flux-fraction path search (analysis/trace_flux_path.py,
traceFluxPath port)."""
import cobra
import pandas as pd
import pytest

from raven_toolbox.analysis import TraceFluxPathResult, trace_flux_path


@pytest.fixture
def junction_model():
    """R1 produces B and ATP together; B splits 70/30 into R2a/R2b, which
    rejoin at D via R3/R4. ATP instead heads straight to a dead end (R5) that
    trace_material should make unreachable by default.

    R1: A -> B + ATP           flux 10
    R2a: B -> C                flux 7   (70% of B)
    R2b: B -> E                flux 3   (30% of B)
    R3: C -> D                 flux 7
    R4: E -> D                 flux 3
    R5: ATP -> D2              flux 10  (currency-only route, no material link)
    """
    m = cobra.Model("junction")
    a = cobra.Metabolite("A", compartment="c")
    b = cobra.Metabolite("B", compartment="c")
    c = cobra.Metabolite("C", compartment="c")
    e = cobra.Metabolite("E", compartment="c")
    d = cobra.Metabolite("D", compartment="c")
    d2 = cobra.Metabolite("D2", compartment="c")
    atp = cobra.Metabolite("atp_c", name="ATP", compartment="c")
    m.add_metabolites([a, b, c, e, d, d2, atp])

    r1 = cobra.Reaction("R1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({a: -1, b: 1, atp: 1})
    r2a = cobra.Reaction("R2a", lower_bound=0, upper_bound=1000)
    r2a.add_metabolites({b: -1, c: 1})
    r2b = cobra.Reaction("R2b", lower_bound=0, upper_bound=1000)
    r2b.add_metabolites({b: -1, e: 1})
    r3 = cobra.Reaction("R3", lower_bound=0, upper_bound=1000)
    r3.add_metabolites({c: -1, d: 1})
    r4 = cobra.Reaction("R4", lower_bound=0, upper_bound=1000)
    r4.add_metabolites({e: -1, d: 1})
    r5 = cobra.Reaction("R5", lower_bound=0, upper_bound=1000)
    r5.add_metabolites({atp: -1, d2: 1})
    m.add_reactions([r1, r2a, r2b, r3, r4, r5])

    fluxes = pd.Series(
        {"R1": 10.0, "R2a": 7.0, "R2b": 3.0, "R3": 7.0, "R4": 3.0, "R5": 10.0}
    )
    return m, fluxes


def test_finds_direct_junction_with_correct_fraction(junction_model):
    model, fluxes = junction_model
    result = trace_flux_path(model, fluxes, "R1", "R2a")
    assert isinstance(result, TraceFluxPathResult)
    assert result.reactions == ["R1", "R2a"]
    assert result.metabolites == ["B"]
    assert result.cumulative_fraction == pytest.approx(0.7)


def test_fractions_compound_across_hops(junction_model):
    model, fluxes = junction_model
    result = trace_flux_path(model, fluxes, "R1", "R3")
    assert result.reactions == ["R1", "R2a", "R3"]
    assert result.metabolites == ["B", "C"]
    assert result.cumulative_fraction == pytest.approx(0.7)


def test_minority_branch_still_found(junction_model):
    model, fluxes = junction_model
    result = trace_flux_path(model, fluxes, "R1", "R4")
    assert result.reactions == ["R1", "R2b", "R4"]
    assert result.cumulative_fraction == pytest.approx(0.3)


def test_trace_material_blocks_the_currency_route_by_default(junction_model):
    model, fluxes = junction_model
    result = trace_flux_path(model, fluxes, "R1", "R5")
    assert result.reactions == []
    assert result.metabolites == []
    assert result.cumulative_fraction == 0.0


def test_disabling_trace_material_finds_the_currency_route(junction_model):
    model, fluxes = junction_model
    result = trace_flux_path(model, fluxes, "R1", "R5", trace_material=False)
    assert result.reactions == ["R1", "R5"]
    assert result.metabolites == ["atp_c"]
    assert result.cumulative_fraction == pytest.approx(1.0)


def test_max_hops_prunes_a_path_that_needs_more(junction_model):
    model, fluxes = junction_model
    result = trace_flux_path(model, fluxes, "R1", "R3", max_hops=1)
    assert result.reactions == []


def test_no_path_between_unconnected_reactions(junction_model):
    model, fluxes = junction_model
    # R2a and R2b are parallel branches from B; neither reaches the other.
    result = trace_flux_path(model, fluxes, "R2a", "R2b")
    assert result.reactions == []
    assert result.cumulative_fraction == 0.0


def test_same_reaction_is_a_trivial_path(junction_model):
    model, fluxes = junction_model
    result = trace_flux_path(model, fluxes, "R1", "R1")
    assert result.reactions == ["R1"]
    assert result.metabolites == []
    assert result.cumulative_fraction == 1.0


def test_exclude_mets_extends_the_currency_list():
    m = cobra.Model("exclude")
    a = cobra.Metabolite("A", compartment="c")
    x = cobra.Metabolite("X", compartment="c")
    d = cobra.Metabolite("D", compartment="c")
    m.add_metabolites([a, x, d])
    r1 = cobra.Reaction("R1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({a: -1, x: 1})
    r2 = cobra.Reaction("R2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({x: -1, d: 1})
    m.add_reactions([r1, r2])
    fluxes = pd.Series({"R1": 5.0, "R2": 5.0})

    assert trace_flux_path(m, fluxes, "R1", "R2").reactions == ["R1", "R2"]
    assert trace_flux_path(m, fluxes, "R1", "R2", exclude_mets=["X"]).reactions == []


def test_carbon_only_drops_a_carbon_free_intermediate():
    m = cobra.Model("carbon")
    a = cobra.Metabolite("A", compartment="c", formula="C6H12O6")
    ion = cobra.Metabolite("ion_c", name="Sodium", compartment="c", formula="Na")
    d = cobra.Metabolite("D", compartment="c", formula="C6H12O6")
    m.add_metabolites([a, ion, d])
    r1 = cobra.Reaction("R1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({a: -1, ion: 1})
    r2 = cobra.Reaction("R2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({ion: -1, d: 1})
    m.add_reactions([r1, r2])
    fluxes = pd.Series({"R1": 5.0, "R2": 5.0})

    assert trace_flux_path(m, fluxes, "R1", "R2").reactions == ["R1", "R2"]
    assert trace_flux_path(m, fluxes, "R1", "R2", carbon_only=True).reactions == []


def test_unknown_reaction_raises(junction_model):
    model, fluxes = junction_model
    with pytest.raises(ValueError, match="not found"):
        trace_flux_path(model, fluxes, "nope", "R1")


def test_missing_flux_entry_raises():
    m = cobra.Model("small")
    a = cobra.Metabolite("A", compartment="c")
    b = cobra.Metabolite("B", compartment="c")
    m.add_metabolites([a, b])
    r1 = cobra.Reaction("R1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({a: -1, b: 1})
    r2 = cobra.Reaction("R2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({b: -1})
    m.add_reactions([r1, r2])

    # Both reactions exist (so resolution succeeds), but the flux vector has
    # nothing for R1, which the search needs as soon as it starts expanding it.
    with pytest.raises(ValueError, match="R1"):
        trace_flux_path(m, {}, "R1", "R2")
