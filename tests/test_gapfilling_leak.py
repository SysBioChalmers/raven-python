"""Tests for stoichiometric leak detection (gapfilling/leak.py, findLeakMetabolite port)."""
import cobra
import pytest

from raven_toolbox.gapfilling import LeakMetaboliteResult, find_leak_metabolite


@pytest.fixture
def leak_model():
    """A produce-only leak (A/B) and a consume-only paired leak (C/D), disjoint.

    R1: A -> B, R2: B -> 2 A. Running the R1/R2 cycle once nets +1 A from
    nothing -- a classic unbalanced-cycle leak. Producing a net unit of B
    instead would need twice the cycle throughput (2 A -> 2 B -> 4 A, net
    +2 A wasted along the way), so the minimal-flux solution always prefers
    A. Neither reaction can support a *consume*-direction leak at all: with
    no relaxation there, running either one forces the other to run at an
    equal rate just to keep the other metabolite balanced, and that leaves
    nothing left over to feed a probe.

    R3: C + D -> (consumes both, nothing produces or removes either
    otherwise), so a consume-direction leak can only be demonstrated by
    probing C and D together -- neither one alone can feed R3 without the
    other.
    """
    m = cobra.Model("leak")
    a = cobra.Metabolite("A", compartment="c")
    b = cobra.Metabolite("B", compartment="c")
    c = cobra.Metabolite("C", compartment="c")
    d = cobra.Metabolite("D", compartment="c")
    m.add_metabolites([a, b, c, d])

    r1 = cobra.Reaction("R1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({a: -1, b: 1})
    r2 = cobra.Reaction("R2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({b: -1, a: 2})
    r3 = cobra.Reaction("R3", lower_bound=0, upper_bound=1000)
    r3.add_metabolites({c: -1, d: -1})
    m.add_reactions([r1, r2, r3])
    return m


def test_produce_finds_the_cheaper_metabolite(leak_model):
    res = find_leak_metabolite(leak_model, "produce")
    assert isinstance(res, LeakMetaboliteResult)
    assert res.status == "optimal"
    assert res.metabolites == ["A"]
    assert res.fluxes["R1"] == pytest.approx(1.0)
    assert res.fluxes["R2"] == pytest.approx(1.0)


def test_consume_reports_a_metabolite_pair(leak_model):
    res = find_leak_metabolite(leak_model, "consume")
    assert res.status == "optimal"
    assert set(res.metabolites) == {"C", "D"}
    assert res.fluxes["R3"] == pytest.approx(0.5)


def test_ignoring_every_metabolite_is_infeasible(leak_model):
    res = find_leak_metabolite(
        leak_model, "produce", ignore_mets=["A", "B", "C", "D"],
    )
    assert res.status == "infeasible"
    assert res.metabolites == []
    assert res.fluxes.empty


def test_ignore_mets_by_name(leak_model):
    leak_model.metabolites.get_by_id("A").name = "alanine"
    res = find_leak_metabolite(
        leak_model, "produce", ignore_mets=["alanine"], is_names=True,
    )
    # With A excluded, the cheapest remaining produce-leak is B itself.
    assert res.status == "optimal"
    assert res.metabolites == ["B"]


def test_invalid_direction_raises(leak_model):
    with pytest.raises(ValueError, match="direction"):
        find_leak_metabolite(leak_model, "sideways")


def test_min_nr_fluxes_reaches_the_same_wiring(leak_model):
    # min_nr_fluxes minimizes active-reaction *count*, not flux sum: A costs
    # less flux (4 vs 5) but both routes activate exactly four reactions, a
    # genuine tie the MILP is free to break either way (same caveat already
    # recorded against get_min_nr_fluxes's own active-set ties) -- this only
    # checks the min_nr_fluxes path is wired up and still finds a real,
    # single-metabolite leak, not which of the tied answers comes back.
    minimal = find_leak_metabolite(leak_model, "produce", min_nr_fluxes=True)
    assert minimal.status == "optimal"
    assert minimal.metabolites in (["A"], ["B"])


def test_original_model_is_not_modified(leak_model):
    n_rxns, n_mets = len(leak_model.reactions), len(leak_model.metabolites)
    find_leak_metabolite(leak_model, "produce")
    assert len(leak_model.reactions) == n_rxns
    assert len(leak_model.metabolites) == n_mets
