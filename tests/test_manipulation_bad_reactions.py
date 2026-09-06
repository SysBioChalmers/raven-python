"""Tests for leak-driven reaction removal (manipulation/bad_reactions.py, removeBadRxns port)."""
import cobra
import pytest

from raven_toolbox.gapfilling import find_leak_metabolite
from raven_toolbox.manipulation import remove_bad_reactions


def _cycle_model(a_formula: str | None, b_formula: str | None) -> cobra.Model:
    """A -> B -> 2A: the same unbalanced-cycle leak used in test_gapfilling_leak.py."""
    m = cobra.Model("bad")
    a = cobra.Metabolite("A", compartment="c", formula=a_formula)
    b = cobra.Metabolite("B", compartment="c", formula=b_formula)
    m.add_metabolites([a, b])
    r1 = cobra.Reaction("R1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({a: -1, b: 1})
    r2 = cobra.Reaction("R2", lower_bound=0, upper_bound=1000)
    r2.add_metabolites({b: -1, a: 2})
    m.add_reactions([r1, r2])
    return m


def test_removes_the_elementally_unbalanced_reaction():
    # R1 (1 C = 1 C) is balanced; R2 (1 C -> 2 C) is not -- the default
    # balance_elements includes carbon, so this is unambiguous.
    model = _cycle_model("C1", "C1")
    removed = remove_bad_reactions(model)
    assert removed == ["R2"]
    assert "R1" in model.reactions
    assert "R2" not in model.reactions
    # Removing the culprit actually fixes the leak.
    assert find_leak_metabolite(model, "produce").status == "infeasible"


def test_rxn_rules_1_gives_up_when_only_an_untracked_element_is_off():
    # Both reactions balance in every element remove_bad_reactions checks by
    # default (C, P, S, N, O) -- the only real imbalance is in hydrogen,
    # which nobody asked to track.
    model = _cycle_model("H1", "H1")
    with pytest.warns(UserWarning, match="No unbalanced reactions"):
        removed = remove_bad_reactions(model, rxn_rules=1)
    assert removed == []
    assert {"R1", "R2"} <= {r.id for r in model.reactions}


def test_rxn_rules_2_also_gives_up_on_an_untracked_element():
    # rxn_rules=2 only additionally catches *unknown* (missing-formula)
    # reactions; an imbalance that is merely outside balance_elements is
    # neither unbalanced-for-our-purposes nor unknown, so this still aborts.
    model = _cycle_model("H1", "H1")
    with pytest.warns(UserWarning, match="No unbalanced or unparsable"):
        removed = remove_bad_reactions(model, rxn_rules=2)
    assert removed == []


def test_rxn_rules_3_removes_one_reaction_arbitrarily_and_fixes_the_leak():
    model = _cycle_model("H1", "H1")
    with pytest.warns(UserWarning, match="error in the metabolite formulas"):
        removed = remove_bad_reactions(model, rxn_rules=3)
    assert removed in (["R1"], ["R2"])
    assert find_leak_metabolite(model, "produce").status == "infeasible"


def test_rxn_rules_1_ignores_a_reaction_with_no_formula():
    # No formula at all -> "unknown", not "unbalanced"; rxn_rules=1 only
    # acts on determinable imbalances, so this aborts without removing R2.
    model = _cycle_model(None, None)
    with pytest.warns(UserWarning, match="No unbalanced reactions"):
        removed = remove_bad_reactions(model, rxn_rules=1)
    assert removed == []


def test_rxn_rules_2_removes_a_reaction_with_no_formula():
    # Both R1 and R2 touch only A/B, and neither has a formula, so both read
    # as "unknown" -- rxn_rules=2's random choice can land on either one,
    # and (as with the rxn_rules=3 tie) removing either fixes the leak.
    model = _cycle_model(None, None)
    removed = remove_bad_reactions(model, rxn_rules=2)
    assert removed in (["R1"], ["R2"])
    assert find_leak_metabolite(model, "produce").status == "infeasible"


def test_consume_direction_reports_a_pair_leak():
    # Same fixture as find_leak_metabolite's own "exchanged in pairs" case:
    # C + D -> consumes both, and nothing else touches either one.
    m = cobra.Model("bad")
    c = cobra.Metabolite("C", compartment="c")
    d = cobra.Metabolite("D", compartment="c")
    m.add_metabolites([c, d])
    r3 = cobra.Reaction("R3", lower_bound=0, upper_bound=1000)
    r3.add_metabolites({c: -1, d: -1})
    m.add_reactions([r3])

    removed = remove_bad_reactions(m, rxn_rules=2)
    assert removed == ["R3"]
    assert find_leak_metabolite(m, "consume").status == "infeasible"


def test_infeasible_model_raises():
    m = cobra.Model("infeasible")
    x = cobra.Metabolite("X", compartment="c")
    m.add_metabolites([x])
    demand = cobra.Reaction("demand", lower_bound=5, upper_bound=5)
    demand.add_metabolites({x: -1})
    m.add_reactions([demand])

    with pytest.raises(ValueError, match="not feasible"):
        remove_bad_reactions(m)


def test_open_exchange_reactions_warn():
    m = cobra.Model("open")
    x = cobra.Metabolite("X", compartment="c")
    m.add_metabolites([x])
    exch = cobra.Reaction("EX_X", lower_bound=-10, upper_bound=10)
    exch.add_metabolites({x: -1})
    m.add_reactions([exch])

    with pytest.warns(UserWarning, match="open exchange reactions"):
        removed = remove_bad_reactions(m)
    assert removed == []


def test_ignore_mets_prevents_the_leak_from_being_found():
    model = _cycle_model("C1", "C1")
    removed = remove_bad_reactions(model, ignore_mets=["A", "B"])
    assert removed == []
    assert {"R1", "R2"} <= {r.id for r in model.reactions}
