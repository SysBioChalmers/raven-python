"""Tests for stoichiometric parameter fitting (biomass/fit_parameters.py,
fitParameters port)."""
import cobra
import numpy as np
import pytest

from raven_toolbox.biomass import FitParametersResult, ParameterPosition, fit_parameters


@pytest.fixture
def chain_model():
    """glc -> (R1) -> atp + x; atp -> (MAINT, forced drain); k*x -> (GROWTH).

    x has exactly one producer (R1, 1:1 with glc) and one consumer (GROWTH,
    at coefficient -k), so mass balance alone -- no objective needed --
    pins GROWTH's flux to glc_uptake / k for any k. That makes the correct
    fitted k analytically exact: no alternate optima, no solver ambiguity.
    """
    m = cobra.Model("fit_toy")
    glc = cobra.Metabolite("glc", compartment="c")
    atp = cobra.Metabolite("atp", compartment="c")
    x = cobra.Metabolite("x", compartment="c")
    m.add_metabolites([glc, atp, x])

    ex_glc = cobra.Reaction("EX_glc", lower_bound=0, upper_bound=1000)
    ex_glc.add_metabolites({glc: 1})
    r1 = cobra.Reaction("R1", lower_bound=0, upper_bound=1000)
    r1.add_metabolites({glc: -1, atp: 1, x: 1})
    maint = cobra.Reaction("MAINT", lower_bound=0, upper_bound=1000)
    maint.add_metabolites({atp: -1})
    growth = cobra.Reaction("GROWTH", lower_bound=0, upper_bound=1000)
    growth.add_metabolites({x: -1.0})  # placeholder; fit_parameters overwrites this
    m.add_reactions([ex_glc, r1, maint, growth])
    return m


def test_recovers_the_true_parameter(chain_model):
    true_k = 2.0
    x_values = [[2.0], [4.0], [6.0]]
    values_to_fit = [[v[0] / true_k] for v in x_values]  # noiseless: [1, 2, 3]
    positions = [ParameterPosition((("GROWTH", "x"),), (True,))]

    result = fit_parameters(
        chain_model, ["EX_glc"], x_values, ["GROWTH"], values_to_fit, positions,
        initial_guess=[1.0],
    )

    assert isinstance(result, FitParametersResult)
    assert result.parameters[0] == pytest.approx(true_k, abs=1e-4)
    assert result.fitness_score == pytest.approx(0.0, abs=1e-6)
    assert result.success
    np.testing.assert_allclose(
        result.resulting_fluxes["GROWTH"].to_numpy(), [1.0, 2.0, 3.0], atol=1e-4,
    )


def test_fitted_model_carries_the_parameter(chain_model):
    positions = [ParameterPosition((("GROWTH", "x"),), (True,))]
    result = fit_parameters(
        chain_model, ["EX_glc"], [[2.0], [4.0]], ["GROWTH"], [[1.0], [2.0]], positions,
    )
    fitted_growth = result.model.reactions.get_by_id("GROWTH")
    fitted_x = result.model.metabolites.get_by_id("x")
    assert fitted_growth.metabolites[fitted_x] == pytest.approx(-2.0, abs=1e-4)
    # The input model is untouched.
    original_growth = chain_model.reactions.get_by_id("GROWTH")
    original_x = chain_model.metabolites.get_by_id("x")
    assert original_growth.metabolites[original_x] == -1.0


def test_fit_to_ratio_false_also_recovers_a_noiseless_fit(chain_model):
    positions = [ParameterPosition((("GROWTH", "x"),), (True,))]
    result = fit_parameters(
        chain_model, ["EX_glc"], [[2.0], [4.0], [6.0]], ["GROWTH"], [[1.0], [2.0], [3.0]],
        positions, fit_to_ratio=False,
    )
    assert result.parameters[0] == pytest.approx(2.0, abs=1e-4)


def test_unknown_x_rxn_raises(chain_model):
    positions = [ParameterPosition((("GROWTH", "x"),), (True,))]
    with pytest.raises(ValueError, match="x_rxns"):
        fit_parameters(chain_model, ["nope"], [[2.0]], ["GROWTH"], [[1.0]], positions)


def test_unknown_rxn_to_fit_raises(chain_model):
    positions = [ParameterPosition((("GROWTH", "x"),), (True,))]
    with pytest.raises(ValueError, match="rxns_to_fit"):
        fit_parameters(chain_model, ["EX_glc"], [[2.0]], ["nope"], [[1.0]], positions)


def test_parameter_position_length_mismatch_raises():
    with pytest.raises(ValueError, match="same length"):
        ParameterPosition((("R", "M"),), (True, False))
