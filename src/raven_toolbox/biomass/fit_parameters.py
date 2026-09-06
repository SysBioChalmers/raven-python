"""Fit stoichiometric parameters (e.g. maintenance ATP) — port of RAVEN's
``fitParameters``.

Fits one or more unknown coefficients in the stoichiometric matrix (a
maintenance ATP cost is the usual example) by nonlinear least squares: for
each experimental data point, a set of reactions is fixed to measured flux
values, the model is re-solved, and the resulting flux at another set of
reactions is compared against its own measured values. The parameters are
whatever values make those predicted fluxes best match the measurements
across every data point at once.

``plotFitting`` has no port: it is a plotting convenience with no bearing on
the fit itself, and pulling in a plotting library for it would be a poor
trade for something a caller can already do with the returned
``resulting_fluxes`` and their own measurements.
"""
from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import cobra
import numpy as np
import pandas as pd
from scipy.optimize import minimize

__all__ = ["ParameterPosition", "FitParametersResult", "fit_parameters"]


@dataclass(frozen=True)
class ParameterPosition:
    """Where one fitted parameter's coefficient goes in the model.

    Attributes
    ----------
    entries:
        ``(reaction_id, metabolite_id)`` pairs this parameter's magnitude is
        written to.
    negative:
        Per entry, whether that coefficient is ``-parameter`` rather than
        ``+parameter`` (to tell production from consumption at the same
        magnitude) — same length and order as ``entries``.
    """

    entries: tuple[tuple[str, str], ...]
    negative: tuple[bool, ...]

    def __post_init__(self) -> None:
        if len(self.entries) != len(self.negative):
            raise ValueError("entries and negative must have the same length.")


@dataclass
class FitParametersResult:
    """Outcome of a :func:`fit_parameters` call.

    Parameters
    ----------
    parameters:
        Fitted parameter magnitudes, same order as ``parameter_positions``.
        Always non-negative — sign is carried entirely by
        :attr:`ParameterPosition.negative`.
    fitness_score:
        The residual sum of squares at the fitted parameters.
    success:
        Whether the optimizer reports convergence
        (``scipy.optimize.OptimizeResult.success``). RAVEN's own ``exitFlag``
        (``fminsearch``'s 1/0/-1 codes) has no equivalent worth carrying
        over: different optimizer, different convergence bookkeeping.
    model:
        Copy of the input model with the fitted parameters applied.
    resulting_fluxes:
        The fitted model's own predicted flux at ``rxns_to_fit``, one row per
        data point, columns in ``rxns_to_fit`` order — what a caller would
        plot against ``values_to_fit`` to see the fit.
    """

    parameters: np.ndarray
    fitness_score: float
    success: bool
    model: cobra.Model
    resulting_fluxes: pd.DataFrame


def _apply_parameters(
    model: cobra.Model,
    parameters: np.ndarray,
    parameter_positions: Sequence[ParameterPosition],
) -> None:
    for value, position in zip(parameters, parameter_positions, strict=True):
        for (rxn_id, met_id), is_negative in zip(position.entries, position.negative, strict=True):
            rxn = model.reactions.get_by_id(rxn_id)
            met = model.metabolites.get_by_id(met_id)
            target = -value if is_negative else value
            current = rxn.metabolites.get(met, 0.0)
            rxn.add_metabolites({met: target - current}, combine=True)


def _evaluate(
    parameters: np.ndarray,
    model: cobra.Model,
    x_rxns: Sequence[str],
    x_values: np.ndarray,
    rxns_to_fit: Sequence[str],
    values_to_fit: np.ndarray,
    parameter_positions: Sequence[ParameterPosition],
    fit_to_ratio: bool,
) -> tuple[float, np.ndarray, cobra.Model]:
    parameters = np.abs(parameters)
    fitted = model.copy()
    _apply_parameters(fitted, parameters, parameter_positions)

    n_points = x_values.shape[0]
    resulting_fluxes = np.empty((n_points, len(rxns_to_fit)))
    rss = 0.0
    for i in range(n_points):
        with fitted:
            for rxn_id, v in zip(x_rxns, x_values[i], strict=True):
                fitted.reactions.get_by_id(rxn_id).bounds = (float(v), float(v))
            solution = fitted.optimize()
            fit_flux = solution.fluxes.loc[list(rxns_to_fit)].to_numpy()
        resulting_fluxes[i] = fit_flux
        residual = fit_flux / values_to_fit[i] - 1.0 if fit_to_ratio else fit_flux - values_to_fit[i]
        rss += float(residual @ residual)

    return rss, resulting_fluxes, fitted


def fit_parameters(
    model: cobra.Model,
    x_rxns: Sequence[str],
    x_values: np.ndarray,
    rxns_to_fit: Sequence[str],
    values_to_fit: np.ndarray,
    parameter_positions: Sequence[ParameterPosition],
    *,
    fit_to_ratio: bool = True,
    initial_guess: Sequence[float] | None = None,
) -> FitParametersResult:
    """Fit stoichiometric parameters against experimental flux data.

    Parameters
    ----------
    model:
        Model to fit against. Never modified — the fitted model is returned
        separately.
    x_rxns:
        Ids of the reactions that get fixed (both bounds) to a measured
        value at each data point.
    x_values:
        Measured values for ``x_rxns``, shaped *data points × len(x_rxns)*,
        columns in ``x_rxns`` order.
    rxns_to_fit:
        Ids of the reactions whose predicted flux is compared against
        measurements.
    values_to_fit:
        Measured values for ``rxns_to_fit``, shaped *data points ×
        len(rxns_to_fit)*, columns in ``rxns_to_fit`` order.
    parameter_positions:
        One :class:`ParameterPosition` per fitted parameter.
    fit_to_ratio:
        Fit the ratio of predicted to measured value instead of the absolute
        difference, so a large flux does not dominate the residual sum of
        squares just by being large (default True).
    initial_guess:
        Starting value for each parameter (default all ones).

    Returns
    -------
    FitParametersResult
    """
    missing_x = [r for r in x_rxns if r not in model.reactions]
    if missing_x:
        raise ValueError(f"x_rxns not found in the model: {missing_x}")
    missing_fit = [r for r in rxns_to_fit if r not in model.reactions]
    if missing_fit:
        raise ValueError(f"rxns_to_fit not found in the model: {missing_fit}")

    x_values = np.atleast_2d(np.asarray(x_values, dtype=float))
    values_to_fit = np.atleast_2d(np.asarray(values_to_fit, dtype=float))
    guess = np.ones(len(parameter_positions)) if initial_guess is None else np.asarray(initial_guess, dtype=float)

    def objective(params: np.ndarray) -> float:
        rss, _, _ = _evaluate(
            params, model, x_rxns, x_values, rxns_to_fit, values_to_fit, parameter_positions, fit_to_ratio
        )
        return rss

    opt = minimize(objective, guess, method="Nelder-Mead")
    parameters = np.abs(opt.x)
    fitness_score, resulting_fluxes, fitted_model = _evaluate(
        parameters, model, x_rxns, x_values, rxns_to_fit, values_to_fit, parameter_positions, fit_to_ratio
    )

    return FitParametersResult(
        parameters=parameters,
        fitness_score=fitness_score,
        success=bool(opt.success),
        model=fitted_model,
        resulting_fluxes=pd.DataFrame(resulting_fluxes, columns=list(rxns_to_fit)),
    )
