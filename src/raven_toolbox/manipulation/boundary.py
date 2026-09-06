"""Close a model's exchange reactions for metabolic-task checking (RAVEN ``closeModel``)."""
from __future__ import annotations

import cobra


def _is_boundary_reaction(rxn: cobra.Reaction) -> bool:
    """Whether ``rxn`` has metabolites on only one side (RAVEN's rule, matching ``getExchangeRxns``)."""
    coeffs = list(rxn.metabolites.values())
    if not coeffs:
        return False
    return all(c >= 0 for c in coeffs) or all(c <= 0 for c in coeffs)


def close_model(model: cobra.Model) -> cobra.Model:
    """Return a copy with the exchange reactions closed, RAVEN ``closeModel``-style.

    RAVEN's ``closeModel`` adds a boundary metabolite to every reaction with metabolites on
    only one side — no substrates or no products, an exchange / sink / demand reaction —
    balancing it against a boundary metabolite so it can no longer carry flux. After it, a
    metabolic task's inputs and outputs are defined *solely* by the task constraints (as
    ``checkTasks`` assumes), not by leftover open exchanges.

    This closes exactly that set of reactions, which for a single-metabolite reaction is the
    same set as cobra's own ``Reaction.boundary``, but also includes a reaction that consumes
    or produces several metabolites at once without a balancing side.

    Modifies a copy; the original is untouched.
    """
    out = model.copy()
    close_model_in_place(out)
    return out


def close_model_in_place(model: cobra.Model) -> list[str]:
    """Close the boundary reactions of ``model`` in place; return their ids.

    In-place variant of :func:`close_model` for callers that already hold a working copy.
    """
    closed: list[str] = []
    for rxn in model.reactions:
        if _is_boundary_reaction(rxn):
            rxn.bounds = (0.0, 0.0)
            closed.append(rxn.id)
    return closed
