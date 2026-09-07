"""Core Routes for the local MagneticTB workbench."""

from __future__ import annotations

from fastapi import APIRouter

from .route_dependencies import (
    Any,
    Body,
    Dict,
    GroupAlgebra,
    JSONResponse,
    Mapping,
    _bounded,
    _compute,
    common_kernel,
    cyclotomic_context,
    geometry,
    guard_json_tree,
    linear_algebra,
    null_space,
    partial,
    representation,
    tight_binding,
)


router = APIRouter()


def _group_dispatch(payload: Mapping[str, Any]) -> Dict[str, Any]:
    action = str(payload.get("action", "inspect"))
    group = GroupAlgebra(payload["multiplication_table"])
    if action == "inspect":
        return {
            "order": group.order,
            "identity": group.identity,
            "inverse_indices": group.inverse_indices,
            "multiplication_table": group.multiplication_table,
        }
    arguments = payload.get("arguments", {})
    if not isinstance(arguments, Mapping):
        raise TypeError("group arguments must be an object")
    calls = {
        "product": lambda: group.product(arguments["left"], arguments["right"]),
        "generated_subgroup": lambda: group.generated_subgroup(arguments["generators"]),
        "find_generator_indices": lambda: group.find_generator_indices(arguments.get("subgroup")),
        "right_coset": lambda: group.right_coset(arguments["subgroup"], arguments["representative"]),
        "left_coset": lambda: group.left_coset(arguments["subgroup"], arguments["representative"]),
        "schreier_decomposition": lambda: group.schreier_decomposition(
            arguments["subgroup"], arguments["representatives"], arguments["operation"], arguments["source"]
        ),
    }
    if action in calls:
        return {"action": action, "result": calls[action]()}
    if action.startswith("action."):
        compiled = group.compile_action(payload["action_table"])
        name = action.split(".", 1)[1]
        action_calls = {
            "inspect": lambda: {
                "action_table": compiled.action_table,
                "object_count": compiled.object_count,
                "faithful": compiled.faithful,
            },
            "image": lambda: compiled.image(arguments["operation"], arguments["source"]),
            "orbit": lambda: compiled.orbit(arguments["source"]),
            "orbits": compiled.orbits,
            "stabilizer": lambda: compiled.stabilizer(arguments["source"]),
            "transporters": lambda: compiled.transporters(arguments["source"], arguments["target"]),
        }
        if name in action_calls:
            return {"action": action, "result": action_calls[name]()}
    raise ValueError(f"unsupported group action {action!r}")


def _core_dispatch(domain: str, payload: Dict[str, Any]) -> Any:
    guard_json_tree(payload)
    request = dict(payload)
    if domain == "cyclotomic-context":
        return cyclotomic_context(request["conductor"], request.get("max_degree", 128))
    if domain == "null-space":
        return null_space(request["context"], request["matrix"])
    if domain == "common-kernel":
        return common_kernel(request)
    if domain == "linear-algebra":
        operation = request.pop("operation")
        context = request.pop("context")
        return linear_algebra(operation, context, **request)
    if domain == "representation":
        operation = request.pop("operation")
        return representation(operation, **request)
    if domain == "geometry":
        operation = request.pop("operation")
        return geometry(operation, **request)
    if domain == "tight-binding":
        operation = request.pop("operation")
        context = request.pop("context")
        return tight_binding(operation, context, **request)
    if domain == "group":
        return _group_dispatch(request)
    raise ValueError("unknown Rust dispatcher domain")


@router.post("/api/core/{domain}")
async def core_dispatch(domain: str, payload: Dict[str, Any] = Body(...)) -> JSONResponse:
    return _bounded(await _compute(partial(_core_dispatch, domain, payload)))
