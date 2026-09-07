"""Stable-style init orchestration for the local browser workbench."""

from __future__ import annotations

import io
from contextlib import redirect_stdout
from typing import Any, Dict

from ..modeling import (
    DirectProductRepresentation,
    MSG,
    PreparedModel,
    Wyckoff,
    _CURRENT_MODEL,
    current_model,
    init,
    initfromrep,
)
from .data_options_service import (
    _bravais_matrix_to_user,
    _bravais_option,
    _named_parameters,
    _parameter_seed,
    _selected_wyckoff_option,
)
from .form_parsing import (
    _field_from_form,
    _operation_from_form,
    _parse_lattice_matrix,
    _parse_lattice_scalar,
    _parse_representation_matrix,
    _parse_vector,
)
from .model_schemas import FriendlyInitRequest
from .runtime import _catalog


def _friendly_init(payload: FriendlyInitRequest) -> tuple[PreparedModel, str]:
    arguments: Dict[str, Any] = {
        "basis_functions": tuple(tuple(orbit.basis_functions) for orbit in payload.orbits),
        "initial_bond_shells": payload.initial_bond_shells,
        "generate_symmetry_group": payload.generate_symmetry_group,
        "representation_mode": payload.representation_mode,
        "exact_field": _field_from_form(payload.exact_field),
    }
    if payload.symmetry.source == "msg":
        if not payload.symmetry.msg_id:
            raise ValueError("MSG symmetry requires one selected stable MSG id")
        if payload.lattice.preset != "msg_bravais":
            raise ValueError("MSG models require the Bravais lattice fixed by the selected group")
        msg_record = _catalog().msg(payload.symmetry.msg_id)
        bravais_id = str(msg_record["bravais_lattice_id"])
        bravais = _bravais_option(bravais_id)
        arguments["lattice"] = _bravais_matrix_to_user(
            _catalog().bravais(bravais_id)["primitive_vectors"]
        )
        arguments["lattpar"] = _named_parameters(
            payload.lattice.parameters,
            bravais["parameter_symbols"],
            "Bravais lattice parameters",
        )
        arguments["symminformation"] = MSG(stable_id=payload.symmetry.msg_id)
        selections = []
        for index, orbit in enumerate(payload.orbits, start=1):
            if orbit.source_ordinal is None:
                raise ValueError(f"orbit {index} requires a Wyckoff source ordinal")
            if orbit.position is not None or orbit.moment is not None:
                raise ValueError(
                    f"orbit {index} is Data-constrained; provide only its named free parameters"
                )
            option = _selected_wyckoff_option(
                payload.symmetry.msg_id,
                orbit.source_ordinal,
                orbit.letter,
            )
            selections.append(
                Wyckoff(
                    orbit.letter,
                    position=_parameter_seed(
                        orbit.coordinate_parameters,
                        option["coordinate_symbols"],
                        ("x", "y", "z"),
                        f"orbit {index} coordinate parameters",
                    ),
                    moment=_parameter_seed(
                        orbit.moment_parameters,
                        option["moment_symbols"],
                        ("mx", "my", "mz"),
                        f"orbit {index} moment parameters",
                    ),
                    source_ordinal=orbit.source_ordinal,
                )
            )
        arguments["wyckoffposition"] = tuple(selections)
    else:
        if not payload.symmetry.operations:
            raise ValueError("explicit symmetry requires at least one operation")
        if payload.lattice.preset == "msg_bravais":
            raise ValueError("explicit symmetry requires an explicitly selected lattice")
        if payload.lattice.preset == "simple_cubic":
            arguments["lattice"] = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
        elif payload.lattice.preset == "custom":
            if payload.lattice.matrix is None:
                raise ValueError("custom lattice requires a 3 by 3 matrix")
            arguments["lattice"] = _parse_lattice_matrix(
                payload.lattice.matrix, "lattice"
            )
        if payload.lattice.parameters:
            arguments["lattpar"] = {
                str(name): _parse_lattice_scalar(value, f"lattice parameter {name}")
                for name, value in payload.lattice.parameters.items()
            }
        for index, orbit in enumerate(payload.orbits, start=1):
            if orbit.coordinate_parameters or orbit.moment_parameters:
                raise ValueError(
                    f"explicit orbit {index} uses position/moment seeds, not Data parameter bindings"
                )
        arguments["symminformation"] = tuple(
            _operation_from_form(operation) for operation in payload.symmetry.operations
        )
        arguments["wyckoffposition"] = tuple(
            (
                _parse_vector(
                    orbit.position
                    if orbit.position is not None
                    else _missing_explicit_seed(index, "position"),
                    f"orbit {index} position",
                ),
                _parse_vector(
                    orbit.moment
                    if orbit.moment is not None
                    else _missing_explicit_seed(index, "moment"),
                    f"orbit {index} moment",
                ),
            )
            for index, orbit in enumerate(payload.orbits, start=1)
        )
    initializer = init
    if payload.direct_representation is not None:
        if payload.representation_mode != "DirectProduct":
            raise ValueError(
                "direct_representation requires representation_mode=DirectProduct"
            )
        if payload.generate_symmetry_group:
            raise ValueError(
                "initfromrep requires an explicit finite operation list, not generated generators"
            )
        representation = payload.direct_representation
        arguments.pop("basis_functions")
        arguments.pop("generate_symmetry_group")
        arguments["repinformation"] = DirectProductRepresentation(
            matrices_by_orbit=tuple(
                tuple(
                    _parse_representation_matrix(
                        matrix,
                        f"direct representation orbit {orbit_index} matrix {matrix_index}",
                    )
                    for matrix_index, matrix in enumerate(matrices, start=1)
                )
                for orbit_index, matrices in enumerate(
                    representation.matrices_by_orbit, start=1
                )
            ),
            continuous_site_matrices=(
                tuple(
                    tuple(
                        _parse_representation_matrix(
                            matrix,
                            f"continuous representation orbit {orbit_index} site {site_index}",
                        )
                        for site_index, matrix in enumerate(matrices, start=1)
                    )
                    for orbit_index, matrices in enumerate(
                        representation.continuous_site_matrices, start=1
                    )
                )
                if representation.continuous_site_matrices is not None
                else None
            ),
        )
        arguments["orbital_labels"] = tuple(
            tuple(labels) for labels in representation.orbital_labels
        )
        initializer = initfromrep
    output = io.StringIO()
    token = _CURRENT_MODEL.set(None)
    try:
        with redirect_stdout(output):
            initializer(**arguments)
            model = current_model()
    finally:
        _CURRENT_MODEL.reset(token)
    return model, output.getvalue()


def _missing_explicit_seed(index: int, field: str) -> List[Any]:
    raise ValueError(f"explicit orbit {index} requires a three-component {field} seed")
