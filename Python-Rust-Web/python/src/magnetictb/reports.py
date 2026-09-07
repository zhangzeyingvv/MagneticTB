"""Stable 2.0.10 presentation/query APIs over Rust-owned model results.

The functions in this module select and format already compiled Rust records.
They never construct groups, representations, bond orbits, constraints, or
Hamiltonians in Python.
"""

from __future__ import annotations

from dataclasses import dataclass, is_dataclass
from fractions import Fraction
from html import escape
from math import sqrt as numeric_sqrt
from types import MappingProxyType
from typing import Any, Mapping, Optional, Sequence, Tuple, Union

from .api import DataCatalog, ModelError
from .data import MSG, MagneticGroupSelector, _catalog_exact_to_user
from .linear_algebra import (
    CyclotomicElement,
    ExactExpr,
    ExactMatrix,
    QuadraticNumber,
    _decode_matrix,
    _decode_quadratic,
)
from .model import PreparedModel, current_model
from .symmetry import SpinSpaceOperation, SymmetryOperation


def _immutable_mapping(value: Mapping[str, Any]) -> Mapping[str, Any]:
    return MappingProxyType(dict(value))


def _to_user(value: Any) -> Any:
    if isinstance(value, ExactMatrix):
        return [[_to_user(item) for item in row] for row in value.rows]
    if isinstance(value, CyclotomicElement):
        return str(value).replace("I", "i")
    if isinstance(value, QuadraticNumber):
        return str(value)
    if isinstance(value, ExactExpr):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _to_user(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_to_user(item) for item in value]
    if isinstance(value, list):
        return [_to_user(item) for item in value]
    if is_dataclass(value):
        return str(value)
    return value


@dataclass(frozen=True)
class TableReport:
    """Immutable, language-neutral replacement for a Mathematica ``Grid``."""

    schema: str
    title: str
    columns: Tuple[str, ...]
    rows: Tuple[Mapping[str, Any], ...]
    metadata: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "columns", tuple(self.columns))
        object.__setattr__(
            self,
            "rows",
            tuple(_immutable_mapping(row) for row in self.rows),
        )
        object.__setattr__(self, "metadata", _immutable_mapping(self.metadata))

    def __len__(self) -> int:
        return len(self.rows)

    def __iter__(self):
        return iter(self.rows)

    def __getitem__(self, index: int) -> Mapping[str, Any]:
        return self.rows[index]

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": self.schema,
            "title": self.title,
            "columns": list(self.columns),
            "rows": [_to_user(row) for row in self.rows],
            "metadata": _to_user(self.metadata),
        }


def _model(value: Optional[PreparedModel]) -> PreparedModel:
    if value is None:
        return current_model()
    if not isinstance(value, PreparedModel):
        raise TypeError("model must be PreparedModel or None")
    return value


def orbitalTable(*, model: Optional[PreparedModel] = None) -> TableReport:
    """Return the stable Hamiltonian basis order prepared by ``init``."""

    prepared = _model(model)
    rows = tuple(
        {
            "HamiltonianIndex": record.orbital_id,
            "SiteIndex": record.site_id,
            "WyckoffOrbitIndex": record.wyckoff_orbit,
            "EquivalentAtomIndex": record.equivalent_atom,
            "FractionalPosition": record.fractional_position,
            "CartesianPosition": record.cartesian_position,
            "LocalOrbitalIndex": record.local_orbital,
            "BasisFunctionInput": record.basis_function_input,
            "BasisFunction": record.basis_function,
            "BasisState": record.basis_state,
            "SpatialOrbital": record.spatial_orbital,
            "SpinStructure": record.spin_structure,
            "SpinState": record.spin_state,
            "BasisConvention": record.basis_convention,
            "ReferenceEquivalentAtomIndex": record.reference_equivalent_atom,
            "TransportOperationIndex": record.transport_operation,
            "TransportOperationLabel": record.transport_operation_label,
        }
        for record in prepared.orbital_table
    )
    return TableReport(
        schema="MagneticTBHamiltonianBasisOrder",
        title=(
            f"Hamiltonian dimension = {len(rows)}; "
            f"representation mode = {prepared.representation_mode}"
        ),
        columns=tuple(rows[0]) if rows else (),
        rows=rows,
        metadata={
            "SchemaVersion": 4,
            "Dimension": len(rows),
            "RepresentationMode": prepared.representation_mode,
            "Convention": "H[i,j] = <basis i|H|basis j>",
        },
    )


def showHamiltonianBasis(
    row: Optional[int] = None,
    column: Optional[int] = None,
    *,
    model: Optional[PreparedModel] = None,
) -> TableReport:
    """Show the full ordered basis or the row/column basis for ``H[row,column]``."""

    report = orbitalTable(model=model)
    if row is None and column is None:
        return report
    dimension = int(report.metadata["Dimension"])
    if (
        isinstance(row, bool)
        or not isinstance(row, int)
        or isinstance(column, bool)
        or not isinstance(column, int)
        or row < 1
        or column < 1
        or row > dimension
        or column > dimension
    ):
        raise ModelError(
            "InvalidHamiltonianBasisIndex",
            f"row and column must lie in 1..{dimension}; received ({row}, {column})",
        )
    rows = (
        _immutable_mapping({"Role": "Row (bra)", **dict(report.rows[row - 1])}),
        _immutable_mapping({"Role": "Column (ket)", **dict(report.rows[column - 1])}),
    )
    return TableReport(
        schema="MagneticTBHamiltonianElementBasis",
        title=f"H[{row},{column}] = <row|H|column>",
        columns=("Role",) + report.columns,
        rows=rows,
        metadata={**dict(report.metadata), "Row": row, "Column": column},
    )


def _quadratic_float(value: QuadraticNumber) -> float:
    coefficient = value.coefficients
    return float(coefficient[0]) + float(coefficient[1]) * numeric_sqrt(2.0) + float(
        coefficient[2]
    ) * numeric_sqrt(3.0) + float(coefficient[3]) * numeric_sqrt(6.0)


def bondTable(shell: int, *, model: Optional[PreparedModel] = None) -> TableReport:
    """Read the complete directed-bond data without displaying it.

    ``showbonds`` is the display-only entry point. This data interface retains
    the Rust record order, exact endpoints, and report metadata for callers.
    """

    prepared = _model(model)
    if (
        isinstance(shell, bool)
        or not isinstance(shell, int)
        or shell < 1
        or shell > prepared.initial_bond_shells
    ):
        raise ModelError(
            "ShellNotPrepared",
            f"shell {shell} is outside 1..{prepared.initial_bond_shells}",
        )
    bonds = prepared.bond_shells[shell - 1]
    raw = prepared.to_canonical_dict(shell)
    site_positions: dict[int, Tuple[CyclotomicElement, ...]] = {}
    site_metadata: dict[int, tuple[int, int]] = {}
    for orbital in prepared.orbital_table:
        site_positions.setdefault(orbital.site_id, orbital.fractional_position)
        site_metadata.setdefault(
            orbital.site_id, (orbital.wyckoff_orbit, orbital.equivalent_atom)
        )
    rows = []
    for site_id in sorted(site_positions):
        selected = tuple(bond for bond in bonds if bond.row_site == site_id)
        orbit, equivalent = site_metadata[site_id]
        rows.append(
            {
                "SiteIndex": site_id,
                "OrbitIndex": orbit,
                "EquivalentIndex": equivalent,
                "AtomPosition": site_positions[site_id],
                "BondCount": len(selected),
                "TargetPositions": tuple(bond.column_endpoint for bond in selected),
            }
        )
    squared = bonds[0].squared_distance if bonds else QuadraticNumber((Fraction(0),) * 4)
    distance = numeric_sqrt(max(0.0, _quadratic_float(squared)))
    return TableReport(
        schema="MagneticTBBondShellDisplayData",
        title=(
            f"Bond shell {shell} ({shell - 1}-th neighbour), bond length = {distance:.12g}; "
            f"{len(bonds)} directed bonds"
        ),
        columns=(
            "SiteIndex",
            "OrbitIndex",
            "EquivalentIndex",
            "AtomPosition",
            "BondCount",
            "TargetPositions",
        ),
        rows=tuple(rows),
        metadata={
            "SchemaVersion": 1,
            "Shell": shell,
            "NeighbourOrder": shell - 1,
            "BondLength": distance,
            "SquaredDistance": squared,
            "SiteCount": len(site_positions),
            "DirectedBondCount": len(bonds),
            "DirectedOrbitCount": len(raw["directed_bond_orbits"]),
        },
    )


_BOND_TABLE_HEADERS = (
    "Site", "Orbit", "Equivalent atom", "Atom position", "Directed bonds", "Target positions",
)


def _bond_cell(value: Any) -> str:
    if isinstance(value, (tuple, list)):
        return "(" + ", ".join(_bond_cell(item) for item in value) + ")"
    return str(_to_user(value))


def _bond_table_cells(report: TableReport) -> Tuple[Tuple[str, ...], ...]:
    # Formatting only: never sort endpoints, collapse reverse bonds, or round
    # their exact coordinates. Each target gets its own visible table line.
    return tuple(
        tuple(
            "\n".join(_bond_cell(position) for position in row[column])
            if column == "TargetPositions" else _bond_cell(row[column])
            for column in report.columns
        )
        for row in report.rows
    )


def _bond_table_text(report: TableReport) -> str:
    cells = _bond_table_cells(report)
    widths = [len(header) for header in _BOND_TABLE_HEADERS]
    for row in cells:
        for index, cell in enumerate(row):
            widths[index] = max(widths[index], *(len(line) for line in cell.split("\n")))

    def line(values: Sequence[str]) -> str:
        return " | ".join(value.ljust(width) for value, width in zip(values, widths)).rstrip()

    lines = [report.title, line(_BOND_TABLE_HEADERS), "-+-".join("-" * width for width in widths)]
    for row in cells:
        parts = [cell.split("\n") for cell in row]
        for index in range(max(map(len, parts))):
            lines.append(line(tuple(part[index] if index < len(part) else "" for part in parts)))
    return "\n".join(lines)


def _bond_table_html(report: TableReport) -> str:
    headings = "".join(f'<th scope="col">{escape(header)}</th>' for header in _BOND_TABLE_HEADERS)
    rows = "".join(
        "<tr>" + "".join(
            '<td style="text-align:left;vertical-align:top;padding:0.35em 0.6em">'
            + escape(cell).replace("\n", "<br>") + "</td>" for cell in row
        ) + "</tr>"
        for row in _bond_table_cells(report)
    )
    return (
        '<div style="overflow-x:auto"><table class="magnetictb-bond-table">'
        f"<caption>{escape(report.title)}</caption><thead><tr>{headings}</tr></thead>"
        f"<tbody>{rows}</tbody></table></div>"
    )


def showbonds(shell: int, *, model: Optional[PreparedModel] = None) -> None:
    """Display a prepared shell's bond table and return ``None``.

    Jupyter displays an HTML table; scripts and terminals print a text table.
    Use ``bondTable(shell, model=...)`` to obtain the complete data object.
    No bond search or mathematical calculation is performed here.
    """

    report = bondTable(shell, model=model)
    try:
        from IPython import get_ipython
    except ImportError:
        notebook = None
    else:
        notebook = get_ipython()
    if notebook is not None and getattr(notebook, "kernel", None) is not None:
        from IPython.display import HTML, display

        display(HTML(_bond_table_html(report)))
    else:
        print(_bond_table_text(report))


def _selection_indices(
    selection: Union[str, int, Sequence[int]], prepared: PreparedModel
) -> Tuple[int, ...]:
    operation_count = len(prepared.symmetry_operations)
    if selection == "Automatic":
        indices = tuple(index + 1 for index in prepared.generator_indices)
    elif selection == "All":
        indices = tuple(range(1, operation_count + 1))
    elif isinstance(selection, int) and not isinstance(selection, bool):
        indices = (selection,)
    elif isinstance(selection, Sequence) and not isinstance(selection, (str, bytes)):
        indices = tuple(dict.fromkeys(selection))
    else:
        raise ModelError(
            "InvalidSymmetryRepresentationSelection",
            "selection must be 'Automatic', 'All', a positive integer, or an integer sequence",
        )
    if not indices or any(
        isinstance(index, bool)
        or not isinstance(index, int)
        or index < 1
        or index > operation_count
        for index in indices
    ):
        raise ModelError(
            "SymmetryRepresentationSelectionOutOfRange",
            f"operation selection must lie in 1..{operation_count}",
        )
    return indices


def showSymmetryRepresentations(
    selection: Union[str, int, Sequence[int]] = "Automatic",
    *,
    model: Optional[PreparedModel] = None,
) -> TableReport:
    """Select exact Rust-compiled representation matrices in stable operation order."""

    prepared = _model(model)
    indices = _selection_indices(selection, prepared)
    rows = []
    for number in indices:
        operation = prepared.symmetry_operations[number - 1]
        if isinstance(operation, SpinSpaceOperation):
            rotation = operation.space_rotation
            translation = tuple(operation.space_translation)
        else:
            rotation = operation.rotation
            translation = tuple(operation.translation)
        rows.append(
            {
                "OperationIndex": number,
                "Label": operation.label,
                "Generator": number - 1 in prepared.generator_indices,
                "Rotation": ExactMatrix(tuple(tuple(row) for row in rotation)),
                "Translation": translation,
                "SpinAction": prepared.spin_actions[number - 1],
                "Antiunitary": operation.antiunitary,
                "Matrix": prepared.representation_matrices[number - 1],
            }
        )
    dimension = prepared.representation_matrices[0].shape[0]
    return TableReport(
        schema="MagneticTBSymmetryRepresentationDisplayData",
        title=(
            f"Representation mode: {prepared.representation_mode}; dimension = {dimension}; "
            f"showing {len(indices)} of {len(prepared.symmetry_operations)} operations"
        ),
        columns=tuple(rows[0]) if rows else (),
        rows=tuple(rows),
        metadata={
            "SchemaVersion": 1,
            "RepresentationMode": prepared.representation_mode,
            "Convention": prepared.representation_method,
            "Dimension": dimension,
            "OperationCount": len(prepared.symmetry_operations),
            "GeneratorIndices": tuple(index + 1 for index in prepared.generator_indices),
            "SelectedIndices": indices,
            "UnitaryVerified": True,
        },
    )


def _zero(value: CyclotomicElement) -> bool:
    return all(coefficient == 0 for coefficient in value.coefficients)


def _site_offsets(dimensions: Sequence[int]) -> Tuple[int, ...]:
    result = []
    offset = 0
    for dimension in dimensions:
        result.append(offset)
        offset += int(dimension)
    return tuple(result)


def _hopping_record(
    *,
    prepared: PreparedModel,
    shell: int,
    parameter: str,
    coefficient: CyclotomicElement,
    orbit_index: int,
    representative_bond_index: int,
    bond_index: int,
    row_block: int,
    column_block: int,
    row_local: int,
    column_local: int,
    raw: Mapping[str, Any],
) -> Mapping[str, Any]:
    bond = raw["bonds"][bond_index]
    offsets = _site_offsets(raw["site_dimensions"])
    row_h = offsets[row_block] + row_local
    column_h = offsets[column_block] + column_local
    return {
        "Parameter": parameter,
        "Coefficient": coefficient,
        "Shell": shell,
        "ConstraintOrbitIndex": orbit_index + 1,
        "BondIndex": bond_index + 1,
        "RepresentativeBondIndex": representative_bond_index + 1,
        "Representative": bond_index == representative_bond_index,
        "RowSiteIndex": row_block + 1,
        "RowCell": (0, 0, 0),
        "RowEndpoint": tuple(_decode_quadratic(item) for item in bond["source_endpoint"]),
        "RowLocalOrbitalIndex": row_local + 1,
        "RowHamiltonianIndex": row_h + 1,
        "RowBasisState": prepared.orbital_table[row_h].basis_state,
        "ColumnSiteIndex": column_block + 1,
        "ColumnCell": tuple(int(item) for item in bond["translation"]),
        "ColumnEndpoint": tuple(_decode_quadratic(item) for item in bond["target_endpoint"]),
        "ColumnLocalOrbitalIndex": column_local + 1,
        "ColumnHamiltonianIndex": column_h + 1,
        "ColumnBasisState": prepared.orbital_table[column_h].basis_state,
        "CellTranslation": tuple(int(item) for item in bond["translation"]),
    }


def showHoppingParameters(
    shell: int,
    parameter: Optional[str] = None,
    *,
    model: Optional[PreparedModel] = None,
    hermitian: bool = True,
    kernel_method: str = "Iterative",
    validation_level: str = "Basic",
) -> TableReport:
    """Show representative sources or all real-space occurrences of a parameter."""

    prepared = _model(model)
    space = prepared.hamiltonian_space(
        shell,
        hermitian=hermitian,
        kernel_method=kernel_method,  # type: ignore[arg-type]
        validation_level=validation_level,  # type: ignore[arg-type]
    )
    raw = space.to_canonical_dict()
    field = prepared._field
    rows = []
    if parameter is None:
        for term in space.basis:
            parameter_record = term.parameter
            orbit_index = parameter_record.constraint_orbit - 1
            bond_index = parameter_record.representative_bond - 1
            matrix = term.representative_hoppings[orbit_index]
            bond = raw["bonds"][bond_index]
            row_block = int(bond["source_site"])
            column_block = int(bond["target_site"])
            for row_index, matrix_row in enumerate(matrix.rows):
                for column_index, coefficient in enumerate(matrix_row):
                    if not _zero(coefficient):
                        rows.append(
                            _hopping_record(
                                prepared=prepared,
                                shell=shell,
                                parameter=parameter_record.name,
                                coefficient=coefficient,
                                orbit_index=orbit_index,
                                representative_bond_index=bond_index,
                                bond_index=bond_index,
                                row_block=row_block,
                                column_block=column_block,
                                row_local=row_index,
                                column_local=column_index,
                                raw=raw,
                            )
                        )
        title = f"Independent parameter sources for shell {shell} (representative hoppings)"
    else:
        if not isinstance(parameter, str) or parameter not in space.parameter_names:
            raise ModelError(
                "UnknownHamiltonianParameter",
                f"parameter {parameter!r} is not one of {space.parameter_names}",
            )
        parameter_index = space.parameter_names.index(parameter)
        parameter_record = space.parameters[parameter_index]
        solution = raw["parameter_space"]["parameter_basis_solutions"][parameter_index]
        for solved_term in solution["terms"]:
            bond_index = int(solved_term["bond_index"])
            matrix = _decode_matrix(solved_term["matrix"], field)
            for row_index, matrix_row in enumerate(matrix.rows):
                for column_index, coefficient in enumerate(matrix_row):
                    if not _zero(coefficient):
                        rows.append(
                            _hopping_record(
                                prepared=prepared,
                                shell=shell,
                                parameter=parameter,
                                coefficient=coefficient,
                                orbit_index=int(solved_term["orbit_index"]),
                                representative_bond_index=parameter_record.representative_bond - 1,
                                bond_index=bond_index,
                                row_block=int(solved_term["row_block"]),
                                column_block=int(solved_term["column_block"]),
                                row_local=row_index,
                                column_local=column_index,
                                raw=raw,
                            )
                        )
        title = f"All real-space occurrences of {parameter} in shell {shell}"
    columns = tuple(rows[0]) if rows else (
        "Parameter",
        "ConstraintOrbitIndex",
        "BondIndex",
        "ColumnSiteIndex",
        "RowSiteIndex",
        "CellTranslation",
        "Coefficient",
        "Representative",
    )
    return TableReport(
        schema=(
            "MagneticTBHoppingParameterSourceData"
            if parameter is None
            else "MagneticTBHoppingParameterOccurrenceData"
        ),
        title=title,
        columns=columns,
        rows=tuple(rows),
        metadata={
            "SchemaVersion": 1,
            "Shell": shell,
            "Parameter": parameter,
            "Parameters": space.parameter_names,
            "Convention": "row=bra/destination, column=ket/source",
        },
    )


def _resolve_msg_id(group: Union[int, MSG, MagneticGroupSelector]) -> str:
    catalog = DataCatalog()
    if isinstance(group, bool):
        raise TypeError("group must be a stable MSG selector")
    if isinstance(group, int):
        return catalog.resolve_msg_source_index(group)
    if isinstance(group, MSG):
        return group.stable_id or catalog.resolve_msg_id(bns=group.bns)  # type: ignore[arg-type]
    if isinstance(group, MagneticGroupSelector):
        return catalog.resolve_msg_key(group.classification, group.source_key)
    raise TypeError("group must be an integer, MSG, or MagneticGroupSelector")


def _catalog_value(value: Any) -> Any:
    if isinstance(value, Mapping) and "kind" in value:
        return _catalog_exact_to_user(value)
    if isinstance(value, Mapping):
        return {key: _catalog_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return tuple(_catalog_value(item) for item in value)
    return value


def showMSGWyckoff(group: Union[int, MSG, MagneticGroupSelector]) -> TableReport:
    """Return the complete stable Wyckoff table for one magnetic space group."""

    catalog = DataCatalog()
    stable_id = _resolve_msg_id(group)
    msg = catalog.msg(stable_id)
    wyckoff = catalog.wyckoff(stable_id)
    rows = tuple(
        {
            "Multiplicity": int(entry["multiplicity"]),
            "WyckoffLetter": str(entry["letter"]),
            "SourceOrdinal": int(entry["source_ordinal"]),
            "AtomicPositionsAndMagnetizationDirections": tuple(
                {
                    "Position": _catalog_value(position["coordinates"]),
                    "Moment": _catalog_value(position["moment"]),
                }
                for position in entry["positions"]
            ),
        }
        for entry in reversed(wyckoff["entries"])
    )
    return TableReport(
        schema="MagneticTBMSGWyckoffTable",
        title=f"MSG: {tuple(msg['bns_key'])} {msg['symbol']}",
        columns=tuple(rows[0]) if rows else (),
        rows=rows,
        metadata={
            "SchemaVersion": 1,
            "StableId": stable_id,
            "BNS": tuple(msg["bns_key"]),
            "Symbol": msg["symbol"],
            "EntryCount": len(rows),
        },
    )


# PEP-8 aliases keep notebooks natural without changing the stable names.
orbital_table = orbitalTable
show_hamiltonian_basis = showHamiltonianBasis
show_bonds = showbonds
bond_table = bondTable
show_symmetry_representations = showSymmetryRepresentations
show_hopping_parameters = showHoppingParameters
show_msg_wyckoff = showMSGWyckoff


__all__ = [
    "TableReport",
    "orbitalTable",
    "orbital_table",
    "showHamiltonianBasis",
    "show_hamiltonian_basis",
    "showbonds",
    "show_bonds",
    "bondTable",
    "bond_table",
    "showSymmetryRepresentations",
    "show_symmetry_representations",
    "showHoppingParameters",
    "show_hopping_parameters",
    "showMSGWyckoff",
    "show_msg_wyckoff",
]
