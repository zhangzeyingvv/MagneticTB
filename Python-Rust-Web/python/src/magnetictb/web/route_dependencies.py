"""Explicit dependency surface shared by the domain HTTP routers."""

from __future__ import annotations

import io
from contextlib import redirect_stdout
from functools import partial
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

from fastapi import Body, HTTPException, Query
from fastapi.responses import FileResponse, JSONResponse

from .. import __version__
from ..abstract_group import GroupAlgebra
from ..core_bindings import (
    common_kernel,
    cyclotomic_context,
    geometry,
    linear_algebra,
    null_space,
    representation,
    tight_binding,
)
from ..data_api import DataCatalog
from ..data import MagneticGroupSelector, mlgop, msgop, mrgop
from ..model import PreparedModel
from ..plotting import bandplot, standardKPath
from ..properties import (
    berryCurvature,
    berryPhase,
    buildBlochHamiltonian,
    buildRealSpaceHamiltonian,
    buildSlabHamiltonian,
    findGaplessPoints,
    pointChernNumber,
    surfaceGreenFunction,
    wilsonLoop,
)
from ..reports import (
    bondTable,
    orbitalTable,
    showHamiltonianBasis,
    showHoppingParameters,
    showMSGWyckoff,
    showSymmetryRepresentations,
)
from ..tight_binding import Hamiltonian
from ..visualization import (
    plotBerryCurvature2D,
    plotBerryCurvature3D,
    plotRealSpaceWavefunction,
    plotSurfaceSpectrum,
    plotWilsonLoop,
    showBrillouinZone,
    showCrystalStructure,
)
from .codec import decode_tagged, guard_json_tree, jsonable
from .crystal_geometry_schemas import (
    BrillouinZonePlotRequest,
    CrystalStructurePlotRequest,
)
from .data_schemas import MsgopRequest, SubperiodicOperationRequest
from .model_examples import INIT_EXAMPLES
from .model_input_service import (
    _bravais_option,
    _friendly_init,
    _msg_family_options,
    _wyckoff_options,
)
from .model_schemas import (
    CombineRequest,
    EvaluateRequest,
    FriendlyInitRequest,
    HamiltonianBasisReportRequest,
    HamiltonianOptions,
    HoppingReportRequest,
    PrepareRequest,
    ShellReportRequest,
    SymmetryReportRequest,
)
from .model_service import (
    _cumulative_hamiltonian,
    _hamiltonian,
    _hamiltonian_payload,
    _model_summary,
    _prepare,
    _property_centers,
    _property_hamiltonian,
    _property_hopping_data,
    _property_parameters,
    _property_shells,
)
from .properties_schemas import (
    BandPlotRequest,
    BerryCurvature2DPlotRequest,
    BerryCurvature3DPlotRequest,
    BerryCurvatureRequest,
    BerryPhaseRequest,
    BlochRequest,
    GaplessRequest,
    PointChernRequest,
    RealSpaceRequest,
    SlabRequest,
    SolvedShellRequest,
    SurfaceRequest,
    SurfaceSpectrumRequest,
    WavefunctionPlotRequest,
    WilsonPlotRequest,
    WilsonRequest,
)
from .runtime import (
    COMPUTE_TIMEOUT,
    MAX_CONCURRENT,
    MAX_REQUEST_BYTES,
    MAX_RESPONSE_BYTES,
    MODEL_LIMIT,
    _REGISTRY,
    _bounded,
    _catalog,
    _compute,
)


STATIC_DIRECTORY = Path(__file__).with_name("static")
