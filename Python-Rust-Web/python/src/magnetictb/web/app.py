"""FastAPI assembly and compatibility exports for the MagneticTB workbench."""

from __future__ import annotations

import os
from typing import Any, Dict

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles

from .. import __version__
from ..errors import MagneticTBError
from .codec import jsonable
from .crystal_geometry_schemas import BrillouinZonePlotRequest, CrystalStructurePlotRequest
from .data_schemas import MsgopRequest, SubperiodicOperationRequest
from .model_examples import INIT_EXAMPLES
from .model_input_service import (
    _bravais_option,
    _friendly_init,
    _parse_lattice_scalar,
    _wyckoff_options,
)
from .model_schemas import (
    CombineRequest,
    EvaluateRequest,
    FriendlyInitRequest,
    HamiltonianBasisReportRequest,
    HamiltonianOptions,
    HoppingReportRequest,
    InitDirectRepresentationRequest,
    InitFieldRequest,
    InitLatticeRequest,
    InitOperationRequest,
    InitOrbitRequest,
    InitSymmetryRequest,
    PrepareRequest,
    ShellReportRequest,
    SymmetryReportRequest,
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
from .route_dependencies import STATIC_DIRECTORY
from .runtime import MAX_REQUEST_BYTES, _REGISTRY


def normalize_base_path(value: str) -> str:
    base_path = value.strip()
    if not base_path or base_path == "/":
        return ""
    if not base_path.startswith("/"):
        raise RuntimeError("MAGNETICTB_WEB_BASE_PATH must start with '/'")
    if "?" in base_path or "#" in base_path or "//" in base_path:
        raise RuntimeError("MAGNETICTB_WEB_BASE_PATH is not a valid URL path")
    segments = base_path.strip("/").split("/")
    if any(segment in {"", ".", ".."} for segment in segments):
        raise RuntimeError("MAGNETICTB_WEB_BASE_PATH is not a valid URL path")
    return "/" + "/".join(segments)


class BodyLimitMiddleware:
    """Reject oversized HTTP bodies before Pydantic or Rust sees them."""

    def __init__(self, app: Any, maximum: int) -> None:
        self.application = app
        self.maximum = maximum

    async def __call__(self, scope: Dict[str, Any], receive: Any, send: Any) -> None:
        if scope.get("type") != "http":
            await self.application(scope, receive, send)
            return
        headers = dict(scope.get("headers", []))
        raw_length = headers.get(b"content-length")
        if raw_length is not None:
            try:
                if int(raw_length) > self.maximum:
                    response = JSONResponse(
                        status_code=413,
                        content={"error": {"tag": "WebRequestTooLarge", "detail": "request body exceeds the configured limit"}},
                    )
                    await response(scope, receive, send)
                    return
            except ValueError:
                pass
        buffered = []
        received = 0
        while True:
            message = await receive()
            buffered.append(message)
            if message.get("type") != "http.request":
                break
            received += len(message.get("body", b""))
            if received > self.maximum:
                response = JSONResponse(
                    status_code=413,
                    content={"error": {"tag": "WebRequestTooLarge", "detail": "request body exceeds the configured limit"}},
                )
                await response(scope, receive, send)
                return
            if not message.get("more_body", False):
                break
        iterator = iter(buffered)

        async def replay_receive() -> Dict[str, Any]:
            try:
                return next(iterator)
            except StopIteration:
                return {"type": "http.request", "body": b"", "more_body": False}

        await self.application(scope, replay_receive, send)


BASE_PATH = normalize_base_path(os.environ.get("MAGNETICTB_WEB_BASE_PATH", ""))

app = FastAPI(
    title="MagneticTB Web",
    description="Local browser API backed by the MagneticTB Rust mathematical core",
    version=__version__,
    root_path=BASE_PATH,
    docs_url="/api/docs",
    redoc_url=None,
)
app.add_middleware(BodyLimitMiddleware, maximum=MAX_REQUEST_BYTES)
app.mount("/static", StaticFiles(directory=STATIC_DIRECTORY), name="static")

# Generated from the authoritative documentation by scripts/package-web-help.py.
# Never depend on a checkout (or expose its development scripts) at runtime.
HELP_DIRECTORY = STATIC_DIRECTORY / "help"
app.mount("/help", StaticFiles(directory=HELP_DIRECTORY, html=True), name="help")


@app.middleware("http")
async def security_headers(request: Request, call_next: Any) -> Any:
    response = await call_next(request)
    response.headers["X-Content-Type-Options"] = "nosniff"
    response.headers["X-Frame-Options"] = "DENY"
    response.headers["Referrer-Policy"] = "no-referrer"
    response.headers["Permissions-Policy"] = "camera=(), microphone=(), geolocation=()"
    response.headers["Cross-Origin-Resource-Policy"] = "same-origin"
    response.headers["Content-Security-Policy"] = (
        "default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; "
        "font-src 'self'; img-src 'self' data:; connect-src 'self'; object-src 'none'; "
        "frame-ancestors 'none'; base-uri 'self'; form-action 'self'"
    )
    return response


@app.exception_handler(MagneticTBError)
async def magnetic_error_handler(request: Request, error: MagneticTBError) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={
            "error": {
                "tag": error.tag,
                "detail": error.detail,
                "domain": type(error).__name__,
            }
        },
    )


@app.exception_handler(RequestValidationError)
async def validation_error_handler(request: Request, error: RequestValidationError) -> JSONResponse:
    return JSONResponse(
        status_code=422,
        content={"error": {"tag": "WebValidationError", "detail": jsonable(error.errors())}},
    )

from .base_routes import (
    router as base_routes_router,
    index,
    health,
    capabilities,
    init_options,
    init_examples,
    init_msg_family,
    init_wyckoff_options,
    init_bravais_option,
)

from .model_routes import (
    router as model_routes_router,
    create_init_model,
    list_models,
    create_model,
    get_model,
    delete_model,
    model_symmetry_operations,
    model_orbitals,
    model_orbital_table_report,
    model_hamiltonian_basis_report,
    model_bond_report,
    model_symmetry_report,
    model_hopping_report,
    model_representations,
    model_identity,
    model_canonical,
    model_bonds,
    model_symham,
    model_hamiltonian,
    model_unsymham,
    model_hamiltonian_space,
    evaluate_hamiltonian,
    combine_hamiltonians,
    model_hopping_data,
)

from .properties_routes import (
    router as properties_routes_router,
    model_property_parameter_names,
    model_bloch_hamiltonian,
    model_band_plot,
    model_real_space_hamiltonian,
    model_slab_hamiltonian,
    model_surface_green_function,
    model_wilson_loop,
    model_berry_phase,
    model_berry_curvature,
    model_point_chern,
    model_gapless_points,
    model_crystal_structure_plot,
    model_brillouin_zone_plot,
    model_wilson_loop_plot,
    model_surface_spectrum_plot,
    model_berry_curvature_2d_plot,
    model_berry_curvature_3d_plot,
    model_real_space_wavefunction_plot,
)

from .data_routes import (
    router as data_routes_router,
    data_msgop,
    data_subperiodic_operations,
    data_msg_wyckoff_report,
    data_summary,
    data_classifications,
    list_data,
    get_data_record,
)

from .core_routes import (
    router as core_routes_router,
    core_dispatch,
)

app.include_router(base_routes_router)
app.include_router(model_routes_router)
app.include_router(properties_routes_router)
app.include_router(data_routes_router)
app.include_router(core_routes_router)
