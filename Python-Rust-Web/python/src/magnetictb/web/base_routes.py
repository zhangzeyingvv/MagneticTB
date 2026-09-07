"""Base Routes for the local MagneticTB workbench."""

from __future__ import annotations

from fastapi import APIRouter

from .route_dependencies import (
    Any,
    COMPUTE_TIMEOUT,
    Dict,
    FileResponse,
    HTTPException,
    INIT_EXAMPLES,
    JSONResponse,
    MAX_CONCURRENT,
    MAX_REQUEST_BYTES,
    MAX_RESPONSE_BYTES,
    MODEL_LIMIT,
    STATIC_DIRECTORY,
    _REGISTRY,
    __version__,
    _bounded,
    _bravais_option,
    _compute,
    _msg_family_options,
    _wyckoff_options,
    partial,
)


router = APIRouter()


@router.api_route("/", methods=["GET", "HEAD"], include_in_schema=False)
async def index() -> FileResponse:
    return FileResponse(STATIC_DIRECTORY / "index.html")


@router.get("/api/health")
async def health() -> Dict[str, Any]:
    return {
        "status": "ok",
        "magnetictb_version": __version__,
        "rust_core": True,
        "katex": "bundled",
        "models": len(_REGISTRY.list()),
        "limits": {
            "request_bytes": MAX_REQUEST_BYTES,
            "response_bytes": MAX_RESPONSE_BYTES,
            "concurrent_computations": MAX_CONCURRENT,
            "response_timeout_seconds": COMPUTE_TIMEOUT,
            "model_capacity": MODEL_LIMIT,
        },
    }


@router.get("/api/capabilities")
async def capabilities() -> Dict[str, Any]:
    return {
        "modeling": [
            "init",
            "initfromrep",
            "CurrentModelSession",
            "CompileMagneticTBInput",
            "prepare_model",
            "prepare_model_from_rep",
        ],
        "model_results": [
            "symham",
            "hamiltonian_object",
            "hamiltonian_space",
            "evaluate",
            "combine_shells",
            "unsymham",
        ],
        "reports": [
            "orbitalTable",
            "showHamiltonianBasis",
            "showbonds",
            "showSymmetryRepresentations",
            "showHoppingParameters",
            "showMSGWyckoff",
        ],
        "utilities": ["texOutput"],
        "properties": [
            "showCrystalStructure",
            "showBrillouinZone",
            "standardKPath",
            "bandplot",
            "bandManipulate",
            "hoppingData",
            "buildBlochHamiltonian",
            "buildRealSpaceHamiltonian",
            "buildSlabHamiltonian",
            "surfaceGreenFunction",
            "plotSurfaceSpectrum",
            "wilsonLoop",
            "wLoop",
            "z2path",
            "plotWilsonLoop",
            "berryPhase",
            "berryph",
            "berryCurvature",
            "plotBerryCurvature2D",
            "plotBerryCurvature3D",
            "plotRealSpaceWavefunction",
            "pointChernNumber",
            "findGaplessPoints",
        ],
        "data": ["bravais", "msg", "wyckoff", "rod", "layer", "msgop", "mlgop", "mrgop"],
        "rust_dispatchers": [
            "cyclotomic_context",
            "null_space",
            "common_kernel",
            "linear_algebra",
            "representation",
            "geometry",
            "tight_binding",
            "group",
        ],
        "math_core": "Rust",
        "python_role": "typed binding, session registry, and serialization only",
    }


@router.get("/api/init/options")
async def init_options() -> Dict[str, Any]:
    """Return the finite UI choices for every ordinary ``init`` option."""

    return {
        "lattice_presets": [
            {"value": "msg_bravais", "label": "Selected MSG Bravais lattice"},
            {"value": "stable_default", "label": "Stable default (hexagonal)"},
            {"value": "simple_cubic", "label": "Simple cubic identity"},
            {"value": "custom", "label": "Custom 3 × 3 lattice"},
        ],
        "basis_groups": [
            {"label": "scalar s", "values": ["s"]},
            {"label": "scalar p", "values": ["px", "py", "pz", "px+ipy", "px-ipy"]},
            {"label": "scalar d", "values": ["dx2-y2", "dz2", "dxy", "dyz", "dxz"]},
            {"label": "spin s", "values": ["sup", "sdn"]},
            {
                "label": "spin p",
                "values": [
                    "pxup", "pxdn", "pyup", "pydn", "pzup", "pzdn",
                    "px+ipy up", "px+ipy dn", "px-ipy up", "px-ipy dn",
                ],
            },
            {
                "label": "spin d",
                "values": [
                    "dx2-y2up", "dx2-y2dn", "dz2up", "dz2dn", "dxyup", "dxydn",
                    "dyzup", "dyzdn", "dxzup", "dxzdn",
                ],
            },
            {"label": "catalog test witnesses", "values": ["ptest3", "ptest4"]},
        ],
        "representation_modes": ["DirectProduct", "Induced"],
        "exact_fields": ["crystallographic", "rational", "gaussian", "custom"],
        "space_group_numbers": list(range(1, 231)),
        "defaults": {
            "initial_bond_shells": 10,
            "generate_symmetry_group": False,
            "representation_mode": "DirectProduct",
            "exact_field": "crystallographic",
        },
    }


@router.get("/api/init/examples")
async def init_examples() -> Dict[str, Any]:
    """Return replayed examples in the same shape accepted by ordinary init."""

    return {"examples": INIT_EXAMPLES}


@router.get("/api/init/msg-family/{space_group_number}")
async def init_msg_family(
    space_group_number: int,
) -> JSONResponse:
    if space_group_number < 1 or space_group_number > 230:
        raise HTTPException(
            status_code=422,
            detail={"tag": "InvalidSpaceGroupNumber", "detail": "space-group number must be 1..230"},
        )
    return _bounded(
        {
            "space_group_number": space_group_number,
            "groups": await _compute(partial(_msg_family_options, space_group_number)),
        }
    )


@router.get("/api/init/msg/{msg_id}/wyckoff-options")
async def init_wyckoff_options(msg_id: str) -> JSONResponse:
    return _bounded({"msg_id": msg_id, "wyckoff": await _compute(partial(_wyckoff_options, msg_id))})


@router.get("/api/init/bravais/{stable_id}")
async def init_bravais_option(stable_id: str) -> JSONResponse:
    return _bounded(await _compute(partial(_bravais_option, stable_id)))
