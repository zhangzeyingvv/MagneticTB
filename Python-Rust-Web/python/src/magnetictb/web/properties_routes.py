"""Properties Routes for the local MagneticTB workbench."""

from __future__ import annotations

from fastapi import APIRouter

from .route_dependencies import (
    Any,
    BandPlotRequest,
    BerryCurvature2DPlotRequest,
    BerryCurvature3DPlotRequest,
    BerryCurvatureRequest,
    BerryPhaseRequest,
    BlochRequest,
    BrillouinZonePlotRequest,
    CrystalStructurePlotRequest,
    GaplessRequest,
    JSONResponse,
    PointChernRequest,
    RealSpaceRequest,
    SlabRequest,
    SolvedShellRequest,
    SurfaceRequest,
    SurfaceSpectrumRequest,
    WavefunctionPlotRequest,
    WilsonPlotRequest,
    WilsonRequest,
    _REGISTRY,
    _bounded,
    _compute,
    _property_centers,
    _property_hamiltonian,
    _property_hopping_data,
    _property_parameters,
    bandplot,
    berryCurvature,
    berryPhase,
    buildBlochHamiltonian,
    buildRealSpaceHamiltonian,
    buildSlabHamiltonian,
    decode_tagged,
    findGaplessPoints,
    jsonable,
    partial,
    plotBerryCurvature2D,
    plotBerryCurvature3D,
    plotRealSpaceWavefunction,
    plotSurfaceSpectrum,
    plotWilsonLoop,
    pointChernNumber,
    showBrillouinZone,
    showCrystalStructure,
    standardKPath,
    surfaceGreenFunction,
    wilsonLoop,
)


router = APIRouter()


@router.post("/api/models/{model_id}/properties/parameter-names")
async def model_property_parameter_names(
    model_id: str, payload: SolvedShellRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    result = await _compute(partial(_property_hamiltonian, model, payload))
    return _bounded(
        {
            "shells": list(result.shells),
            "parameter_names": list(result.parameter_names),
        }
    )


@router.post("/api/models/{model_id}/properties/bloch")
async def model_bloch_hamiltonian(
    model_id: str, payload: BlochRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        return buildBlochHamiltonian(
            _property_hopping_data(model, payload),
            decode_tagged(payload.momentum),
            Output="Data",
        )

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/properties/bands")
async def model_band_plot(
    model_id: str, payload: BandPlotRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        path = standardKPath(
            BravaisType=payload.bravais_type,
            Tolerance=payload.tolerance,
            model=model,
        )
        return bandplot(
            path,
            payload.npoint,
            _property_hamiltonian(model, payload),
            _property_parameters(payload),
        ).to_dict()

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/properties/real-space")
async def model_real_space_hamiltonian(
    model_id: str, payload: RealSpaceRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        return buildRealSpaceHamiltonian(
            _property_hopping_data(model, payload),
            decode_tagged(payload.geometry),
            BoundaryConditions=payload.boundary_conditions,
            Output="Data",
        )

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/properties/slab")
async def model_slab_hamiltonian(
    model_id: str, payload: SlabRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        return buildSlabHamiltonian(
            _property_hopping_data(model, payload),
            payload.size,
            decode_tagged(payload.momentum),
            PeriodicDirections=payload.periodic_directions,
            CellMatrix=("Automatic" if payload.cell_matrix is None else payload.cell_matrix),
            Output="Data",
        )

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/properties/surface")
async def model_surface_green_function(
    model_id: str, payload: SurfaceRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        return surfaceGreenFunction(
            _property_hopping_data(model, payload),
            decode_tagged(payload.momentum),
            decode_tagged(payload.energy),
            CellMatrix=("Automatic" if payload.cell_matrix is None else payload.cell_matrix),
            Broadening=payload.broadening,
            Tolerance=payload.tolerance,
            MaxIterations=payload.max_iterations,
            Surface=payload.surface,
            Output="Data",
        )

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/properties/wilson-loop")
async def model_wilson_loop(
    model_id: str, payload: WilsonRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        hopping = _property_hopping_data(model, payload)
        return wilsonLoop(
            _property_hamiltonian(model, payload),
            _property_centers(hopping),
            payload.occupied,
            decode_tagged(payload.start),
            decode_tagged(payload.end),
            parameters=_property_parameters(payload),
            PathSubdivisions=payload.path_subdivisions,
            Output="Data",
        )

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/properties/berry-phase")
async def model_berry_phase(
    model_id: str, payload: BerryPhaseRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        hopping = _property_hopping_data(model, payload)
        return berryPhase(
            _property_hamiltonian(model, payload),
            _property_centers(hopping),
            payload.occupied,
            decode_tagged(payload.path),
            parameters=_property_parameters(payload),
            Output="Data",
        )

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/properties/berry-curvature")
async def model_berry_curvature(
    model_id: str, payload: BerryCurvatureRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        hopping = _property_hopping_data(model, payload)
        return berryCurvature(
            _property_hamiltonian(model, payload),
            _property_centers(hopping),
            payload.occupied,
            decode_tagged(payload.point),
            parameters=_property_parameters(payload),
            Directions=payload.directions,
            StepSize=decode_tagged(payload.step_size),
            Output="Data",
        )

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/properties/point-chern")
async def model_point_chern(
    model_id: str, payload: PointChernRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        return pointChernNumber(
            _property_hamiltonian(model, payload),
            payload.occupied,
            decode_tagged(payload.point),
            decode_tagged(payload.radius),
            parameters=_property_parameters(payload),
            SurfaceSubdivisions=payload.surface_subdivisions,
            RequireGaplessCenter=payload.require_gapless_center,
            Output="Data",
        )

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/properties/gapless-points")
async def model_gapless_points(
    model_id: str, payload: GaplessRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        return findGaplessPoints(
            _property_hamiltonian(model, payload),
            payload.occupied,
            parameters=_property_parameters(payload),
            BrillouinZone=decode_tagged(payload.brillouin_zone),
            GridSize=decode_tagged(payload.grid_size),
            CandidateCount=payload.candidate_count,
            GapTolerance=payload.gap_tolerance,
            MergeTolerance=payload.merge_tolerance,
            MaxIterations=payload.max_iterations,
            RefinementMethod=payload.refinement_method,
            Output="Data",
        )

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/plots/crystal-structure")
async def model_crystal_structure_plot(
    model_id: str, payload: CrystalStructurePlotRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    result = await _compute(
        partial(
            showCrystalStructure,
            CellRange=payload.cell_range,
            MomentScale=("Automatic" if payload.moment_scale is None else payload.moment_scale),
            AtomRadius=("Automatic" if payload.atom_radius is None else payload.atom_radius),
            ShowAtomLabels=payload.show_atom_labels,
            model=model,
        )
    )
    return _bounded(jsonable(result))


@router.post("/api/models/{model_id}/plots/brillouin-zone")
async def model_brillouin_zone_plot(
    model_id: str, payload: BrillouinZonePlotRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model
    result = await _compute(
        partial(
            showBrillouinZone,
            KPath=decode_tagged(payload.k_path),
            ShowKPath=payload.show_k_path,
            TranslationRange=payload.translation_range,
            Tolerance=payload.tolerance,
            model=model,
        )
    )
    return _bounded(jsonable(result))


@router.post("/api/models/{model_id}/plots/wilson-loop")
async def model_wilson_loop_plot(
    model_id: str, payload: WilsonPlotRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        hopping = _property_hopping_data(model, payload)
        return plotWilsonLoop(
            _property_hamiltonian(model, payload),
            _property_centers(hopping),
            payload.occupied,
            decode_tagged(payload.start),
            decode_tagged(payload.end),
            decode_tagged(payload.parameter_path),
            parameters=_property_parameters(payload),
            ParameterSubdivisions=payload.parameter_subdivisions,
            LoopSubdivisions=payload.path_subdivisions,
        )

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/plots/surface-spectrum")
async def model_surface_spectrum_plot(
    model_id: str, payload: SurfaceSpectrumRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        return plotSurfaceSpectrum(
            _property_hopping_data(model, payload),
            decode_tagged(payload.momentum_path),
            decode_tagged(payload.energy_range),
            MomentumSubdivisions=payload.momentum_subdivisions,
            EnergyPoints=payload.energy_points,
            Broadening=payload.broadening,
        )

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/plots/berry-curvature-2d")
async def model_berry_curvature_2d_plot(
    model_id: str, payload: BerryCurvature2DPlotRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        hopping = _property_hopping_data(model, payload)
        return plotBerryCurvature2D(
            _property_hamiltonian(model, payload),
            _property_centers(hopping),
            payload.occupied,
            decode_tagged(payload.ranges),
            parameters=_property_parameters(payload),
            Directions=payload.directions,
            FixedCoordinates=decode_tagged(payload.fixed_coordinates),
            GridSize=decode_tagged(payload.grid_size),
            StepSize=decode_tagged(payload.step_size),
        )

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/plots/berry-curvature-3d")
async def model_berry_curvature_3d_plot(
    model_id: str, payload: BerryCurvature3DPlotRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        hopping = _property_hopping_data(model, payload)
        return plotBerryCurvature3D(
            _property_hamiltonian(model, payload),
            _property_centers(hopping),
            payload.occupied,
            decode_tagged(payload.ranges),
            parameters=_property_parameters(payload),
            GridSize=decode_tagged(payload.grid_size),
            StepSize=decode_tagged(payload.step_size),
            MagnitudeThreshold=payload.magnitude_threshold,
        )

    return _bounded(await _compute(run))


@router.post("/api/models/{model_id}/plots/real-space-wavefunction")
async def model_real_space_wavefunction_plot(
    model_id: str, payload: WavefunctionPlotRequest
) -> JSONResponse:
    model = _REGISTRY.get(model_id).model

    def run() -> Any:
        finite = buildRealSpaceHamiltonian(
            _property_hopping_data(model, payload),
            decode_tagged(payload.geometry),
            Output="Data",
        )
        return plotRealSpaceWavefunction(
            finite,
            decode_tagged(payload.state),
            Aggregation=payload.aggregation,
            Normalize=payload.normalize,
            WeightThreshold=payload.weight_threshold,
            PhaseColoring=payload.phase_coloring,
        )

    return _bounded(await _compute(run))
