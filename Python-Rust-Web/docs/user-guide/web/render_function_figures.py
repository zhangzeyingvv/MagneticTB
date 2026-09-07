"""Render selected help examples through their public result.savefig() methods.

Run from this worktree with its built extension on PYTHONPATH. Only this help
directory is written. No fixture replay, alternate renderer, or model formula
is used. Select one or more page keys on the command line for incremental work.
"""

from __future__ import annotations

import argparse
import ast
import contextlib
import hashlib
import io
import json
import platform
import re
import subprocess
from pathlib import Path

import bs4
import matplotlib
import magnetictb
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parent
ASSETS = ROOT / "assets" / "function-plots"
MANIFEST = ASSETS / "manifest.json"

# Block indices refer only to input code, not printed output. Each tuple after
# the scope is (last executed block, variable name, image name, view options).
PAGES = {
    "bandplot": ("reference/bandplot.html", "main", None, {
        1: [("result", "bandplot-basic", {})],
        2: [("automatic", "bandplot-automatic", {})],
        **{i: [("candidate", f"bandplot-option-{i + 1}", {})] for i in range(3, 8)},
    }),
    "showband": ("reference/showband.html", "main", None, {
        1: [("result", "showband-basic", {})],
        **{i: [("candidate", f"showband-option-{i + 1}", {})] for i in range(2, 7)},
    }),
    "crystal": ("reference/showCrystalStructure.html", "main", None, {
        **{i: [("crystal", f"crystal-magnetism-{i}", {})] for i in range(1, 5)},
        5: [("expanded", "crystal-cell-range", {})],
        **{i: [("result", f"crystal-option-{i}", {})] for i in range(7, 12)},
    }),
    "zone": ("reference/showBrillouinZone.html", "main", None, {
        1: [("bz", "zone-hexagonal", {"elev": 60, "azim": 30})],
        2: [("custom", "zone-custom-path", {"elev": 90, "azim": -90})],
        **{i: [("result", f"zone-option-{i + 1}", {"elev": 60, "azim": 30})] for i in range(3, 8)},
        8: [("bz", "zone-tetragonal", {})],
        9: [("bz", "zone-orthorhombic", {})],
    }),
    "graphene": ("tutorial/GettingStarted.html", "#inputs, #hamiltonian, #bands", [0, *range(2, 8)], {
        7: [("result", "graphene-getting-started", {})],
    }),
    "kagome-crystal": ("tutorial/CrystalAndKPaths.html", "#model, #crystal", None, {
        1: [("crystal", "kagome-crystal", {})],
        2: [("top", "kagome-top-view", {"elev": 90, "azim": -90})],
    }),
    "kagome": ("tutorial/CrystalAndKPaths.html", "main", None, {
        5: [("bz", "kagome-zone", {})],
        7: [("bands", "kagome-nearest-bands", {})],
    }),
    "kagome-induced": ("reference/init.html", "#kagome-induced", None, {
        2: [("bands", "kagome-induced-bands", {})],
    }),
    "continuous": ("tutorial/ContinuousSymmetry.html", "main", list(range(6)), {
        2: [("structure", "continuous-crystal", {})],
        5: [("soc_bands", "continuous-soc-bands", {}),
            ("no_soc_bands", "continuous-no-soc-bands", {})],
    }),
    "materials": ("tutorial/GeneralExamples.html", "#graphene", None, {
        3: [("bands", "graphene-material-bands", {})],
    }),
    "mos2": ("tutorial/GeneralExamples.html", "#mos2", None, {
        1: [("mo_s2_bands", "mos2-material-bands", {})],
    }),
    "banddata": ("reference/banddata.html", "main", [0, 3], {
        3: [("plot", "banddata-graphene", {})],
    }),
    "comparison": ("reference/compareBand.html", "#fit, #comparison", None, {
        2: [("full", "compareband-full", {}), ("low", "compareband-low", {})],
    }),
    # Render from the tutorial's saved-parameter alternative, not fittingTB().
    "fitting-nearest": ("tutorial/BandFitting.html", "#data article:first-of-type, #model, #redraw-nearest", None, {
        2: [("nearest_comparison", "bandfitting-nearest", {})],
    }),
    "fitting-energy": ("tutorial/BandFitting.html", "#data article:first-of-type, #example-reference, #model, #energy-window", None, {
        3: [("fermi_comparison", "bandfitting-energy-window", {})],
    }),
    "fitting-k": ("tutorial/BandFitting.html", "#data article:first-of-type, #example-reference, #model, #k-neighborhood", None, {
        3: [("k_comparison", "bandfitting-k-neighborhood", {})],
    }),
}

REFERENCE_DATA_PAGES = {"comparison", "fitting-nearest", "fitting-energy", "fitting-k"}


def blocks_for(soup, scope):
    return [pre for part in soup.select(scope)
            for pre in part.select('.code-block:not(.output) pre')]


def source_for(pre):
    snippet = BeautifulSoup(str(pre), "html.parser")
    for item in snippet.select('.locale-en, .figure-export'):
        item.decompose()
    return snippet.get_text().strip()


def sha(value):
    return hashlib.sha256(value).hexdigest()


def gamma_label_syntax(source):
    """Compare code while allowing only G-to-Gamma string-label changes."""
    class Labels(ast.NodeTransformer):
        def visit_Constant(self, node):
            if isinstance(node.value, str) and node.value in ("G", "\\Gamma", "Γ"):
                return ast.copy_location(ast.Constant(value="Γ"), node)
            return node
    return ast.dump(Labels().visit(ast.parse(source)))


def check_mos2_bands(result, hamiltonian):
    """Compare this one plot with the stored immutable reference, before export.

    Read only the three MoS2 input cells and its nine numeric LineBox series;
    never execute Mathematica, regenerate fixtures, or replay other models.
    """
    tag = "mathematica-stable-2.0.10"
    commit = "4661f13f1dd0335dbe5f54c420fc4e1b6e1c7c70"
    source_path = "Developer/Documentation/Sources/Tutorials/GeneralExamples.source.wl"
    expected_sha = "b7d3bdece699b0a07e6dd1348f0df81054f190ea21701f4526b3835ed191f711"
    resolved = subprocess.check_output(["git", "rev-parse", tag + "^{commit}"], cwd=ROOT, text=True).strip()
    if resolved != commit:
        raise RuntimeError("STOP: MoS2 reference tag changed")
    raw = subprocess.check_output(["git", "show", commit + ":" + source_path], cwd=ROOT)
    if sha(raw) != expected_sha:
        raise RuntimeError("STOP: MoS2 reference source hash changed")
    source = raw.decode()
    cells = list(re.finditer(r'Cell\[BoxData\[("(?:\\.|[^"\\])*")\], "Input"\]', source))
    inputs = [json.loads(cell.group(1)) for cell in cells]
    index = inputs.index("bandplot[moS2Path, 60, moS2Hamiltonian, moS2Parameters]")
    output = source[cells[index].end():cells[index + 1].start()]
    series = re.findall(r'LineBox\[(\{\{.*?\}\})\]', output)
    lines = [ast.literal_eval(value.replace("{", "[").replace("}", "]")) for value in series]
    if len(lines) != 9 or any(len(line) != 61 for line in lines):
        raise RuntimeError("STOP: unexpected MoS2 reference curve shape")
    expected = [[lines[3 * segment + band][point][1] for band in range(3)]
                for segment in range(3) for point in range(61)]
    actual = result.eigenvalues
    if len(actual) != 183 or any(len(row) != 3 for row in actual):
        raise RuntimeError("STOP: unexpected MoS2 Python curve shape")
    residual = max(abs(value - reference) for row, ref in zip(actual, expected)
                   for value, reference in zip(row, ref))
    if residual > 1e-10:
        raise RuntimeError(f"STOP: MoS2 stable band mismatch: max absolute error {residual:.17g} eV")
    if not hamiltonian.covariance_verified:
        raise RuntimeError("STOP: MoS2 Hamiltonian covariance is not verified")
    return {
        "status": "pass", "baseline_tag": tag, "baseline_commit": commit,
        "source_path": source_path, "source_sha256": expected_sha,
        "input_cells": inputs[index - 2:index + 1],
        "reference_series_sha256": sha(json.dumps(lines).encode()),
        "sample_count": 183, "band_count": 3, "absolute_tolerance_eV": 1e-10,
        "max_absolute_error_eV": residual, "covariance_verified": True,
        "actual_eigenvalues_eV": actual,
        "boundary_eigenvalues_eV": [actual[i] for i in (0, 60, 121, 182)],
        "scope": "Only this fixed parameter set and its 183-by-3 numeric bandplot output are compared. This does not prove exact Hamiltonian parameter-space equality for arbitrary parameters or identical raw symham parameter numbering. No other model or full gate was rerun.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Check unchanged code/runtime and image hashes; do not execute examples")
    parser.add_argument("--gamma-labels-only", action="store_true",
                        help="Allow only G-to-Gamma source changes and export only band figures containing G labels")
    parser.add_argument("--view-only", action="store_true",
                        help="Allow only camera/label changes; keep unchanged figures and skip later unaffected models")
    parser.add_argument("pages", nargs="*", choices=list(PAGES))
    args = parser.parse_args()
    ASSETS.mkdir(parents=True, exist_ok=True)
    previous = json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {"pages": {}}
    source_root = Path(magnetictb.__file__).parent
    # Hash runtime sources, not bundles. Unchanged image groups are reused.
    sources = {
        str(path.relative_to(source_root)): sha(path.read_bytes())
        for path in sorted(source_root.rglob('*.py'))
        if 'web' not in path.relative_to(source_root).parts
    }
    sources.update({path.name: sha(path.read_bytes()) for path in source_root.glob('_core*.so')})
    runtime = {"python": platform.python_version(), "matplotlib": matplotlib.__version__,
               "beautifulsoup": bs4.__version__, "sources": sources}
    for key in args.pages or PAGES:
        page, scope, selected, captures = PAGES[key]
        soup = BeautifulSoup((ROOT / page).read_text(), "html.parser")
        blocks = blocks_for(soup, scope)
        selected = list(range(len(blocks))) if selected is None else selected
        source = {i: source_for(blocks[i]) for i in selected}
        dependencies = {}
        if key in REFERENCE_DATA_PAGES:
            data_path = ROOT / "assets/examples/compareband-reference.json"
            dependencies[str(data_path.relative_to(ROOT))] = sha(data_path.read_bytes())
        fingerprint_items = [source, captures, runtime]
        if dependencies:
            fingerprint_items.append(dependencies)
        fingerprint = sha(json.dumps(fingerprint_items, sort_keys=True).encode())
        old = previous["pages"].get(key, {})
        old_source = old.get("executed_source", {})
        syntax_matches = set(map(str, source)) == set(old_source) and all(
            ast.dump(ast.parse(code)) == ast.dump(ast.parse(old_source[str(index)]))
            for index, code in source.items()
        )
        expected_captures = [(index, variable, name + '.svg', view)
                             for index, items in captures.items()
                             for variable, name, view in items]
        old_captures = [(f['block'], f['variable'], f['file'], f['view'])
                        for f in old.get('figures', [])]
        if old.get("input_sha256") == fingerprint or (
            syntax_matches and old.get('runtime') == runtime and old_captures == expected_captures
            and old.get('dependencies', {}) == dependencies
        ):
            for asset in old["figures"]:
                path = ASSETS / asset["file"]
                if not path.is_file() or sha(path.read_bytes()) != asset["sha256"]:
                    raise RuntimeError(f"STOP: changed or missing generated asset: {path}")
            print(f"REUSE {key}: {len(old['figures'])} verified images", flush=True)
            continue
        if args.check:
            raise RuntimeError(f"STOP: {key} code or runtime changed; no examples were rerun")
        if args.gamma_labels_only or args.view_only:
            if old.get("runtime") != runtime or set(map(str, source)) != set(old_source) or any(
                gamma_label_syntax(code) != gamma_label_syntax(old_source[str(index)])
                for index, code in source.items()
            ):
                raise RuntimeError(f"STOP: {key} has changes beyond display labels/views")
            for asset in old["figures"]:
                if sha((ASSETS / asset["file"]).read_bytes()) != asset["sha256"]:
                    raise RuntimeError(f"STOP: changed cached figure: {asset['file']}")
        execute_indices = selected
        if args.view_only:
            if [item[:3] for item in expected_captures] != [item[:3] for item in old_captures]:
                raise RuntimeError("STOP: view-only update cannot change capture identities")
            changed_blocks = [new[0] for new, prior in zip(expected_captures, old_captures)
                              if new[3] != prior[3]]
            if not changed_blocks:
                raise RuntimeError("STOP: no changed camera views selected")
            execute_indices = [index for index in selected if index <= max(changed_blocks)]
        print(f"RUN {key}: {page}", flush=True)
        namespace = {"__name__": "__help_example__"}
        figures = []
        outputs = dict(old.get("outputs", {})) if args.view_only else {}
        validation = None
        for index in execute_indices:
            stream = io.StringIO()
            data_directory = contextlib.chdir(ROOT / "assets/examples") if key in REFERENCE_DATA_PAGES else contextlib.nullcontext()
            with contextlib.redirect_stdout(stream), data_directory:
                exec(compile(source[index], f"{page}:input-block-{index}", "exec"), namespace)
            outputs[str(index)] = stream.getvalue()
            if key == "kagome-crystal":
                # This display-only update reuses the documented geometry and
                # never reaches the tutorial's Hamiltonian/band input blocks.
                expected = previous["pages"]["kagome"]["outputs"][str(index)]
                if outputs[str(index)] != expected:
                    raise RuntimeError("STOP: Kagome geometry/report output changed")
            for variable, name, view in captures.get(index, []):
                result = namespace[variable]
                if args.view_only:
                    cached = next(asset for asset in old["figures"] if asset["file"] == name + ".svg")
                    if cached["view"] == view:
                        figures.append(cached)
                        print(f"  KEEP {cached['file']}", flush=True)
                        continue
                if args.gamma_labels_only:
                    cached = next(asset for asset in old["figures"] if asset["file"] == name + ".svg")
                    if cached["result_type"] != "BandPlotResult" or b"<!-- G -->" not in (ASSETS / cached["file"]).read_bytes():
                        figures.append(cached)
                        print(f"  KEEP {cached['file']}", flush=True)
                        continue
                if key == "mos2":
                    validation = check_mos2_bands(result, namespace["mo_s2"])
                if key in REFERENCE_DATA_PAGES:
                    if len(result.model.eigenvalues) != 93 or len(result.reference_bands) != 90:
                        raise RuntimeError("STOP: unexpected saved-parameter comparison dimensions")
                    expected_samples = {
                        "fitting-nearest": ((0, (-0.41620768255750884, 2.4537921207343953, 2.4537921207343953)),),
                        "comparison": ((0, (-0.30885635496703756, 2.6028739005231274, 2.6028739005231274)),
                                       (60, (-0.10881749148798939, 0.15671209576878298, 2.419268293156152))),
                        "fitting-energy": ((60, (-0.13931254335422114, 0.13763489474515225, 1.5049173141544039)),),
                        "fitting-k": (),
                    }[key]
                    for sample_index, expected in expected_samples:
                        if max(abs(a - b) for a, b in zip(result.model.eigenvalues[sample_index], expected)) > 1e-10:
                            raise RuntimeError(f"STOP: saved-parameter comparison changed at model index {sample_index}")
                if key in ("fitting-energy", "fitting-k"):
                    fit = namespace["fermi_fit" if key == "fitting-energy" else "k_fit"]
                    expected_iterations, initial_loss, fitted_loss = (
                        (24, 16.772988964303096, 0.00017166800200888645)
                        if key == "fitting-energy" else
                        (20, 29.378270550141963, 0.0024843950327300217)
                    )
                    if not fit.converged or fit.iterations != expected_iterations:
                        raise RuntimeError("STOP: focused-fit convergence changed from the recorded example")
                    if abs(fit.initial_loss - initial_loss) > 1e-10 or abs(fit.fitted_loss - fitted_loss) > 1e-10:
                        raise RuntimeError("STOP: focused-fit loss changed from the recorded example")
                    if key == "fitting-energy" and (
                        len(fit.used_k_point_indices) != 28 or fit.residual_count != 41
                        or fit.energy_window != (-0.5, 0.5)
                    ):
                        raise RuntimeError("STOP: energy-window selection changed")
                    if key == "fitting-k":
                        if fit.used_k_point_indices != tuple(range(49, 71)):
                            raise RuntimeError("STOP: K-neighbourhood selection changed")
                        expected = (-0.5253262755254161, 0.6362529410871435, 2.278716163645225)
                        if max(abs(a - b) for a, b in zip(fit.comparison.fitted_bands[0], expected)) > 1e-10:
                            raise RuntimeError("STOP: K-neighbourhood fitted sample changed")
                    validation = {
                        "status": "pass", "absolute_tolerance": 1e-10,
                        "scope": "Only this targeted tutorial fit and its compareBand plot; checks against the previously recorded Python result, not a new cross-language or full-suite gate.",
                        "fit_result": fit.to_dict(),
                    }
                if not callable(getattr(result, "savefig", None)):
                    raise TypeError(f"{page}: {variable} has no public savefig()")
                path = ASSETS / f"{name}.svg"
                # The public API draws every artist and performs file export.
                result.savefig(path, **view, metadata={"Date": None})
                content = path.read_bytes()
                if b'<svg' not in content or b'Matplotlib' not in content:
                    raise RuntimeError(f"Invalid public figure export: {path}")
                figures.append({"block": index, "variable": variable, "file": path.name,
                                "result_type": type(result).__name__, "view": view,
                                "sha256": sha(content)})
                if key == "kagome-crystal":
                    figures[-1]["result_repr"] = repr(result)
                if key in REFERENCE_DATA_PAGES:
                    figures[-1]["comparison_data"] = result.to_dict()
                print(f"  SAVED {path.name} ({type(result).__name__})", flush=True)
        if args.view_only:
            updated = {figure["file"]: figure for figure in figures}
            figures = [updated.get(figure["file"], figure) for figure in old["figures"]]
            if outputs != old["outputs"]:
                raise RuntimeError("STOP: printed geometry changed during a view-only update")
        previous["pages"][key] = {"page": page, "scope": scope,
            "input_sha256": fingerprint, "figures": figures, "outputs": outputs,
            "executed_source": source, "runtime": runtime}
        if key == "kagome-crystal":
            previous["pages"][key]["validation_scope"] = (
                "Only one tutorial initialization and its two showCrystalStructure "
                "calls. Their printed geometry/order/counts match the historical "
                "record. Each image uses the public savefig method; no Hamiltonian, "
                "band, Brillouin-zone, fixture, or full-suite replay."
            )
        if dependencies:
            previous["pages"][key]["dependencies"] = dependencies
            previous["pages"][key]["validation_scope"] = (
                "First figure generation for this focused-fit example: one targeted fittingTB call, recorded-result checks, and compareBand.savefig. Complete fit parameters and comparison arrays are cached for reuse. No other fit or full cross-language gate was rerun."
                if key in ("fitting-energy", "fitting-k") else
                "Display and recorded sample-value checks using saved parameters; no fitting optimizer or full cross-language gate was rerun."
            )
        if validation is not None:
            previous["pages"][key]["validation"] = validation
            # Preserve the previous withheld result as history, not as an active gap.
            withheld = previous.get("withheld", {}).pop("MoS2", None)
            if withheld is not None:
                previous["pages"][key]["previous_withheld"] = withheld
        MANIFEST.write_text(json.dumps(previous, ensure_ascii=False, indent=2) + '\n')


if __name__ == "__main__":
    main()
