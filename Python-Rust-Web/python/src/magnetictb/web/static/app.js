(() => {
  "use strict";
  const i18n = window.MagneticTBI18n;
  const state = { models: [], modelId: "", hamiltonian: null, resultTab: "matrix", locale: i18n.start(), initOptions: null, initExamples: [], symmetrySource: "msg", wyckoffOptions: [], bravais: null, loadedDirectRepresentation: null, modelParameterDefaults: {}, modelBandShells: {}, modelParameterSources: {} };
  const $ = selector => document.querySelector(selector);
  const $$ = selector => [...document.querySelectorAll(selector)];
  const toast = (message, error = false) => {
    const node = $("#toast"); node.textContent = message; node.className = `toast show${error ? " error" : ""}`;
    clearTimeout(toast.timer); toast.timer = setTimeout(() => node.className = "toast", 3800);
  };
  const errorMessage = payload => {
    const detail = payload && (payload.error || payload.detail);
    if (detail && typeof detail === "object") { const tag = detail.tag || "Error"; let message = typeof detail.detail === "string" ? detail.detail : JSON.stringify(detail.detail); while (message.startsWith(`${tag}: `)) message = message.slice(tag.length + 2); return `${tag}: ${message}`; }
    return typeof detail === "string" ? detail : "Request failed";
  };
  const request = async (path, options = {}) => {
    const response = await fetch(path, { headers: { "Content-Type": "application/json", ...(options.headers || {}) }, ...options });
    let payload; try { payload = await response.json(); } catch (_) { payload = { detail: response.statusText }; }
    if (!response.ok) throw new Error(errorMessage(payload));
    return payload;
  };
  const json = value => JSON.stringify(value, null, 2);
  const escapeHtml = value => String(value).replace(/[&<>"]/g, character => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[character]);
  const parse = (value, label) => { try { return JSON.parse(value); } catch (error) { throw new Error(`${label} JSON: ${error.message}`); } };
  const withLoading = async (node, work) => { node.classList.add("loading"); try { return await work(); } finally { node.classList.remove("loading"); } };
  const typeset = root => {
    if (!window.katex) return;
    (root || document).querySelectorAll("[data-tex]").forEach(node => {
      const tex = node.dataset.tex || node.textContent;
      try { window.katex.render(tex, node, { throwOnError: false, displayMode: node.classList.contains("math-display") || node.classList.contains("ham-matrix"), strict: "ignore" }); }
      catch (_) { /* Accessible plain-text fallback stays visible. */ }
    });
  };
  const numericComplexText = value => { const real = Number(value.real), imaginary = Number(value.imaginary); if (!Number.isFinite(real) || !Number.isFinite(imaginary)) throw new Error("Rust numeric result contains NaN or Infinity"); return `${real.toPrecision(7)}${imaginary < 0 ? "" : "+"}${imaginary.toPrecision(7)}i`; };
  const scalarText = value => {
    if (value === null || value === undefined) return "0";
    if (typeof value === "object") return value.text || (value.type === "complex" ? numericComplexText(value) : JSON.stringify(value));
    return String(value);
  };
  const parameterTex = value => String(value).replace(/^([A-Za-z]+)(\d+)$/, "$1_{$2}");
  const rootOfUnityText = (orderValue, powerValue) => {
    const order = Number(orderValue), rawPower = Number(powerValue);
    const power = ((rawPower % order) + order) % order;
    if (power === 0) return "1";
    const gcd = (left, right) => { while (right) [left, right] = [right, left % right]; return Math.abs(left); };
    const divisor = gcd(2 * power, order), numerator = 2 * power / divisor, denominator = order / divisor;
    return `e^(${numerator === 1 ? "" : `${numerator}*`}pi*i${denominator === 1 ? "" : `/${denominator}`})`;
  };
  const rootExponentialTex = phase => {
    const match = String(phase).match(/^(-?)(?:(\d+)\*)?pi\*i(?:\/(\d+))?$/);
    if (!match) return `e^{${phase}}`;
    const [, sign, magnitude, denominator] = match;
    const coefficient = !magnitude || magnitude === "1" ? "" : magnitude;
    const numerator = `${sign}${coefficient}\\pi i`;
    return denominator ? `e^{\\frac{${numerator}}{${denominator}}}` : `e^{${numerator}}`;
  };
  const expressionText = text => String(text)
    .replace(/root_of_unity\((\d+),\s*(-?\d+)\)/g, (_, order, power) => rootOfUnityText(order, power));
  const exponentialTex = source => {
    const marker = "exp(I*("; let cursor = 0, output = "";
    while (true) {
      const start = source.indexOf(marker, cursor); if (start < 0) return output + source.slice(cursor);
      output += source.slice(cursor, start); let depth = 1, end = start + marker.length;
      for (; end < source.length && depth; end += 1) { if (source[end] === "(") depth += 1; else if (source[end] === ")") depth -= 1; }
      if (depth) return output + source.slice(start);
      const phase = source.slice(start + marker.length, end - 1); output += `e^{i\\left(${phase}\\right)}`;
      cursor = source[end] === ")" ? end + 1 : end;
    }
  };
  const expressionTex = text => {
    let value = expressionText(text).replaceAll("+ -", "- ");
    value = value.replace(/e\^\(([^()]*)\)/g, (_, phase) => rootExponentialTex(phase));
    value = exponentialTex(value);
    value = value.replace(/\b(kx|ky|kz)\b/g, symbol => `k_{${symbol[1]}}`);
    value = value.replace(/\b([A-Za-z]+)(\d+)\b/g, (_, name, index) => `${name}_{${index}}`);
    value = value.replace(/\bI\b/g, "i");
    value = value.replace(/sqrt\((\d+)\)/g, "\\sqrt{$1}");
    value = value.replace(/(?<![\w}])(\-?\d+)\/(\d+)/g, "\\frac{$1}{$2}");
    value = value.replaceAll("*", "\\,");
    return value;
  };
  const matrixTex = rows => `\\begin{pmatrix}${rows.map(row => row.map(expressionTex).join(" & ")).join(" \\\\ ")}\\end{pmatrix}`;
  const hamiltonianTex = value => `H(\\mathbf{k})=${value.tex || matrixTex(value.matrix_text)}`;
  const hoppingText = value => {
    const lines = ["# MagneticTB symbolic hopping table", "# Rust-generated exact real-space Fourier coefficients", "# d1 d2 d3  row column  parameter  coefficient"];
    (value.terms || []).forEach(term => (term.fourier_coefficients || []).forEach(record => (record.matrix_text || []).forEach((row, rowIndex) => row.forEach((coefficient, columnIndex) => { const text = expressionText(coefficient).replace(/\bI\b/g, "i"); if (text !== "0") lines.push(`${record.displacement_text.join(" ")}  ${rowIndex + 1} ${columnIndex + 1}  ${term.parameter}  ${text}`); }))));
    lines.push("# Parameters remain symbolic."); return lines.join("\n");
  };
  const downloadText = (name, content) => { const url = URL.createObjectURL(new Blob([content], { type: "text/plain;charset=utf-8" })); const link = document.createElement("a"); link.href = url; link.download = name; link.click(); setTimeout(() => URL.revokeObjectURL(url), 0); };
  const matrixFromNumeric = rows => `\\begin{pmatrix}${rows.map(row => row.map(cell => {
    if (cell && cell.type === "complex") return numericComplexText(cell);
    return scalarText(cell);
  }).join(" & ")).join(" \\\\ ")}\\end{pmatrix}`;
  const finiteNumber = (value, label) => { const number = Number(value); if (!String(value).trim() || !Number.isFinite(number)) throw new Error(`${label} 必须是有限数`); return number; };
  const vectorValues = (root, label) => [...root.querySelectorAll("input")].map((input, index) => finiteNumber(input.value, `${label}[${index + 1}]`));
  const commaVector = (value, label) => String(value).split(",").map((item, index) => finiteNumber(item.trim(), `${label}[${index + 1}]`));
  const amplitudeValues = value => String(value).split(",").map((item, index) => {
    const source = item.replaceAll(" ", "");
    if (!source) throw new Error(`振幅[${index + 1}] 不能为空`);
    if (!source.endsWith("i")) return finiteNumber(source, `振幅[${index + 1}]`);
    const core = source.slice(0, -1); let split = -1;
    for (let position = 1; position < core.length; position += 1) if ((core[position] === "+" || core[position] === "-") && !"eE".includes(core[position - 1])) split = position;
    const realSource = split < 0 ? "0" : core.slice(0, split);
    const imaginarySource = split < 0 ? core : core.slice(split);
    const imaginary = imaginarySource === "" || imaginarySource === "+" ? 1 : imaginarySource === "-" ? -1 : finiteNumber(imaginarySource, `振幅[${index + 1}] 虚部`);
    return { $type: "complex_number", real: finiteNumber(realSource, `振幅[${index + 1}] 实部`), imaginary };
  });
  const shellValues = () => { const values = $("#propertyShells").value.split(",").map(value => Number(value.trim())); if (!values.length || values.some(value => !Number.isInteger(value) || value < 1) || new Set(values).size !== values.length) throw new Error("壳层必须是无重复的正整数列表"); return values; };
  const defaultNumericParameter = name => /^e\d+$/.test(name) ? "0" : "1";
  const activeParameterDefaults = () => state.modelParameterDefaults[state.modelId] || {};
  const updateParameterHint = () => {
    const source = state.modelParameterSources[state.modelId];
    $("#propertyParameterHint").textContent = source
      ? `已载入 ${source} 的固定参数；可直接修改。`
      : "默认 onsite=0、hopping=1，仅用于避免全零平带预览；请按实际模型修改。";
  };
  const renderNumericParameters = (root, names, defaults = {}, preserve = true) => {
    const previous = {}; if (preserve) root.querySelectorAll("input[data-parameter]").forEach(input => { previous[input.dataset.parameter] = input.value; }); root.replaceChildren();
    if (!names.length) { const note = document.createElement("span"); note.className = "fixed-constraint"; note.textContent = "此结果没有自由参数"; root.appendChild(note); return; }
    names.forEach(name => { const label = document.createElement("label"); label.textContent = name; const input = document.createElement("input"); input.type = "number"; input.step = "any"; input.value = previous[name] ?? defaults[name] ?? defaultNumericParameter(name); input.dataset.parameter = name; label.appendChild(input); root.appendChild(label); });
  };
  const readNumericParameters = root => Object.fromEntries([...root.querySelectorAll("input[data-parameter]")].map(input => [input.dataset.parameter, finiteNumber(input.value, input.dataset.parameter)]));
  const propertyBase = () => ({ shells: shellValues(), parameters: readNumericParameters($("#propertyParameters")), hermitian: $("#hermitianInput").value === "true", validation_level: $("#validationInput").value, kernel_method: $("#kernelInput").value });
  const renderPropertyResult = (root, value) => {
    if (value.Hamiltonian) root.innerHTML = `<div class="result-meta"><span>${escapeHtml(value.Schema || "Hamiltonian")}</span><span>dimension ${value.Dimension || value.NumWannier || value.Hamiltonian.length}</span>${value.HermitianResidual === undefined ? "" : `<span>residual ${Number(value.HermitianResidual).toExponential(2)}</span>`}</div><div class="ham-matrix" data-tex="${escapeHtml(matrixFromNumeric(value.Hamiltonian))}"></div>`;
    else if (value.GreenFunction) root.innerHTML = `<div class="result-meta"><span>spectral weight ${Number(value.SpectralWeight).toPrecision(8)}</span><span>${value.Iterations} iterations</span><span>residual ${Number(value.CouplingResidual).toExponential(2)}</span></div><div class="ham-matrix" data-tex="${escapeHtml(matrixFromNumeric(value.GreenFunction))}"></div>`;
    else if (value.PhasesOverPi) root.innerHTML = `<div class="result-meta"><span>phases / π</span><span>${value.PhasesOverPi.map(item => Number(item).toPrecision(8)).join(" · ")}</span><span>min gap ${Number(value.MinimumGap).toExponential(2)}</span></div>`;
    else if (value.ChernNumber !== undefined) root.innerHTML = `<div class="result-meta"><span>Chern ${value.ChernNumber}</span><span>raw ${Number(value.RawChernNumber).toPrecision(8)}</span><span>error ${Number(value.QuantizationError).toExponential(2)}</span></div>`;
    else if (value.Curvature !== undefined) root.innerHTML = `<div class="result-meta"><span>curvature ${Number(value.Curvature).toPrecision(9)}</span><span>flux ${Number(value.Flux).toPrecision(9)}</span></div>`;
    else if (value.Points) root.innerHTML = `<div class="result-meta"><span>${value.Points.length} gapless point(s)</span><span>evaluations ${value.HamiltonianEvaluations}</span></div><pre class="code-result">${escapeHtml(value.Points.map((point, index) => `${index + 1}: (${point.map(item => Number(item).toPrecision(9)).join(", ")})`).join("\n"))}</pre>`;
    else root.innerHTML = `<pre class="code-result">${escapeHtml(json(value))}</pre>`;
    typeset(root);
  };
  const reportCell = value => {
    if (value === null || value === undefined) return "—";
    if (typeof value === "boolean") return value ? "True" : "False";
    if (typeof value !== "object") return String(value);
    if (Array.isArray(value)) return `[${value.map(reportCell).join(", ")}]`;
    if (typeof value.text === "string" && ["cyclotomic", "quadratic", "rational", "exact_expression"].includes(value.type)) return value.text.replaceAll("I", "i");
    return JSON.stringify(value);
  };
  const renderTableReport = (root, value) => {
    const columns = value.columns || [], rows = value.rows || [];
    root.className = "report-result";
    root.innerHTML = `<div class="result-meta"><span>${escapeHtml(value.title || value.schema)}</span><span>${rows.length} rows</span><span>${escapeHtml(value.schema || "report")}</span></div><div class="report-table-wrap"><table class="report-table"><thead><tr>${columns.map(column => `<th>${escapeHtml(column)}</th>`).join("")}</tr></thead><tbody>${rows.map(row => `<tr>${columns.map(column => `<td><code>${escapeHtml(reportCell(row[column]))}</code></td>`).join("")}</tr>`).join("")}</tbody></table></div>`;
  };
  const renderOperations = (root, value) => {
    root.innerHTML = `<p>${value.count} operations${value.og_key ? ` · OG ${value.og_key.join(".")}` : ""}</p><div class="report-table-wrap"><table class="report-table"><thead><tr><th>#</th><th>Label</th><th>Rotation</th><th>Translation</th><th>F/T</th></tr></thead><tbody>${value.operations.map((operation, index) => `<tr><td>${index + 1}</td><td>${escapeHtml(operation[0])}</td><td><code>${escapeHtml(JSON.stringify(operation[1]))}</code></td><td><code>${escapeHtml(JSON.stringify(operation[2]))}</code></td><td>${escapeHtml(operation[3])}</td></tr>`).join("")}</tbody></table></div>`;
  };
  const renderBandPlot = (root, value) => {
    const width = 920, height = 380, left = 62, right = 20, top = 20, bottom = 50;
    const energies = (value.bands || []).flat().map(Number).filter(Number.isFinite);
    if (!energies.length) throw new Error("能带结果没有有限本征值");
    let ymin = Math.min(...energies), ymax = Math.max(...energies);
    const span = Math.max(ymax - ymin, 1e-9), padding = Math.max(0.05 * span, 0.05);
    ymin -= padding; ymax += padding;
    const xmax = Math.max(...value.boundary_positions, 1);
    const x = item => left + Number(item) / xmax * (width - left - right);
    const y = item => top + (ymax - Number(item)) / (ymax - ymin) * (height - top - bottom);
    const segmentLength = Number(value.npoint) + 1;
    const paths = (value.bands || []).flatMap((band, bandIndex) =>
      (value.path || []).map((_, segmentIndex) => {
        const start = segmentIndex * segmentLength, stop = start + segmentLength;
        const points = band.slice(start, stop).map((energy, localIndex) => `${x(value.x_coordinates[start + localIndex]).toFixed(2)},${y(energy).toFixed(2)}`).join(" ");
        return `<polyline class="band-line band-${bandIndex % 8}" points="${points}"/>`;
      })
    ).join("");
    const vertical = value.boundary_positions.map(position => `<line class="band-grid" x1="${x(position)}" x2="${x(position)}" y1="${top}" y2="${height - bottom}"/>`).join("");
    const labels = value.boundary_positions.map((position, index) => `<text class="band-label" x="${x(position)}" y="${height - 18}" text-anchor="middle">${escapeHtml(value.boundary_labels[index] || "")}</text>`).join("");
    const zero = ymin <= 0 && ymax >= 0 ? `<line class="band-zero" x1="${left}" x2="${width - right}" y1="${y(0)}" y2="${y(0)}"/>` : "";
    root.className = "band-result";
    root.innerHTML = `<div class="result-meta"><span>${value.bands.length} bands</span><span>${value.path.length} segments</span><span>${value.eigenvalues.length} sampled points</span><span>${ymin.toPrecision(5)}…${ymax.toPrecision(5)}</span></div><svg class="band-svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="Band structure"><line class="band-axis" x1="${left}" x2="${left}" y1="${top}" y2="${height - bottom}"/><line class="band-axis" x1="${left}" x2="${width - right}" y1="${height - bottom}" y2="${height - bottom}"/>${vertical}${zero}${paths}${labels}<text class="band-y-label" transform="translate(18 ${height / 2}) rotate(-90)" text-anchor="middle">Energy</text><text class="band-y-tick" x="${left - 8}" y="${y(ymax) + 4}" text-anchor="end">${ymax.toPrecision(4)}</text><text class="band-y-tick" x="${left - 8}" y="${y(ymin) + 4}" text-anchor="end">${ymin.toPrecision(4)}</text></svg>`;
  };
  const renderBranchPlot = (root, value) => {
    const width = 760, height = 300, left = 52, right = 18, top = 18, bottom = 42;
    const rows = value.values || [], branchCount = rows[0]?.length || 0, coordinates = value.path_coordinate || [];
    if (!branchCount || !coordinates.length) throw new Error("Wilson 分支结果为空");
    const all = rows.flat().map(Number), ymin = Math.min(...all, -1), ymax = Math.max(...all, 1), xmin = Math.min(...coordinates), xmax = Math.max(...coordinates);
    const x = item => left + (Number(item) - xmin) / Math.max(xmax - xmin, 1e-12) * (width - left - right);
    const y = item => top + (ymax - Number(item)) / Math.max(ymax - ymin, 1e-12) * (height - top - bottom);
    const lines = Array.from({ length: branchCount }, (_, branch) => `<polyline class="band-line band-${branch % 8}" points="${rows.map((row, index) => `${x(coordinates[index]).toFixed(2)},${y(row[branch]).toFixed(2)}`).join(" ")}"/>`).join("");
    const boundaries = (value.boundary_coordinates || []).map((position, index) => `<line class="band-grid" x1="${x(position)}" x2="${x(position)}" y1="${top}" y2="${height - bottom}"/><text class="band-label" x="${x(position)}" y="${height - 14}" text-anchor="middle">${escapeHtml(value.boundary_labels?.[index] || "")}</text>`).join("");
    root.className = "plot-result"; root.innerHTML = `<div class="result-meta"><span>${branchCount} branches</span><span>${rows.length} samples</span><span>${escapeHtml(value.phase_convention)}</span></div><svg class="plot-svg" viewBox="0 0 ${width} ${height}">${boundaries}${lines}</svg>`;
  };
  const colorScale = (value, maximum) => { const ratio = Math.max(-1, Math.min(1, Number(value) / Math.max(maximum, 1e-15))); const strength = Math.round(55 + 165 * Math.abs(ratio)); return ratio >= 0 ? `rgb(${strength},${88 - Math.round(30 * ratio)},${112 - Math.round(40 * ratio)})` : `rgb(${88 - Math.round(30 * -ratio)},${130 + Math.round(45 * -ratio)},${strength})`; };
  const renderHeatmap = (root, xValues, yValues, rows, title, symmetric = false) => {
    const width = 760, height = 360, left = 56, right = 20, top = 20, bottom = 42, nx = xValues.length, ny = yValues.length;
    if (!nx || !ny || rows.length !== nx) throw new Error("网格结果为空或形状不一致");
    const flat = rows.flat().map(Number), maximum = Math.max(...flat.map(Math.abs), 1e-15), minimum = Math.min(...flat), upper = Math.max(...flat);
    const cellWidth = (width - left - right) / nx, cellHeight = (height - top - bottom) / ny;
    const cells = rows.map((row, ix) => row.map((value, iy) => `<rect x="${left + ix * cellWidth}" y="${top + (ny - iy - 1) * cellHeight}" width="${cellWidth + .35}" height="${cellHeight + .35}" fill="${colorScale(symmetric ? value : Number(value) - (minimum + upper) / 2, symmetric ? maximum : Math.max((upper - minimum) / 2, 1e-15))}"/>`).join("")).join("");
    root.className = "plot-result"; root.innerHTML = `<div class="result-meta"><span>${escapeHtml(title)}</span><span>${nx} × ${ny}</span><span>${minimum.toPrecision(5)}…${upper.toPrecision(5)}</span></div><svg class="plot-svg heatmap-svg" viewBox="0 0 ${width} ${height}">${cells}<line class="band-axis" x1="${left}" x2="${left}" y1="${top}" y2="${height - bottom}"/><line class="band-axis" x1="${left}" x2="${width - right}" y1="${height - bottom}" y2="${height - bottom}"/></svg>`;
  };
  const renderWavefunction = (root, value) => {
    const records = value.site_records || []; if (!records.length) throw new Error("波函数结果没有位置记录");
    const width = 760, height = 300, padding = 34, xs = records.map(record => Number(record.cartesian_position[0])), ys = records.map(record => Number(record.cartesian_position[1]));
    const xmin = Math.min(...xs), xmax = Math.max(...xs), ymin = Math.min(...ys), ymax = Math.max(...ys), span = Math.max(xmax - xmin, ymax - ymin, 1);
    const x = item => padding + (Number(item) - xmin) / span * (width - 2 * padding), y = item => height - padding - (Number(item) - ymin) / span * (height - 2 * padding);
    const circles = records.map(record => { const radius = 4 + 26 * Math.cbrt(record.weight / Math.max(value.maximum_weight, 1e-15)); const hue = ((Number(record.phase) / (2 * Math.PI)) % 1 + 1) % 1; return `<circle cx="${x(record.cartesian_position[0])}" cy="${y(record.cartesian_position[1])}" r="${radius}" fill="${value.phase_coloring ? `hsl(${hue * 360} 78% 55%)` : "#5969e8"}" fill-opacity="${record.weight > 0 ? .82 : .18}"/><text x="${x(record.cartesian_position[0])}" y="${y(record.cartesian_position[1]) + 4}" text-anchor="middle">${record.site_index}</text>`; }).join("");
    root.className = "plot-result"; root.innerHTML = `<div class="result-meta"><span>${records.length} sites</span><span>${value.visible_site_count} visible</span><span>norm ${Number(value.input_norm).toPrecision(6)}</span><span>${escapeHtml(value.aggregation)}</span></div><svg class="plot-svg wavefunction-svg" viewBox="0 0 ${width} ${height}">${circles}</svg>`;
  };
  const mount3dViewer = (root, points, draw, ariaLabel) => {
    const canvas = document.createElement("canvas");
    canvas.className = "interactive-3d"; canvas.setAttribute("role", "img"); canvas.setAttribute("aria-label", ariaLabel); canvas.tabIndex = 0;
    const hint = document.createElement("div"); hint.className = "viewer-hint"; hint.textContent = "拖动旋转 · 滚轮缩放 · 双击复位";
    root.append(canvas, hint);
    const numeric = points.map(point => point.map(Number)).filter(point => point.length === 3 && point.every(Number.isFinite));
    const center = [0, 1, 2].map(axis => numeric.reduce((sum, point) => sum + point[axis], 0) / Math.max(numeric.length, 1));
    const radius = Math.max(...numeric.map(point => Math.hypot(...point.map((value, axis) => value - center[axis]))), 1e-9);
    const initial = { yaw: -.72, pitch: .48, zoom: 1 }, view = { ...initial };
    let dragging = false, previous = [0, 0];
    const render = () => {
      const bounds = canvas.getBoundingClientRect(), width = Math.max(bounds.width, 320), height = Math.max(bounds.height, 360), ratio = Math.min(window.devicePixelRatio || 1, 2);
      canvas.width = Math.round(width * ratio); canvas.height = Math.round(height * ratio);
      const context = canvas.getContext("2d"); context.setTransform(ratio, 0, 0, ratio, 0, 0);
      const background = context.createLinearGradient(0, 0, 0, height); background.addColorStop(0, "#fbfcff"); background.addColorStop(1, "#f2f5fb"); context.fillStyle = background; context.fillRect(0, 0, width, height);
      const rotate = point => {
        const x = (Number(point[0]) - center[0]) / radius, y = (Number(point[1]) - center[1]) / radius, z = (Number(point[2]) - center[2]) / radius;
        const cy = Math.cos(view.yaw), sy = Math.sin(view.yaw), cp = Math.cos(view.pitch), sp = Math.sin(view.pitch);
        const x1 = cy * x - sy * z, z1 = sy * x + cy * z;
        return [x1, cp * y - sp * z1, sp * y + cp * z1];
      };
      const project = point => { const rotated = rotate(point), scale = Math.min(width, height) * .38 * view.zoom; return [width / 2 + rotated[0] * scale, height / 2 - rotated[1] * scale, rotated[2], 1]; };
      draw(context, project, rotate, width, height);
    };
    canvas.addEventListener("pointerdown", event => { dragging = true; previous = [event.clientX, event.clientY]; canvas.setPointerCapture(event.pointerId); });
    canvas.addEventListener("pointermove", event => { if (!dragging) return; view.yaw += (event.clientX - previous[0]) * .009; view.pitch = Math.max(-1.5, Math.min(1.5, view.pitch + (event.clientY - previous[1]) * .009)); previous = [event.clientX, event.clientY]; render(); });
    canvas.addEventListener("pointerup", event => { dragging = false; canvas.releasePointerCapture(event.pointerId); });
    canvas.addEventListener("pointercancel", () => { dragging = false; });
    canvas.addEventListener("wheel", event => { event.preventDefault(); view.zoom = Math.max(.35, Math.min(5, view.zoom * Math.exp(-event.deltaY * .001))); render(); }, { passive: false });
    canvas.addEventListener("dblclick", () => { Object.assign(view, initial); render(); });
    new ResizeObserver(render).observe(canvas); render();
  };
  const drawSegment3d = (context, project, left, right, color, width = 1.4, dashed = false) => {
    const a = project(left), b = project(right); context.beginPath(); context.moveTo(a[0], a[1]); context.lineTo(b[0], b[1]); context.strokeStyle = color; context.lineWidth = width; context.setLineDash(dashed ? [5, 5] : []); context.stroke(); context.setLineDash([]);
    return [a, b];
  };
  const drawArrow3d = (context, project, left, right, color, label = "") => {
    const [a, b] = drawSegment3d(context, project, left, right, color, 2.6); const angle = Math.atan2(b[1] - a[1], b[0] - a[0]);
    context.beginPath(); context.moveTo(b[0], b[1]); context.lineTo(b[0] - 9 * Math.cos(angle - .45), b[1] - 9 * Math.sin(angle - .45)); context.lineTo(b[0] - 9 * Math.cos(angle + .45), b[1] - 9 * Math.sin(angle + .45)); context.closePath(); context.fillStyle = color; context.fill();
    if (label) { context.fillStyle = color; context.font = "700 13px Inter, sans-serif"; context.fillText(label, b[0] + 7, b[1] - 7); }
  };
  const primitiveCellEdges = lattice => {
    const cartesian = fractional => [0, 1, 2].map(column => fractional.reduce((sum, value, row) => sum + value * Number(lattice[row][column]), 0));
    const edges = [];
    for (let axis = 0; axis < 3; axis += 1) for (let first = 0; first <= 1; first += 1) for (let second = 0; second <= 1; second += 1) { const corner = [0, 0, 0], other = [0, 1, 2].filter(value => value !== axis); corner[other[0]] = first; corner[other[1]] = second; const end = [...corner]; end[axis] = 1; edges.push([cartesian(corner), cartesian(end)]); }
    return edges;
  };
  const renderCrystalStructure = (root, value) => {
    const records = value.atom_records || [], edges = value.cell_edges || [];
    if (!records.length) throw new Error("晶体结构没有原子记录");
    const lattice = value.lattice, primitive = primitiveCellEdges(lattice), origin = [0, 0, 0];
    const points = [...records.map(record => record.cartesian_position), ...edges.flat(), ...primitive.flat(), ...records.filter(record => record.arrow_end).map(record => record.arrow_end), origin, ...lattice];
    root.className = "plot-result"; root.innerHTML = `<div class="result-meta"><span>${value.atom_count} atoms</span><span>${value.magnetic_atom_count} magnetic</span><span>${value.translations.length} cells</span><span>加粗框 = 原胞</span></div>`;
    mount3dViewer(root, points, (context, project) => {
      edges.forEach(edge => drawSegment3d(context, project, edge[0], edge[1], "rgba(111,122,148,.28)", 1));
      primitive.forEach(edge => drawSegment3d(context, project, edge[0], edge[1], "#35415e", 2.2));
      ["#d34c5b", "#288c78", "#4d68c8"].forEach((color, index) => drawArrow3d(context, project, origin, lattice[index], color, `a${index + 1}`));
      records.filter(record => record.arrow_end).forEach(record => drawArrow3d(context, project, record.cartesian_position, record.arrow_end, "#d9475b"));
      records.map(record => ({ record, point: project(record.cartesian_position) })).sort((left, right) => left.point[2] - right.point[2]).forEach(({ record, point }) => { const radius = Math.max(5, 8 * point[3]), hue = (Number(record.orbit_index) * 67) % 360; context.beginPath(); context.arc(point[0], point[1], radius, 0, 2 * Math.PI); context.fillStyle = `hsl(${hue} 68% 55%)`; context.fill(); context.strokeStyle = "rgba(255,255,255,.95)"; context.lineWidth = 1.5; context.stroke(); if (value.show_atom_labels) { context.fillStyle = "#39445f"; context.font = "600 12px Inter, sans-serif"; context.textAlign = "center"; context.fillText(`${record.orbit_index}.${record.equivalent_index}`, point[0], point[1] - radius - 5); } });
    }, "Interactive three-dimensional crystal structure and primitive cell");
  };
  const renderBrillouinZone = (root, value) => {
    const vertices = value.vertices || [], facets = value.facets || [];
    if (vertices.length < 4 || !facets.length) throw new Error("第一布里渊区结果为空");
    const pathSegments = value.cartesian_k_path || [], pathPoints = pathSegments.flat(), reciprocal = value.reciprocal_lattice, origin = [0, 0, 0];
    const labels = value.show_k_path ? (value.displayed_k_path || []).flatMap(segment => [[segment[1][0], segment[0][0]], [segment[1][1], segment[0][1]]]).filter((entry, index, all) => all.findIndex(candidate => candidate[0] === entry[0] && JSON.stringify(candidate[1]) === JSON.stringify(entry[1])) === index).map(entry => ({ label: entry[0].replace("\\Gamma", "Γ"), point: [0, 1, 2].map(column => entry[1].reduce((sum, item, row) => sum + Number(item) * Number(reciprocal[row][column]), 0)) })) : [];
    const edgeKeys = new Set(), zoneEdges = []; facets.forEach(facet => facet.forEach((left, index) => { const right = facet[(index + 1) % facet.length], key = [left, right].sort((a, b) => a - b).join(":"); if (!edgeKeys.has(key)) { edgeKeys.add(key); zoneEdges.push([vertices[left], vertices[right]]); } }));
    root.className = "plot-result"; root.innerHTML = `<div class="result-meta"><span>${escapeHtml(value.bz_type || "BZ")}</span><span>${vertices.length} vertices</span><span>${value.facet_count} facets</span><span>${pathSegments.length} path segments</span></div>`;
    mount3dViewer(root, [...vertices, ...pathPoints, origin, ...reciprocal], (context, project) => {
      facets.map((facet, index) => ({ index, points: facet.map(vertex => project(vertices[vertex])) })).sort((left, right) => left.points.reduce((sum, point) => sum + point[2], 0) - right.points.reduce((sum, point) => sum + point[2], 0)).forEach(facet => { context.beginPath(); facet.points.forEach((point, index) => index ? context.lineTo(point[0], point[1]) : context.moveTo(point[0], point[1])); context.closePath(); context.fillStyle = `hsla(${205 + facet.index * 7},72%,65%,.17)`; context.fill(); });
      zoneEdges.forEach(edge => drawSegment3d(context, project, edge[0], edge[1], "#4e668f", 1.45));
      ["#d34c5b", "#288c78", "#4d68c8"].forEach((color, index) => drawArrow3d(context, project, origin, reciprocal[index], color, `b${index + 1}`));
      if (value.show_k_path) pathSegments.forEach(segment => drawSegment3d(context, project, segment[0], segment[1], "#d9475b", 3));
      labels.forEach(record => { const point = project(record.point); context.beginPath(); context.arc(point[0], point[1], 4, 0, 2 * Math.PI); context.fillStyle = "#d9475b"; context.fill(); context.fillStyle = "#39445f"; context.font = "600 12px Inter, sans-serif"; context.fillText(record.label, point[0] + 7, point[1] - 7); });
    }, "Interactive three-dimensional first Brillouin zone");
  };
  const summary = model => `
    <div class="summary-grid">
      <div class="summary-item"><small>群阶</small><strong>${model.group_order}</strong></div>
      <div class="summary-item"><small>操作 / 轨道</small><strong>${model.operation_count} / ${model.orbital_count}</strong></div>
      <div class="summary-item"><small>准备键层</small><strong>${model.initial_bond_shells}</strong></div>
      <div class="summary-item"><small>逐壳键数</small><strong>${model.bond_shell_counts.join(" · ")}</strong></div>
    </div>`;
  const refreshModels = async preferred => {
    const payload = await request("api/models"); state.models = payload.models || [];
    if (preferred) state.modelId = preferred;
    if (!state.models.some(model => model.model_id === state.modelId)) state.modelId = state.models[0]?.model_id || "";
    const picker = $("#modelPicker"); picker.replaceChildren();
    if (!state.models.length) { const option = new Option("尚无模型", ""); picker.add(option); }
    state.models.forEach(model => picker.add(new Option(`${model.label} · ${model.model_id.slice(0, 12)}`, model.model_id)));
    picker.value = state.modelId; $("#modelCount").textContent = String(state.models.length);
    const current = state.models.find(model => model.model_id === state.modelId);
    $("#modelResult").className = current ? "" : "empty-state";
    $("#modelResult").innerHTML = current ? summary(current) : "<span>◇</span><p>创建模型后，这里显示群阶、轨道和键层。</p>";
    if (current && preferred) $("#modelSummaryCard").open = true;
  };
  const createModel = async (entrypoint, argumentsValue, label) => {
    const result = await request("api/models", { method: "POST", body: JSON.stringify({ entrypoint, arguments: argumentsValue, label }) });
    await refreshModels(result.model_id); toast(`模型已准备：${result.label}`); return result;
  };
  const renderHamiltonian = () => {
    const root = $("#hamiltonianResult"), value = state.hamiltonian;
    $$("[data-result-tab]").forEach(button => button.classList.toggle("active", button.dataset.resultTab === state.resultTab));
    if (!value) return;
    const meta = `<div class="result-meta"><span>shape ${value.shape.join(" × ")}</span><span>shell ${value.shells.join(" + ")}</span><span>${value.parameter_names.map(parameterTex).join(", ") || "no parameters"}</span><span>covariance ${value.covariance_verified}</span></div>`;
    if (state.resultTab === "matrix") root.innerHTML = `${meta}<div class="ham-matrix" data-tex="${escapeHtml(matrixTex(value.matrix_text))}"></div>`;
    else if (state.resultTab === "hopping") { const content = hoppingText(value); root.innerHTML = `${meta}<div class="export-toolbar"><strong>Symbolic hopping</strong><button class="mini-button" type="button">下载 hopping.dat</button></div><pre class="code-result">${escapeHtml(content)}</pre>`; root.querySelector("button").onclick = () => downloadText("hopping.dat", content); }
    else { const content = hamiltonianTex(value); root.innerHTML = `${meta}<div class="export-toolbar"><strong>TeX</strong><button class="mini-button" type="button">下载 Hamiltonian.tex</button></div><pre class="code-result">${escapeHtml(content)}</pre>`; root.querySelector("button").onclick = () => downloadText("Hamiltonian.tex", content); }
    root.className = "result-surface"; typeset(root);
  };
  const hamiltonianOptions = () => ({ shell: Number($("#shellInput").value), hermitian: $("#hermitianInput").value === "true", cartesian_coordinates: false, validation_level: $("#validationInput").value, kernel_method: $("#kernelInput").value });
  const requireModel = () => { if (!state.modelId) throw new Error("请先建立并选择模型"); return state.modelId; };
  const groupSelector = () => {
    const kind = $("#msgopKind").value;
    const values = $("#msgopValue").value.split(",").map(item => Number(item.trim()));
    if (!values.length || values.some(value => !Number.isInteger(value) || value < 1)) throw new Error("MSG 选择器必须是正整数或逗号分隔的正整数 key");
    if (kind === "source_index") { if (values.length !== 1) throw new Error("source index 只接受一个整数"); return { kind, value: values[0] }; }
    if (kind === "bns") { if (values.length !== 2) throw new Error("BNS key 必须是 n,m"); return { kind: "msg", bns: values }; }
    const expected = kind === "gray" ? 1 : 2;
    if (values.length !== expected) throw new Error(`${kind} 需要 ${expected} 个整数`);
    return { kind: "selector", classification: kind, source_key: values };
  };

  const makeMatrixEditor = (root, values, className = "matrix-value") => {
    root.replaceChildren();
    values.flat().forEach((value, index) => { const input = document.createElement("input"); input.value = value; input.className = className; input.dataset.row = String(Math.floor(index / 3)); input.dataset.column = String(index % 3); input.setAttribute("aria-label", `matrix ${Math.floor(index / 3) + 1},${index % 3 + 1}`); root.appendChild(input); });
  };
  const readMatrix = (root, selector = ".matrix-value") => {
    const values = [...root.querySelectorAll(selector)].map(input => input.value.trim());
    if (values.length !== 9) throw new Error("矩阵必须包含 3 × 3 个标量");
    return [values.slice(0, 3), values.slice(3, 6), values.slice(6, 9)];
  };
  const makeVector = (className, values = ["0", "0", "0"]) => {
    const root = document.createElement("div"); root.className = "vector-editor";
    values.forEach(value => { const input = document.createElement("input"); input.className = className; input.value = value; root.appendChild(input); }); return root;
  };
  const readVector = (root, className) => [...root.querySelectorAll(`.${className}`)].map(input => input.value.trim());
  const addLatticeParameter = (name = "", value = "") => {
    const row = document.createElement("div"); row.className = "parameter-row";
    const key = document.createElement("input"); key.className = "parameter-name"; key.placeholder = "a"; key.value = name;
    const scalar = document.createElement("input"); scalar.className = "parameter-value"; scalar.placeholder = "1 或 sqrt(3)/2"; scalar.value = value;
    const remove = document.createElement("button"); remove.type = "button"; remove.textContent = "×"; remove.onclick = () => row.remove();
    row.append(key, scalar, remove); $("#latticeParameters").appendChild(row);
  };
  const renderBravaisParameters = value => {
    const root = $("#latticeParameters"); root.replaceChildren();
    value.parameter_symbols.forEach(name => { const label = document.createElement("label"); label.textContent = name; const input = document.createElement("input"); input.className = "bravais-parameter-value"; input.dataset.name = name; input.value = value.parameter_defaults[name]; input.placeholder = name; label.appendChild(input); root.appendChild(label); });
  };
  const renderExplicitLatticeParameters = () => { const root = $("#latticeParameters"); root.replaceChildren(); addLatticeParameter("a", "1"); addLatticeParameter("c", "3"); };
  const renderBravaisPreview = value => {
    const root = $("#bravaisPreview"); root.replaceChildren(); const title = document.createElement("strong"); title.textContent = `${value.stable_id} · primitive vectors`; root.appendChild(title); const matrix = document.createElement("div"); matrix.className = "bravais-matrix";
    value.primitive_vectors.forEach(row => { const line = document.createElement("code"); row.forEach(cell => { const item = document.createElement("span"); item.textContent = cell; line.appendChild(item); }); matrix.appendChild(line); }); root.appendChild(matrix); const note = document.createElement("div"); note.textContent = `自由参数：${value.parameter_symbols.join(", ") || "无"}`; root.appendChild(note);
  };
  const basisPicker = selectedValues => {
    const details = document.createElement("details"); details.className = "basis-picker"; details.open = true;
    details._basisOrder = [...selectedValues];
    const summary = document.createElement("summary"); summary.textContent = "选择有序局域 basis functions";
    const groups = document.createElement("div"); groups.className = "basis-groups";
    (state.initOptions?.basis_groups || []).forEach(group => {
      const section = document.createElement("div"); section.className = "basis-group"; const title = document.createElement("strong"); title.textContent = group.label;
      const chips = document.createElement("div"); chips.className = "basis-chips";
      group.values.forEach(value => { const label = document.createElement("label"); label.className = "basis-chip"; const input = document.createElement("input"); input.type = "checkbox"; input.value = value; input.checked = selectedValues.includes(value); input.onchange = () => { if (input.checked && !details._basisOrder.includes(value)) details._basisOrder.push(value); if (!input.checked) details._basisOrder = details._basisOrder.filter(item => item !== value); updateBasisSummary(details); }; const text = document.createElement("span"); text.textContent = value; label.append(input, text); chips.appendChild(label); });
      section.append(title, chips); groups.appendChild(section);
    });
    const selected = document.createElement("div"); selected.className = "basis-summary"; details.append(summary, groups, selected); updateBasisSummary(details); return details;
  };
  const updateBasisSummary = details => {
    const summary = details.querySelector(".basis-summary"); if (!summary) return; summary.replaceChildren();
    if (!details._basisOrder.length) { summary.textContent = "至少选择一个 basis"; return; }
    const label = document.createElement("span"); label.textContent = "ordered:"; summary.appendChild(label);
    details._basisOrder.forEach((value, index) => { const item = document.createElement("span"); item.className = "ordered-basis-item"; const text = document.createElement("code"); text.textContent = value; const up = document.createElement("button"); up.type = "button"; up.textContent = "↑"; up.disabled = index === 0; up.onclick = () => { [details._basisOrder[index - 1], details._basisOrder[index]] = [details._basisOrder[index], details._basisOrder[index - 1]]; updateBasisSummary(details); }; const down = document.createElement("button"); down.type = "button"; down.textContent = "↓"; down.disabled = index === details._basisOrder.length - 1; down.onclick = () => { [details._basisOrder[index], details._basisOrder[index + 1]] = [details._basisOrder[index + 1], details._basisOrder[index]]; updateBasisSummary(details); }; item.append(text, up, down); summary.appendChild(item); });
  };
  const wyckoffKey = item => `${item.source_ordinal}:${item.letter}`;
  const selectedWyckoff = card => state.wyckoffOptions.find(item => wyckoffKey(item) === card.querySelector(".orbit-wyckoff").value);
  const fillWyckoffSelect = (select, preferred) => {
    const previous = preferred || select.value; select.replaceChildren();
    if (state.symmetrySource !== "msg") { select.add(new Option("显式 site seed", "explicit")); return; }
    state.wyckoffOptions.forEach(item => select.add(new Option(`${item.letter} · multiplicity ${item.multiplicity} · ordinal ${item.source_ordinal}`, wyckoffKey(item))));
    if ([...select.options].some(option => option.value === String(previous))) select.value = String(previous);
  };
  const renderParameterBindings = (root, names, className, defaults = {}) => {
    const previous = {}; root.querySelectorAll("input").forEach(input => previous[input.dataset.name] = input.value); root.replaceChildren();
    if (!names.length) { const fixed = document.createElement("span"); fixed.className = "fixed-constraint"; fixed.textContent = "无自由参数（由 Data 完全固定）"; root.appendChild(fixed); return; }
    names.forEach(name => { const label = document.createElement("label"); label.textContent = name; const input = document.createElement("input"); input.className = className; input.dataset.name = name; input.value = previous[name] ?? defaults[name] ?? "0"; input.placeholder = `exact ${name}`; label.appendChild(input); root.appendChild(label); });
  };
  const updateOrbitHint = card => {
    const constrained = card.querySelector(".data-wyckoff-fields"), explicit = card.querySelector(".explicit-site-fields"), hint = card.querySelector(".selection-preview");
    if (state.symmetrySource !== "msg") { constrained.hidden = true; explicit.hidden = false; hint.textContent = "显式对称模式：这里是用户给定的代表 site seed，不是 Data Wyckoff 选择。"; return; }
    constrained.hidden = false; explicit.hidden = true; const item = selectedWyckoff(card); hint.replaceChildren();
    if (!item) { hint.textContent = "请选择 Wyckoff 位置"; return; }
    const title = document.createElement("strong"); title.textContent = `Wyckoff ${item.letter} · multiplicity ${item.multiplicity}`; hint.appendChild(title);
    const formulas = document.createElement("div"); formulas.className = "constraint-formulas";
    item.positions.forEach((position, index) => { const line = document.createElement("code"); line.textContent = `r${index + 1} = (${position.coordinates.join(", ")})   m${index + 1} = (${position.moment.join(", ")})`; formulas.appendChild(line); });
    hint.appendChild(formulas);
    constrained.querySelectorAll(".parameter-note").forEach(note => note.remove());
    renderParameterBindings(card.querySelector(".coordinate-parameter-fields"), item.coordinate_symbols, "orbit-coordinate-parameter", { ...(item.coordinate_parameter_defaults || {}), ...(card._coordinateDefaults || {}) });
    renderParameterBindings(card.querySelector(".moment-parameter-fields"), item.moment_symbols, "orbit-moment-parameter", { ...(item.moment_parameter_defaults || {}), ...(card._momentDefaults || {}) });
    if (item.coordinate_symbols.length) { const note = document.createElement("small"); note.className = "parameter-note"; note.textContent = "自由坐标先填入可见的非特殊演示值；若改到特殊位置使等价点重合，请改选对应的低 multiplicity Wyckoff 位。"; constrained.appendChild(note); }
    card._coordinateDefaults = {}; card._momentDefaults = {};
  };
  const addOrbit = (defaults = {}) => {
    const card = document.createElement("div"); card.className = "orbit-card";
    card._coordinateDefaults = defaults.coordinate_parameters || {}; card._momentDefaults = defaults.moment_parameters || {};
    const head = document.createElement("div"); head.className = "repeat-head"; const title = document.createElement("strong"); const remove = document.createElement("button"); remove.type = "button"; remove.className = "remove-button"; remove.textContent = "删除"; remove.onclick = () => { card.remove(); renumberOrbits(); }; head.append(title, remove);
    const preferred = defaults.source_ordinal && defaults.letter ? `${defaults.source_ordinal}:${defaults.letter}` : defaults.source_ordinal;
    const selectLabel = document.createElement("label"); selectLabel.textContent = "Wyckoff"; const select = document.createElement("select"); select.className = "orbit-wyckoff"; fillWyckoffSelect(select, preferred); select.onchange = () => updateOrbitHint(card); selectLabel.appendChild(select);
    const explicit = document.createElement("div"); explicit.className = "field-grid two explicit-site-fields"; const positionLabel = document.createElement("label"); positionLabel.textContent = "代表 position seed"; positionLabel.appendChild(makeVector("orbit-position", defaults.position || ["0", "0", "0"])); const momentLabel = document.createElement("label"); momentLabel.textContent = "代表 moment seed"; momentLabel.appendChild(makeVector("orbit-moment", defaults.moment || ["0", "0", "0"])); explicit.append(positionLabel, momentLabel);
    const constrained = document.createElement("div"); constrained.className = "data-wyckoff-fields"; const coordinateTitle = document.createElement("strong"); coordinateTitle.textContent = "位置自由参数"; const coordinateFields = document.createElement("div"); coordinateFields.className = "field-grid parameter-bindings coordinate-parameter-fields"; const momentTitle = document.createElement("strong"); momentTitle.textContent = "磁矩自由参数"; const momentFields = document.createElement("div"); momentFields.className = "field-grid parameter-bindings moment-parameter-fields"; constrained.append(coordinateTitle, coordinateFields, momentTitle, momentFields);
    const hint = document.createElement("div"); hint.className = "selection-preview";
    card.append(head, selectLabel, hint, constrained, explicit, basisPicker(defaults.basis || ["s"])); $("#orbitList").appendChild(card); renumberOrbits(); updateOrbitHint(card);
  };
  const renumberOrbits = () => $$("#orbitList .orbit-card").forEach((card, index) => card.querySelector(".repeat-head strong").textContent = `Orbit ${index + 1}`);
  const refreshOrbitSelectors = () => $$("#orbitList .orbit-card").forEach(card => { const select = card.querySelector(".orbit-wyckoff"); fillWyckoffSelect(select, select.value); updateOrbitHint(card); });
  const addOperation = (defaults = {}) => {
    const card = document.createElement("div"); card.className = "operation-card";
    const head = document.createElement("div"); head.className = "repeat-head"; const title = document.createElement("strong"); const remove = document.createElement("button"); remove.type = "button"; remove.className = "remove-button"; remove.textContent = "删除"; remove.onclick = () => { card.remove(); renumberOperations(); }; head.append(title, remove);
    const kind = document.createElement("select"); kind.className = "operation-kind"; kind.add(new Option("空间操作", "spatial")); kind.add(new Option("Spin-space 操作", "spin_space")); kind.value = defaults.kind || "spatial";
    const label = document.createElement("input"); label.className = "operation-label"; label.value = defaults.label || "E";
    const anti = document.createElement("input"); anti.type = "checkbox"; anti.className = "operation-antiunitary"; anti.checked = Boolean(defaults.antiunitary);
    const fields = document.createElement("div"); fields.className = "field-grid two"; const kindLabel = document.createElement("label"); kindLabel.textContent = "类型"; kindLabel.appendChild(kind); const nameLabel = document.createElement("label"); nameLabel.textContent = "Label"; nameLabel.appendChild(label); const antiLabel = document.createElement("label"); antiLabel.className = "check-row"; antiLabel.append(anti, document.createTextNode("反幺正 T")); fields.append(kindLabel, nameLabel, antiLabel);
    const rotationLabel = document.createElement("label"); rotationLabel.textContent = "空间 rotation"; const rotation = document.createElement("div"); rotation.className = "matrix-editor operation-rotation"; makeMatrixEditor(rotation, defaults.rotation || [["1","0","0"],["0","1","0"],["0","0","1"]], "operation-rotation-value"); rotationLabel.appendChild(rotation);
    const translationLabel = document.createElement("label"); translationLabel.textContent = "translation"; translationLabel.appendChild(makeVector("operation-translation", defaults.translation || ["0","0","0"]));
    const spin = document.createElement("div"); spin.className = "spin-operation-fields"; const spinLabel = document.createElement("label"); spinLabel.textContent = "spin rotation"; const spinMatrix = document.createElement("div"); spinMatrix.className = "matrix-editor operation-spin"; makeMatrixEditor(spinMatrix, defaults.spin_rotation || [["1","0","0"],["0","1","0"],["0","0","1"]], "operation-spin-value"); spinLabel.appendChild(spinMatrix); const parameterLabel = document.createElement("label"); parameterLabel.textContent = "continuous parameter（可空）"; const parameter = document.createElement("input"); parameter.className = "operation-continuous"; parameter.placeholder = "theta"; parameter.value = defaults.continuous_parameter || ""; parameterLabel.appendChild(parameter); spin.append(spinLabel, parameterLabel);
    const updateKind = () => spin.hidden = kind.value !== "spin_space"; kind.onchange = updateKind;
    card.append(head, fields, rotationLabel, translationLabel, spin); $("#operationList").appendChild(card); updateKind(); renumberOperations();
  };
  const renumberOperations = () => $$("#operationList .operation-card").forEach((card, index) => card.querySelector(".repeat-head strong").textContent = `Operation ${index + 1}`);
  const setSymmetrySource = source => {
    state.symmetrySource = source; $$('[data-symmetry-source]').forEach(button => button.classList.toggle("active", button.dataset.symmetrySource === source)); $("#msgSymmetryPanel").hidden = source !== "msg"; $("#explicitSymmetryPanel").hidden = source !== "explicit"; $("#msgLatticePanel").hidden = source !== "msg"; $("#explicitLatticePanel").hidden = source !== "explicit"; $("#addLatticeParameter").hidden = source !== "explicit"; $("#generateSymmetryGroup").disabled = source === "msg";
    if (source === "explicit") { if (!$("#operationList").children.length) addOperation(); renderExplicitLatticeParameters(); }
    else if ($("#msgSelector").value) loadBravais().catch(error => toast(error.message, true));
    refreshOrbitSelectors();
  };
  const loadWyckoff = async () => { const id = $("#msgSelector").value; if (!id) { state.wyckoffOptions = []; refreshOrbitSelectors(); return; } const payload = await request(`api/init/msg/${encodeURIComponent(id)}/wyckoff-options`); state.wyckoffOptions = payload.wyckoff || []; refreshOrbitSelectors(); };
  const loadBravais = async () => { const option = $("#msgSelector").selectedOptions[0]; if (!option?.dataset.record) return; const group = JSON.parse(option.dataset.record); const value = await request(`api/init/bravais/${encodeURIComponent(group.bravais_lattice_id)}`); state.bravais = value; renderBravaisPreview(value); renderBravaisParameters(value); };
  const loadMsgFamily = async () => {
    const sg = Number($("#spaceGroupNumber").value), payload = await request(`api/init/msg-family/${sg}`), select = $("#msgSelector"); select.replaceChildren();
    payload.groups.forEach(group => { const option = new Option(`${group.bns.join(".")} · ${group.symbol} · ${group.classification}`, group.stable_id); option.dataset.record = JSON.stringify(group); select.add(option); });
    updateMsgPreview(); await Promise.all([loadWyckoff(), loadBravais()]);
  };
  const updateMsgPreview = () => { const option = $("#msgSelector").selectedOptions[0]; if (!option?.dataset.record) { $("#msgPreview").textContent = "没有可用 MSG"; return; } const record = JSON.parse(option.dataset.record); $("#msgPreview").innerHTML = `<strong>BNS ${record.bns.join(".")} · ${escapeHtml(record.symbol)}</strong><br>${record.classification} · ${record.bravais_lattice_id} · ${record.operation_count} operations · stable ${record.stable_id}`; };
  const collectInitPayload = () => {
    const parameters = {}; if (state.symmetrySource === "msg") $$("#latticeParameters .bravais-parameter-value").forEach(input => parameters[input.dataset.name] = input.value.trim()); else $$("#latticeParameters .parameter-row").forEach(row => { const name = row.querySelector(".parameter-name").value.trim(), value = row.querySelector(".parameter-value").value.trim(); if (name && value) parameters[name] = value; });
    const lattice = { preset: state.symmetrySource === "msg" ? "msg_bravais" : $("#latticePreset").value, parameters }; if (state.symmetrySource === "explicit" && lattice.preset === "custom") lattice.matrix = readMatrix($("#customLattice"));
    const symmetry = { source: state.symmetrySource, msg_id: null, operations: [] };
    if (state.symmetrySource === "msg") symmetry.msg_id = $("#msgSelector").value;
    else symmetry.operations = $$("#operationList .operation-card").map(card => { const kind = card.querySelector(".operation-kind").value, operation = { kind, label: card.querySelector(".operation-label").value.trim(), rotation: readMatrix(card.querySelector(".operation-rotation"), ".operation-rotation-value"), translation: readVector(card, "operation-translation"), antiunitary: card.querySelector(".operation-antiunitary").checked }; if (kind === "spin_space") { operation.spin_rotation = readMatrix(card.querySelector(".operation-spin"), ".operation-spin-value"); operation.continuous_parameter = card.querySelector(".operation-continuous").value.trim() || null; } return operation; });
    const orbits = $$("#orbitList .orbit-card").map(card => { const selected = selectedWyckoff(card); const basis = [...card.querySelector(".basis-picker")._basisOrder]; if (!basis.length) throw new Error("每个 orbit 至少选择一个 basis function"); if (state.symmetrySource === "msg") { if (!selected) throw new Error("请选择有效的 Wyckoff 条目"); const coordinateParameters = {}, momentParameters = {}; card.querySelectorAll(".orbit-coordinate-parameter").forEach(input => coordinateParameters[input.dataset.name] = input.value.trim()); card.querySelectorAll(".orbit-moment-parameter").forEach(input => momentParameters[input.dataset.name] = input.value.trim()); return { letter: selected.letter, source_ordinal: selected.source_ordinal, coordinate_parameters: coordinateParameters, moment_parameters: momentParameters, basis_functions: basis }; } return { letter: "site", source_ordinal: null, position: readVector(card, "orbit-position"), moment: readVector(card, "orbit-moment"), basis_functions: basis }; });
    const field = { preset: $("#exactField").value }; if (field.preset === "custom") { field.conductor = Number($("#fieldConductor").value); field.cyclotomic_polynomial = $("#fieldPolynomial").value.split(",").map(value => Number(value.trim())); }
    return { label: $("#initLabel").value.trim(), lattice, symmetry, orbits, initial_bond_shells: Number($("#initialBondShells").value), generate_symmetry_group: $("#generateSymmetryGroup").checked, representation_mode: $("#representationMode").value, exact_field: field, direct_representation: state.loadedDirectRepresentation };
  };
  const applyInitExample = async example => {
    const payload = example.request;
    $("#initLabel").value = payload.label;
    $("#initialBondShells").value = String(payload.initial_bond_shells);
    $("#representationMode").value = payload.representation_mode;
    $("#generateSymmetryGroup").checked = Boolean(payload.generate_symmetry_group);
    $("#exactField").value = payload.exact_field.preset;
    $("#customField").hidden = payload.exact_field.preset !== "custom";
    state.loadedDirectRepresentation = payload.direct_representation || null;
    $("#directRepresentationNotice").hidden = !state.loadedDirectRepresentation;
    if (payload.symmetry.source === "msg") {
      setSymmetrySource("msg");
      $("#spaceGroupNumber").value = String(example.space_group_number);
      await loadMsgFamily();
      $("#msgSelector").value = payload.symmetry.msg_id;
      if (!$("#msgSelector").value) throw new Error(`示例 MSG ${payload.symmetry.msg_id} 不在当前稳定 Data 中`);
      updateMsgPreview();
      await Promise.all([loadWyckoff(), loadBravais()]);
      $$("#latticeParameters .bravais-parameter-value").forEach(input => {
        if (Object.hasOwn(payload.lattice.parameters, input.dataset.name)) input.value = payload.lattice.parameters[input.dataset.name];
      });
    } else {
      setSymmetrySource("explicit");
      $("#latticePreset").value = payload.lattice.preset;
      $("#customLattice").hidden = payload.lattice.preset !== "custom";
      if (payload.lattice.preset === "custom") makeMatrixEditor($("#customLattice"), payload.lattice.matrix);
      $("#latticeParameters").replaceChildren();
      Object.entries(payload.lattice.parameters || {}).forEach(([name, value]) => addLatticeParameter(name, value));
      $("#operationList").replaceChildren();
      payload.symmetry.operations.forEach(operation => addOperation(operation));
    }
    $("#orbitList").replaceChildren();
    payload.orbits.forEach(orbit => addOrbit(payload.symmetry.source === "msg" ? {
      letter: orbit.letter, source_ordinal: orbit.source_ordinal,
      coordinate_parameters: orbit.coordinate_parameters,
      moment_parameters: orbit.moment_parameters, basis: orbit.basis_functions,
    } : { position: orbit.position, moment: orbit.moment, basis: orbit.basis_functions }));
    toast(`已载入：${example.title}`);
  };
  const runInitExample = async example => {
    await applyInitExample(example);
    const result = await request("api/models/init", { method: "POST", body: JSON.stringify(collectInitPayload()) });
    state.modelParameterDefaults[result.model_id] = { ...(example.band_parameters || {}) };
    state.modelBandShells[result.model_id] = [...(example.band_shells || [1])];
    state.modelParameterSources[result.model_id] = example.band_parameter_source || "";
    await refreshModels(result.model_id);
    $("#propertyShells").value = state.modelBandShells[result.model_id].join(",");
    $("#shellInput").value = String(Math.max(...state.modelBandShells[result.model_id]));
    await refreshPropertyParameters(false);
    updateParameterHint();
    toast(`示例已运行：${result.label}`);
  };
  const renderInitExamples = () => {
    const root = $("#exampleList"); root.replaceChildren();
    state.initExamples.forEach(example => {
      const card = document.createElement("article"); card.className = "example-card";
      const title = document.createElement("strong"); title.textContent = example.title;
      const actions = document.createElement("div"); actions.className = "example-actions";
      const load = document.createElement("button"); load.type = "button"; load.className = "mini-button"; load.textContent = "载入"; load.setAttribute("aria-label", `载入 ${example.title}`); load.onclick = () => withLoading(card, () => applyInitExample(example)).catch(error => toast(error.message, true));
      const run = document.createElement("button"); run.type = "button"; run.className = "secondary-button"; run.textContent = "运行"; run.setAttribute("aria-label", `运行 ${example.title}`); run.onclick = () => withLoading(card, () => runInitExample(example)).catch(error => toast(error.message, true));
      actions.append(load, run); card.append(title, actions); root.appendChild(card);
    });
  };
  const initializeInitBuilder = async () => {
    const [options, examples] = await Promise.all([request("api/init/options"), request("api/init/examples")]); state.initOptions = options; state.initExamples = examples.examples || []; const sg = $("#spaceGroupNumber"); sg.replaceChildren(); state.initOptions.space_group_numbers.forEach(number => sg.add(new Option(String(number), String(number)))); renderInitExamples();
    makeMatrixEditor($("#customLattice"), [["1","0","0"],["0","1","0"],["0","0","1"]]); addOrbit(); await loadMsgFamily();
  };

  const boot = async () => {
    typeset(document);
    try {
      const health = await request("api/health");
      $("#healthBadge").innerHTML = `<i></i>Rust ${health.magnetictb_version}`;
      try { const data = await request("api/data/summary"); const counts = data.counts || {}; $("#dataCount").textContent = counts.msg_groups || "ready"; } catch (_) { $("#dataCount").textContent = "ready"; }
      await refreshModels(); await initializeInitBuilder();
    } catch (error) { $("#healthBadge").classList.add("error"); $("#healthBadge").innerHTML = "<i></i>连接失败"; toast(error.message, true); }
  };

  $("#quickStartButton").onclick = () => withLoading($("#quickStartButton"), async () => createModel("init", { initial_bond_shells: 2 }, "Default exact model")).catch(error => toast(error.message, true));
  $("#initForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => { const result = await request("api/models/init", { method: "POST", body: JSON.stringify(collectInitPayload()) }); await refreshModels(result.model_id); toast(`init 完成：${result.label}`); }).catch(error => toast(error.message, true)); };
  $("#modelForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => createModel($("#entrypoint").value, parse($("#modelArguments").value, "模型参数"), $("#modelLabel").value.trim())).catch(error => toast(error.message, true)); };
  $("#modelPicker").onchange = event => { state.modelId = event.target.value; refreshModels(state.modelId).catch(error => toast(error.message, true)); };
  $("#latticePreset").onchange = event => $("#customLattice").hidden = event.target.value !== "custom";
  $("#addLatticeParameter").onclick = () => addLatticeParameter();
  $$('[data-symmetry-source]').forEach(button => button.onclick = () => { state.loadedDirectRepresentation = null; $("#directRepresentationNotice").hidden = true; setSymmetrySource(button.dataset.symmetrySource); });
  $("#spaceGroupNumber").onchange = () => withLoading($("#msgSymmetryPanel"), loadMsgFamily).catch(error => toast(error.message, true));
  $("#msgSelector").onchange = () => { updateMsgPreview(); withLoading($("#msgSymmetryPanel"), async () => Promise.all([loadWyckoff(), loadBravais()])).catch(error => toast(error.message, true)); };
  $("#addOrbit").onclick = () => addOrbit();
  $("#addOperation").onclick = () => addOperation({ label: `g${$("#operationList").children.length + 1}` });
  $("#exactField").onchange = event => $("#customField").hidden = event.target.value !== "custom";

  $("#hamiltonianForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => {
    const modelId = encodeURIComponent(requireModel()), options = hamiltonianOptions();
    const shells = Array.from({ length: options.shell }, (_, index) => index + 1);
    const advanced = shells.length === 1
      ? await request(`api/models/${modelId}/hamiltonian`, { method: "POST", body: JSON.stringify(options) })
      : await request(`api/models/${modelId}/combine`, { method: "POST", body: JSON.stringify({ shells, hermitian: options.hermitian, validation_level: options.validation_level, kernel_method: options.kernel_method }) });
    state.hamiltonian = advanced;
    state.resultTab = "matrix"; renderHamiltonian();
    $("#propertyShells").value = shells.join(",");
    renderNumericParameters($("#parameterInputs"), state.hamiltonian.parameter_names || [], activeParameterDefaults());
    renderNumericParameters($("#propertyParameters"), state.hamiltonian.parameter_names || [], activeParameterDefaults());
    updateParameterHint();
  }).catch(error => toast(error.message, true)); };
  $$("[data-result-tab]").forEach(button => button.onclick = () => { state.resultTab = button.dataset.resultTab; renderHamiltonian(); });
  $("#evaluateForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => { const payload = { ...hamiltonianOptions(), cumulative: true, parameters: readNumericParameters($("#parameterInputs")), momentum: vectorValues($("#momentumInputs"), "k") }; const result = await request(`api/models/${encodeURIComponent(requireModel())}/evaluate`, { method: "POST", body: JSON.stringify(payload) }); $("#evaluateResult").innerHTML = `<div class="ham-matrix" data-tex="${escapeHtml(matrixFromNumeric(result.matrix))}"></div>`; typeset($("#evaluateResult")); }).catch(error => toast(error.message, true)); };
  $("#spaceButton").onclick = () => withLoading($("#spaceButton"), async () => { const result = await request(`api/models/${encodeURIComponent(requireModel())}/hamiltonian-space`, { method: "POST", body: JSON.stringify(hamiltonianOptions()) }); $("#spaceResult").innerHTML = `<pre class="json-result">${escapeHtml(json(result))}</pre>`; }).catch(error => toast(error.message, true));
  $("#reportForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => {
    const modelId = encodeURIComponent(requireModel()), kind = $("#reportKind").value; let endpoint, body = null;
    if (kind === "orbitals") endpoint = "orbital-table";
    else if (kind === "basis") { endpoint = "hamiltonian-basis"; body = { row: Number($("#reportRow").value), column: Number($("#reportColumn").value) }; }
    else if (kind === "bonds") { endpoint = "bonds"; body = { shell: Number($("#reportShell").value) }; }
    else if (kind === "symmetry") { endpoint = "symmetry-representations"; body = { selection: $("#reportSelection").value }; }
    else if (kind === "unsymham") { const result = await request(`api/models/${modelId}/unsymham`, { method: "POST", body: JSON.stringify({ shell: Number($("#reportShell").value) }) }); const root = $("#reportResult"); root.className = "report-result"; root.innerHTML = `<div class="result-meta"><span>unsymham shell ${result.shell}</span><span>${result.parameter_names.map(parameterTex).join(", ")}</span></div><div class="ham-matrix" data-tex="${escapeHtml(matrixTex(result.matrix_text))}"></div>`; typeset(root); return; }
    else { endpoint = "hopping-parameters"; body = { ...hamiltonianOptions(), shell: Number($("#reportShell").value), parameter: $("#reportParameter").value.trim() || null }; }
    const result = await request(`api/models/${modelId}/reports/${endpoint}`, body === null ? {} : { method: "POST", body: JSON.stringify(body) }); renderTableReport($("#reportResult"), result);
  }).catch(error => toast(error.message, true)); };

  const propertyRequest = async (endpoint, payload, root) => { const result = await request(`api/models/${encodeURIComponent(requireModel())}/properties/${endpoint}`, { method: "POST", body: JSON.stringify(payload) }); renderPropertyResult(root, result); return result; };
  const plotRequest = async (endpoint, payload) => request(`api/models/${encodeURIComponent(requireModel())}/plots/${endpoint}`, { method: "POST", body: JSON.stringify(payload) });
  const refreshPropertyParameters = async (preserve = true) => { const body = { ...propertyBase(), parameters: {} }; const result = await request(`api/models/${encodeURIComponent(requireModel())}/properties/parameter-names`, { method: "POST", body: JSON.stringify(body) }); renderNumericParameters($("#propertyParameters"), result.parameter_names || [], activeParameterDefaults(), preserve); updateParameterHint(); };
  $("#propertyShells").onchange = () => withLoading($("#propertyShells"), refreshPropertyParameters).catch(error => toast(error.message, true));
  $("#crystalForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => { const result = await plotRequest("crystal-structure", { cell_range: Number($("#crystalCellRange").value), show_atom_labels: $("#crystalLabels").checked }); renderCrystalStructure($("#crystalResult"), result); }).catch(error => toast(error.message, true)); };
  $("#brillouinForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => { const result = await plotRequest("brillouin-zone", { k_path: "Automatic", show_k_path: $("#brillouinPath").checked, translation_range: Number($("#brillouinRange").value), tolerance: finiteNumber($("#brillouinTolerance").value, "Tolerance") }); renderBrillouinZone($("#brillouinResult"), result); }).catch(error => toast(error.message, true)); };
  $("#bandForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => { const result = await request(`api/models/${encodeURIComponent(requireModel())}/properties/bands`, { method: "POST", body: JSON.stringify({ ...propertyBase(), npoint: Number($("#bandNpoint").value), bravais_type: "Automatic", tolerance: 1e-6 }) }); renderBandPlot($("#bandResult"), result); }).catch(error => toast(error.message, true)); };
  $("#hoppingDataButton").onclick = event => withLoading(event.currentTarget, async () => { const result = await propertyRequest("hopping-data", propertyBase(), $("#blochResult")); $("#blochResult").innerHTML = `<div class="result-meta"><span>${result.NumWannier} orbitals</span><span>${result.NumTranslations} translations</span><span>residual ${Number(result.HermitianResidual).toExponential(2)}</span></div><pre class="code-result">${escapeHtml(result.Translations.map((translation, index) => `${translation.join(" ")}  H(R) ${result.HoppingMatrices[index].length}×${result.HoppingMatrices[index][0].length}`).join("\n"))}</pre>`; }).catch(error => toast(error.message, true));
  $("#blochForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => propertyRequest("bloch", { ...propertyBase(), momentum: vectorValues($("#blochMomentum"), "q") }, $("#blochResult"))).catch(error => toast(error.message, true)); };
  $("#realSpaceForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => propertyRequest("real-space", { ...propertyBase(), geometry: vectorValues($("#finiteSize"), "size").map(Number), boundary_conditions: [...$("#finiteBoundary").querySelectorAll("select")].map(select => select.value) }, $("#realSpaceResult"))).catch(error => toast(error.message, true)); };
  $("#slabForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => propertyRequest("slab", { ...propertyBase(), size: [1, 1, Number($("#slabLayers").value)], momentum: vectorValues($("#slabMomentum"), "q"), periodic_directions: [1, 2] }, $("#slabResult"))).catch(error => toast(error.message, true)); };
  $("#surfaceForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => propertyRequest("surface", { ...propertyBase(), momentum: vectorValues($("#surfaceMomentum"), "q"), energy: finiteNumber($("#surfaceEnergy").value, "energy"), broadening: finiteNumber($("#surfaceBroadening").value, "broadening") }, $("#surfaceResult"))).catch(error => toast(error.message, true)); };
  const wilsonPayload = () => ({ ...propertyBase(), occupied: Number($("#wilsonOccupied").value), start: vectorValues($("#wilsonStart"), "start"), end: vectorValues($("#wilsonEnd"), "end"), path_subdivisions: Number($("#wilsonSubdivisions").value) });
  $("#wilsonForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => propertyRequest("wilson-loop", wilsonPayload(), $("#wilsonResult"))).catch(error => toast(error.message, true)); };
  $("#berryPhaseButton").onclick = event => withLoading(event.currentTarget, async () => { const payload = wilsonPayload(), path = Array.from({ length: payload.path_subdivisions + 1 }, (_, index) => payload.start.map((value, axis) => value + (payload.end[axis] - value) * index / payload.path_subdivisions)); await propertyRequest("berry-phase", { ...propertyBase(), occupied: payload.occupied, path }, $("#wilsonResult")); }).catch(error => toast(error.message, true));
  $("#curvatureForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => propertyRequest("berry-curvature", { ...propertyBase(), occupied: Number($("#curvatureOccupied").value), point: vectorValues($("#curvaturePoint"), "k"), directions: [1, 2], step_size: finiteNumber($("#curvatureStep").value, "step") }, $("#curvatureResult"))).catch(error => toast(error.message, true)); };
  $("#chernForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => propertyRequest("point-chern", { ...propertyBase(), occupied: Number($("#chernOccupied").value), point: vectorValues($("#chernPoint"), "k"), radius: finiteNumber($("#chernRadius").value, "radius"), surface_subdivisions: Number($("#chernSubdivisions").value) }, $("#chernResult"))).catch(error => toast(error.message, true)); };
  $("#gaplessForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => propertyRequest("gapless-points", { ...propertyBase(), occupied: Number($("#gaplessOccupied").value), grid_size: Number($("#gaplessGrid").value), candidate_count: Number($("#gaplessCandidates").value), refinement_method: $("#gaplessMethod").value }, $("#gaplessResult"))).catch(error => toast(error.message, true)); };
  $("#surfaceSpectrumForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => {
    const result = await plotRequest("surface-spectrum", { ...propertyBase(), momentum_path: [[[[0, 0], [0.5, 0]], ["\\Gamma", "X"]]], energy_range: [finiteNumber($("#surfaceSpectrumMin").value, "energy min"), finiteNumber($("#surfaceSpectrumMax").value, "energy max")], momentum_subdivisions: Number($("#surfaceSpectrumSubdivisions").value), energy_points: Number($("#surfaceSpectrumPoints").value), broadening: finiteNumber($("#surfaceBroadening").value, "broadening") });
    renderHeatmap($("#surfaceSpectrumResult"), result.path_coordinate, result.energies, result.spectral_weight, "surface spectral weight");
  }).catch(error => toast(error.message, true)); };
  $("#wilsonPlotForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => {
    const offset = commaVector($("#wilsonOffsetEnd").value, "参数路径终点"); if (offset.length !== 3) throw new Error("参数路径终点必须包含三个坐标");
    const result = await plotRequest("wilson-loop", { ...wilsonPayload(), parameter_path: [[[[0, 0, 0], offset], ["0", "λ"]]], parameter_subdivisions: Number($("#wilsonParameterSubdivisions").value) });
    renderBranchPlot($("#wilsonPlotResult"), result);
  }).catch(error => toast(error.message, true)); };
  $("#berryGridForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => {
    const minimum = finiteNumber($("#berryGridMin").value, "k min"), maximum = finiteNumber($("#berryGridMax").value, "k max"), size = Number($("#berryGridSize").value);
    const result = await plotRequest("berry-curvature-2d", { ...propertyBase(), occupied: Number($("#curvatureOccupied").value), ranges: [[minimum, maximum], [minimum, maximum]], directions: [1, 2], fixed_coordinates: "Automatic", grid_size: size, step_size: finiteNumber($("#berryGridStep").value, "step") });
    renderHeatmap($("#berryGridResult"), result.grid_values[0], result.grid_values[1], result.curvature, "Berry curvature F12", true);
  }).catch(error => toast(error.message, true)); };
  $("#berryVectorButton").onclick = event => withLoading(event.currentTarget, async () => {
    const minimum = finiteNumber($("#berryGridMin").value, "k min"), maximum = finiteNumber($("#berryGridMax").value, "k max"), size = Math.min(Number($("#berryGridSize").value), 7);
    const result = await plotRequest("berry-curvature-3d", { ...propertyBase(), occupied: Number($("#curvatureOccupied").value), ranges: [[minimum, maximum], [minimum, maximum], [minimum, maximum]], grid_size: size, step_size: finiteNumber($("#berryGridStep").value, "step") });
    $("#berryGridResult").className = "plot-result"; $("#berryGridResult").innerHTML = `<div class="result-meta"><span>${result.visible_vector_count} vectors</span><span>grid ${result.grid_size.join(" × ")}</span><span>max |F| ${Number(result.maximum_magnitude).toPrecision(7)}</span><span>components F23, F31, F12</span></div><pre class="code-result">${escapeHtml(json({ component_directions: result.component_directions, step_size: result.step_size, vector_scale: result.vector_scale }))}</pre>`;
  }).catch(error => toast(error.message, true));
  $("#wavefunctionForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => {
    const result = await plotRequest("real-space-wavefunction", { ...propertyBase(), geometry: vectorValues($("#wavefunctionSize"), "cell size").map(Number), state: amplitudeValues($("#wavefunctionState").value), aggregation: $("#wavefunctionAggregation").value, normalize: true, phase_coloring: $("#wavefunctionPhase").checked });
    renderWavefunction($("#wavefunctionResult"), result);
  }).catch(error => toast(error.message, true)); };

  $("#msgopForm").onsubmit = event => { event.preventDefault(); withLoading(event.currentTarget, async () => { const kind = $("#msgopKind").value; let result; if (["rod_og", "layer_og", "grayrod", "graylayer"].includes(kind)) { const values = $("#msgopValue").value.split(",").map(item => Number(item.trim())); const gray = kind.startsWith("gray"), family = kind.includes("rod") ? "rod" : "layer"; if (values.some(value => !Number.isInteger(value) || value < 1) || values.length !== (gray ? 1 : 3)) throw new Error(gray ? `${kind} 需要一个正整数` : `${kind} 需要 n,m,k`); result = await request("api/data/subperiodic-operations", { method: "POST", body: JSON.stringify({ kind: family, selector: gray ? "gray" : "og", key: values }) }); } else result = await request("api/msgop", { method: "POST", body: JSON.stringify({ group: groupSelector() }) }); renderOperations($("#msgopResult"), result); }).catch(error => toast(error.message, true)); };
  $("#wyckoffReportButton").onclick = event => withLoading(event.currentTarget, async () => { if (["rod_og", "layer_og", "grayrod", "graylayer"].includes($("#msgopKind").value)) throw new Error("Wyckoff 表只适用于 MSG 选择器"); const result = await request("api/data/msg-wyckoff-report", { method: "POST", body: JSON.stringify({ group: groupSelector() }) }); renderTableReport($("#msgopResult"), result); }).catch(error => toast(error.message, true));

  $("#menuButton").onclick = () => $(".sidebar").classList.toggle("open");
  $$(".nav a").forEach(link => link.onclick = () => { $$(".nav a").forEach(item => item.classList.remove("active")); link.classList.add("active"); $(".sidebar").classList.remove("open"); });
  $("#languageButton").onclick = () => { state.locale = i18n.toggle(); };
  boot();
})();
