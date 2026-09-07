<|
  "Options" -> "选项",
  "Basic Examples" -> "基本范例",
  "Define the two-band Qi-Wu-Zhang Hamiltonian and place both basis orbitals at the unit-cell origin:" ->
    "定义两带 Qi-Wu-Zhang Hamiltonian，并把两个基轨道都放在原胞原点：",
  "Initialize a one-orbital simple-cubic model and define a two-segment reciprocal-space path:" ->
    "初始化单轨道简单立方模型，并定义一条包含两段的倒空间路径：",
  "Draw the one-band dispersion along the explicit path:" ->
    "沿显式路径绘制单带色散：",
  "Initialize the one-orbital simple-cubic model used by the automatic path:" ->
    "初始化自动路径所使用的单轨道简单立方模型：",
  "Let showband choose the conventional simple-cubic path and draw the dispersion:" ->
    "让 showband 选择约定的简单立方路径并绘制色散：",
  "Initialize a one-orbital simple-cubic model before requesting its conventional path:" ->
    "先初始化单轨道简单立方模型，再读取其约定路径：",
  "Initialize a one-orbital simple-cubic model before drawing its reciprocal cell:" ->
    "先初始化单轨道简单立方模型，再绘制其倒空间原胞：",
  StringJoin[
    "Inspect the plaquette flux, area, and local curvature together at Gamma. ",
    "The value near 0.5 is a curvature density at one momentum point, not a Chern number:"
  ] -> StringJoin[
    "在 Gamma 点同时查看小回路的 Berry 通量、面积和局域曲率。",
    "接近 0.5 的数值是单个动量点上的曲率密度，不是 Chern 数："
  ],
  StringJoin[
    "Plot the same occupied-band curvature over the full Brillouin zone. ",
    "Only its Brillouin-zone integral divided by 2 Pi is the Chern number; ",
    "for this oriented Qi-Wu-Zhang model that integral is 1:"
  ] -> StringJoin[
    "在整个 Brillouin 区绘制同一占据带的 Berry 曲率。",
    "只有全 Brillouin 区积分再除以 2 Pi 才是 Chern 数；",
    "按照这里的取向约定，Qi-Wu-Zhang 模型的积分结果为 1："
  ],
  "Applications" -> "应用",
  "Properties & Relations" -> "属性和关系",
  "Possible Issues" -> "可能存在的问题",
  "All ten public options are listed below with their actual defaults. lattice, lattpar, wyckoffposition, symminformation, and basisFunctions specify the physical model; the remaining options control compilation without changing those physical inputs." -> "下面列出全部十个公开选项及其实际默认值。lattice、lattpar、wyckoffposition、symminformation 和 basisFunctions 指定物理模型，其余选项控制编译过程而不改变这些物理输入。",
  "Primitive direct-lattice vectors stored by rows. After applying lattpar, the result must be a real, nonsingular numerical 3 by 3 matrix." -> "按行存储的原胞正格矢。代入 lattpar 后，结果必须是实的、非奇异的 3×3 数值矩阵。",
  "Replacement rules that assign numerical values to every symbolic lattice parameter needed by lattice." -> "为 lattice 中用到的每个符号晶格参数指定数值的替换规则。",
  "One {fractional seed position, magnetic moment} pair for each Wyckoff orbit. The outer-list order also fixes the order of basisFunctions." -> "每个 Wyckoff 轨道对应一个 {分数坐标种子位置, 磁矩}。外层列表顺序同时固定 basisFunctions 的顺序。",
  "Either ordered magnetic-group records {label,R,t,\"F\"|\"T\"} or a discrete spin-space-group Association. The antiunitary flag belongs to each operation record." -> "可以是有序磁群记录 {label,R,t,\"F\"|\"T\"}，也可以是离散自旋空间群 Association。反幺正标志属于相应对称操作记录。",
  "A nonempty local scalar-orbital or two-component spinor basis for every Wyckoff seed. Named orbitals and exact symbolic functions are accepted." -> "为每个 Wyckoff 种子指定非空的局域标量轨道基或二分量旋量基；既可使用命名轨道，也可使用精确符号函数。",
  "Compatibility debug flag retained in the public init syntax. It never enables a hidden fallback or changes the mathematical model." -> "为保持 init 输入兼容而保留的调试标志；它不会启用隐藏 fallback，也不会改变数学模型。",
  "Positive integer giving the exact number of complete bond shells searched and statically compiled by init. Shell 1 is the zero-distance onsite shell." -> "正整数，指定 init 搜索并静态编译的完整键层数。第 1 层是零距离 onsite 键层。",
  "False treats the supplied ordered records as the complete finite group. True closes supplied generators using the operation-specific multiplication law." -> "False 把输入的有序记录直接视为完整有限群；True 使用相应对称操作的乘法规则闭合输入生成元。",
  "\"DirectProduct\" combines the site permutation with full local actions. \"Induced\" constructs the representation from a reference-site stabilizer and coset representatives." -> "\"DirectProduct\" 组合格点置换与完整局域作用；\"Induced\" 从参考点稳定子及陪集代表构造表示。",
  "In Induced mode, Automatic derives each reference-site stabilizer representation from basisFunctions. An explicit list may instead supply ReferenceSiteIndex, SiteSymmetryOperationIndices, SiteSymmetryMatrices, and optionally CosetRepresentativeIndices." -> "在 Induced 模式中，Automatic 从 basisFunctions 推导各参考点的稳定子表示；也可显式提供列表，其中包含 ReferenceSiteIndex、SiteSymmetryOperationIndices、SiteSymmetryMatrices，并可选给出 CosetRepresentativeIndices。",
  "InitialBondShells fixes the prepared range. This complete graphene input prepares only the onsite shell, whose actual Hamiltonian is displayed:" -> "InitialBondShells 固定初始化准备的范围。这个完整石墨烯输入只准备 onsite 键层，并直接显示其实际 Hamiltonian：",
  "GenerateSymmetryGroup -> True accepts generators instead of a complete operation list. These two spin-space-group generators close to the validated eight-element collinear group; the displayed matrices are the ones used by the Hamiltonian constraints:" -> "GenerateSymmetryGroup -> True 接受生成元而不是完整操作列表。下面两个自旋空间群生成元闭合为已经验证的八元共线群；所显示的矩阵正是 Hamiltonian 约束实际使用的矩阵：",
  "In cubic symmetry, RepresentationMode -> \"DirectProduct\" requires the complete local px, py, pz basis because rotations mix all three orbitals. The orbital table therefore contains all three orbitals on every symmetry-equivalent site:" -> "在立方对称性下，旋转会混合 px、py、pz，因此 RepresentationMode -> \"DirectProduct\" 必须输入完整的三个局域轨道。轨道表会在每个对称等价格点列出这三个轨道：",
  "For the same cubic Wyckoff orbit, RepresentationMode -> \"Induced\" needs only the reference px orbital at {1/2, 0, 0}. Cubic symmetry transports it to the corresponding py and pz states at the other equivalent sites, as the orbital table shows:" -> "对于同一个立方 Wyckoff 轨道，RepresentationMode -> \"Induced\" 只需输入 {1/2, 0, 0} 参考点上的 px 轨道。立方对称性会把它搬运为其他等价格点上相应的 py、pz 态，轨道表直接显示了这个结果：",
  "SiteLocalData can replace Automatic only in Induced mode. In this E plus PT model, E stabilizes the reference site and the antiunitary PT operation transports it to the partner site; the antiunitary flag remains in symminformation, while SiteSymmetryMatrices contains only the local matrix:" -> "SiteLocalData 只在 Induced 模式中用于替代 Automatic。在这个 E+PT 模型中，E 稳定参考点，反幺正 PT 操作把它搬运到伙伴点；反幺正标志仍保存在 symminformation 中，而 SiteSymmetryMatrices 只包含局域矩阵：",
  "Wannier90 input and output" -> "Wannier90 输入与输出",
  "Symmetry databases" -> "对称性数据库",
  "Load MagneticTB" -> "加载 MagneticTB",
  "Related functions and tutorials" -> "相关函数与教程",
  "Band fitting and VASP comparison" -> "能带拟合与 VASP 对比",

  "MagneticTB constructs symmetry-constrained tight-binding Hamiltonians for magnetic, nonmagnetic, and discrete spin-space-group models." -> "MagneticTB 为磁性、非磁性及离散自旋空间群模型构造满足对称性约束的紧束缚 Hamiltonian。",
  " — obtain magnetic-space-group operations from the bundled database." -> " — 从随包数据库读取磁空间群操作。",
  " — initialize a model from basis functions." -> " — 从基函数初始化模型。",
  " — initialize a model from exact representation matrices." -> " — 从精确表示矩阵初始化模型。",
  " — inspect the current compiled model session." -> " — 查看当前已编译模型会话。",
  " — inspect the fixed Hamiltonian orbital order." -> " — 查看固定的 Hamiltonian 轨道顺序。",
  " — display full matrices for selected symmetry operations." -> " — 显示所选对称操作的完整表示矩阵。",
  " — display the Hamiltonian bra-ket order." -> " — 显示 Hamiltonian 的左、右基矢顺序。",
  " — prepare explicit subgroup initialization rules." -> " — 生成显式子群初始化规则。",
  " — solve a prepared bond shell and construct its Hamiltonian contribution." -> " — 求解已准备的键层并构造相应 Hamiltonian。",
  " — change the Bloch phase convention." -> " — 改变 Bloch 相位约定。",
  " — display prepared bond geometry." -> " — 显示已准备的键几何。",
  " — trace each parameter to real-space hoppings." -> " — 追踪每个参数对应的实空间跃迁。",
  " — explore the band structure interactively." -> " — 交互式查看能带结构。",
  " — read selected bands from a complete numeric VASP EIGENVAL file." -> " — 从完整数值 VASP EIGENVAL 文件读取指定能带。",
  " — tune a symbolic tight-binding Hamiltonian interactively against reference bands." -> " — 交互调节符号紧束缚 Hamiltonian，使其匹配参考能带。",
  " — fit an affine-linear tight-binding Hamiltonian to reference bands." -> " — 把仿射线性紧束缚 Hamiltonian 拟合到参考能带。",
  " — overlay fitted tight-binding bands with supplied reference bands." -> " — 叠加显示拟合后的紧束缚能带与输入的参考能带。",
  " — make a static numerical band plot." -> " — 绘制固定参数的数值能带图。",
  " — export numerical band eigenvalues." -> " — 导出数值能带本征值。",
  " — export solved shells or an exponential H(k) to wannier90_hr.dat." -> " — 将已求解键层或指数形式 H(k) 导出为 wannier90_hr.dat。",
  " — read and validate a Wannier90 HR file." -> " — 读取并验证 Wannier90 HR 文件。",
  " — display magnetic Wyckoff data." -> " — 显示磁性 Wyckoff 数据。",
  " — read magnetic layer-group operations." -> " — 读取磁层群操作。",
  " — read magnetic rod-group operations." -> " — 读取磁杆群操作。",
  " — initialize a model with the archived backend." -> " — 使用归档后端初始化模型。",
  " — construct a Hamiltonian with the archived backend." -> " — 使用归档后端构造 Hamiltonian。",
  " — inspect legacy bands interactively." -> " — 交互式查看旧版能带。",

  "exports the numerical band eigenvalues of H along path to file." -> "将 H 沿 path 的数值能带本征值导出到文件。",
  "The output file contains one row per band and one column per sampled momentum." -> "输出文件每条能带占一行，每个采样动量占一列。",
  "Use a portable path such as FileNameJoin[{$TemporaryDirectory,...}] in reusable notebooks." -> "可复用 notebook 应使用 FileNameJoin[{$TemporaryDirectory,...}] 这类可移植路径。",
  "Initialize graphene, export a short band table, and display its first five numerical rows:" -> "初始化石墨烯，导出一份简短能带表，并直接显示前五行数值数据：",

  "plots numerical eigenvalues of H along a reciprocal-fractional momentum path after applying parameter rules." -> "代入参数规则后，沿倒空间分数坐标路径绘制 H 的数值本征值。",
  "Path coordinates are multiplied by 2 Pi internally." -> "路径坐标在函数内部乘以 2 Pi。",
  "The plotRange option defaults to All and changes only the displayed vertical range." -> "plotRange 选项默认为 All，只改变图中显示的纵轴范围。",
  "Construct the complete graphene Hamiltonian and plot a fixed numerical parameter choice:" -> "构造完整的石墨烯 Hamiltonian，并绘制一组固定数值参数对应的能带：",

  "returns complete init rules for the subgroup generated by selected current operation indices or unambiguous labels." -> "返回由当前操作编号或无歧义标签生成的子群所对应的完整 init 规则。",
  "The function never changes the current model; inspect the returned rules and call init@@rules explicitly." -> "该函数不会修改当前模型；请先检查返回规则，再显式调用 init@@rules。",
  "Wyckoff orbits are split without dropping atoms when the retained subgroup no longer relates all sites." -> "当保留的子群不再联系全部格点时，Wyckoff 轨道会被拆分，但不会丢失原子。",
  "Use {} for the identity-only subgroup." -> "仅保留单位元子群时使用 {}。",
  "Initialize graphene and return the initialization rules for the subgroup generated by operation 2:" -> "初始化石墨烯，并直接返回由第 2 个操作生成的子群所对应的初始化规则：",

  "converts cached real-space hopping data for shells 1 through n into Wannier90 HR format." -> "将第 1 到 n 键层的缓存实空间跃迁转换为 Wannier90 HR 格式。",
  "extracts real-space coefficients from a hand-written exponential Bloch Hamiltonian." -> "从手写的指数形式 Bloch Hamiltonian 中提取实空间系数。",
  "The selected shells must already have been solved by symham with matching Hermitian, KernelMethod, and ValidationLevel options." -> "所选键层必须已由 symham 使用一致的 Hermitian、KernelMethod 和 ValidationLevel 选项求解。",
  "Set \"hrExport\" to a path to write wannier90_hr.dat; None returns the generated text." -> "将 \"hrExport\" 设为路径可写出 wannier90_hr.dat；设为 None 则返回生成的文本。",
  "Set \"wcc\" -> Automatic to use the orbital centers prepared by init." -> "设定 \"wcc\" -> Automatic 可使用 init 准备的轨道中心。",
  "Initialize graphene, solve three shells, write a real Wannier90 HR file, and inspect its header and first hopping records:" -> "初始化石墨烯、求解三个键层、写出真实的 Wannier90 HR 文件，并检查文件头和前几条跃迁记录：",

  "returns ordered magnetic layer-group operations from the bundled database." -> "从随包数据库返回有序磁层群操作。",
  "Use graylayer[n] to select a gray magnetic layer group." -> "使用 graylayer[n] 选择灰磁层群。",
  "Each record contains a label, point matrix, translation, and antiunitary flag." -> "每条记录包含标签、点操作矩阵、平移和反幺正标志。",
  "Inspect the first four operations of gray layer group 25:" -> "查看第 25 号灰磁层群的前四个操作：",
  "returns ordered magnetic rod-group operations from the bundled database." -> "从随包数据库返回有序磁杆群操作。",
  "Use grayrod[n] to select a gray magnetic rod group." -> "使用 grayrod[n] 选择灰磁杆群。",
  "Inspect the first four operations of gray rod group 25:" -> "查看第 25 号灰磁杆群的前四个操作：",

  "reads a complete Wannier90 HR file and returns structured translation vectors and hopping matrices." -> "读取完整的 Wannier90 HR 文件，并返回结构化的平移向量和跃迁矩阵。",
  "Malformed or incomplete files are rejected rather than partially interpreted." -> "格式错误或不完整的文件会被拒绝，不会被部分解释。",
  "The ncell option can impose an explicit cell cutoff after parsing." -> "ncell 选项可在解析后施加显式晶胞截断。",
  "Create a graphene HR file with hop, read it back, and display the actual translations and real-space matrices:" -> "用 hop 创建石墨烯 HR 文件，再读回并显示实际平移和实空间矩阵：",

  "displays the fixed Hamiltonian row and column basis prepared by init or initfromrep." -> "显示由 init 或 initfromrep 准备的固定 Hamiltonian 行、列基底。",
  "identifies the bra and ket states associated with Hamiltonian element H[[i,j]]." -> "标识 Hamiltonian 元素 H[[i,j]] 对应的左矢和右矢。",
  "The displayed order is identical to orbitalTable and to the full representation matrices." -> "显示顺序与 orbitalTable 及完整表示矩阵完全一致。",
  "Induced mode additionally displays the transported basis state and transport operation." -> "Induced 模式还会显示搬运后的基态和搬运操作。",
  "Initialize graphene and display the complete two-state Hamiltonian basis:" -> "初始化石墨烯并显示完整的两态 Hamiltonian 基底：",
  "Identify the bra and ket states of the nearest-neighbour matrix element H[[1,2]]:" -> "标识最近邻矩阵元 H[[1,2]] 的左、右基态：",

  "displays the representative real-space matrix element from which every independent parameter of cached symham[n] originates." -> "显示缓存 symham[n] 中每个独立参数所来源的代表性实空间矩阵元。",
  "displays all symmetry-propagated real-space occurrences of parameter p in shell n." -> "显示参数 p 在第 n 键层中由对称性传播得到的全部实空间项。",
  "Evaluate symham[n] first; this inspection function never solves a missing shell implicitly." -> "必须先计算 symham[n]；该检查函数不会隐式求解缺失键层。",
  "Rows identify source and destination orbitals, cell translations, coefficients, and representative bonds." -> "各行给出起点与终点轨道、晶胞平移、系数和代表键。",
  "Initialize graphene, solve its nearest-neighbour shell, and display the actual origin of t1:" -> "初始化石墨烯、求解最近邻键层，并显示 t1 的实际来源：",
  "Display every real-space occurrence generated from the representative t1 hopping:" -> "显示由代表性 t1 跃迁生成的每个实空间项：",

  "displays the Wyckoff positions and allowed magnetic-moment directions stored for a magnetic space group." -> "显示磁空间群数据库中的 Wyckoff 位置和允许磁矩方向。",
  "group can be a database identifier such as gray[n] or the corresponding magnetic-group number." -> "group 可以是 gray[n] 等数据库标识，也可以是对应的磁群编号。",
  "Symbolic x, y, and z entries denote free Wyckoff coordinates." -> "符号 x、y、z 表示自由的 Wyckoff 坐标。",
  "Display the actual Wyckoff table of gray magnetic space group 193:" -> "显示第 193 号灰磁空间群的实际 Wyckoff 表：",

  "displays prepared full representation matrices and their spatial, spin, and antiunitary operation data." -> "显示已准备的完整表示矩阵及其空间、自旋和反幺正操作信息。",
  "The default selection Automatic displays the compact generators; All displays every prepared operation." -> "默认选择 Automatic 显示紧凑生成元；All 显示全部已准备操作。",
  "The antiunitary flag belongs to the symmetry operation, so the display shows the matrix itself rather than appending a formal K symbol." -> "反幺正标志属于对称操作本身，因此这里只显示矩阵，不另加形式上的 K 符号。",
  "This function reads the prepared representation and performs no group generation or Hamiltonian solving." -> "该函数只读取已准备表示，不生成群，也不求解 Hamiltonian。",
  "Initialize graphene and display the matrices of its four compact generators:" -> "初始化石墨烯并显示四个紧凑生成元的矩阵：",
  "Select one prepared operation by index. Operation 25 is pure time reversal in this ordered gray group:" -> "按编号选择一个已准备操作。在该有序灰群中，第 25 个操作是纯时间反演：",

  "rewrites a Bloch Hamiltonian from convention I to the alternative phase convention II using the prepared orbital centers." -> "利用已准备的轨道中心，将 Bloch Hamiltonian 从约定 I 改写到相位约定 II。",
  "Initialize the model before calling symhamII so the orbital-center convention is defined." -> "调用 symhamII 前必须初始化模型，以确定轨道中心约定。",
  "The transformation changes phase placement but not the physical spectrum." -> "该变换只改变相位所在位置，不改变物理能谱。",
  "Construct graphene and display the same Hamiltonian in convention II:" -> "构造石墨烯并以约定 II 显示同一个 Hamiltonian：",

  "Initialize graphene and construct the Hamiltonian before plotting:" -> "绘图前先初始化石墨烯并构造 Hamiltonian：",
  "Define the G-M-K-G path in reciprocal fractional coordinates:" -> "用倒空间分数坐标定义 G-M-K-G 路径：",
  "Open the live interactive panel along the explicit G-M-K-G path. Its three sliders are the onsite, nearest-neighbour, and next-nearest-neighbour amplitudes; no hidden energy shift is added:" -> "沿显式 G-M-K-G 路径打开实时交互面板。三个滑块分别对应 onsite、最近邻和次近邻振幅；程序不加入隐藏能量平移：",
  "Display the conventional path selected automatically from the initialized graphene lattice:" -> "显示根据当前已初始化的石墨烯晶格自动选出的约定路径：",
  "Use the automatically selected path to open the live Manipulate panel:" -> "使用自动选出的路径打开实时 Manipulate 面板：",
  "A fixed parameter choice gives the following static band plot, which is useful for printed documentation and direct comparison:" -> "固定参数后得到下列静态能带图，便于打印和直接比较：",
  "Define the graphene path and open the archived interactive band panel:" -> "定义石墨烯路径并打开归档版交互能带面板：",
  "Fix the archived parameters to obtain a reproducible static band plot:" -> "固定归档版参数，得到可复现的静态能带图：",

  "Read gray magnetic space group 191 and inspect representative unitary and antiunitary records. Each record is {label, point matrix, translation, antiunitary flag}:" -> "读取第 191 号灰磁空间群，并查看代表性的幺正、反幺正记录。每条记录为 {标签, 点操作矩阵, 平移, 反幺正标志}：",
  "Display the complete Grid. OrbitalID is the row and column index used by both the Hamiltonian and every full representation matrix:" -> "显示完整 Grid。OrbitalID 是 Hamiltonian 和每个完整表示矩阵共同使用的行、列编号：",
  "Display the first nonzero-distance shell. The table gives the source site, translated target site, displacement, and operation that produces every directed bond:" -> "显示第一个非零距离键层。表格给出每条有向键的起点、平移后终点、位移和生成它的操作：",

  "Continuous C-infinity Symmetry with initfromrep" -> "使用 initfromrep 处理连续 C-infinity 对称性",
  "Consider two local states that acquire different phases under a rotation by theta. A hopping between them is allowed only when their SO(2) weights are compatible. initfromrep accepts the continuous symmetry as a one-parameter matrix D(theta), differentiates it at theta=0, and applies the resulting generator constraint before the remaining finite symmetries." -> "考虑两个在转角 theta 下取得不同相位的局域态。只有 SO(2) 权相容时，它们之间才允许跃迁。initfromrep 把连续对称性接收为单参数矩阵 D(theta)，在 theta=0 处求导，并先施加所得生成元约束，再处理其余有限对称性。",
  "MagneticTB accepts one pure-internal one-parameter SO(2) generator through initfromrep. The user supplies D(theta); initfromrep differentiates it exactly to obtain the Hermitian generator J and applies the continuous linear constraint before the finite discrete constraints. C_infinity never acts on atomic positions, never participates in finite-group closure, and never creates bond orbits." -> "MagneticTB 通过 initfromrep 接受一个纯内部的一参数 SO(2) 生成元。用户输入 D(theta)，initfromrep 对它精确求导得到厄米生成元 J，并先施加连续线性约束，再施加有限离散约束。C_infinity 不作用于原子位置，不参与有限群闭包，也不生成键轨道。",
  "A minimal two-state continuous representation" -> "最小二态连续表示",
  "The continuous operation is internal: it acts on the local states but does not move atoms or create new bonds. The two states use D(theta)=diag(exp(-i theta/2),exp(i theta/2))." -> "连续操作是内部操作：它作用于局域态，但不移动原子，也不产生新键。两个态使用 D(theta)=diag(exp(-i theta/2),exp(i theta/2))。",
  "The following source is self-contained. Its first model has only the finite identity and D(theta)=diag(exp(-i theta/2),exp(i theta/2)). The cubic lattice is written with the physical lattice constant a even though the continuous operation is internal." -> "下面的源代码可以独立运行。第一个模型的有限部分只有单位元，且 D(theta)=diag(exp(-i theta/2),exp(i theta/2))。虽然连续操作是内部操作，立方晶格仍显式写出物理晶格常数 a。",
  "The off-diagonal onsite matrix elements vanish because their SO(2) weights differ. The displayed finite representation contains only E; C_infinity is stored separately as a continuous constraint." -> "两个态的 SO(2) 权不同，因此 onsite 非对角矩阵元消失。显示的有限表示只有 E；C_infinity 作为连续约束单独保存。",
  "The off-diagonal onsite elements vanish because the SO(2) weights differ. The finite representation contains only E; C_infinity is stored separately as a continuous constraint." -> "两个态的 SO(2) 权不同，因此 onsite 非对角矩阵元消失。有限表示只有 E；C_infinity 作为连续约束单独保存。",
  "How the continuous constraint is applied" -> "连续约束如何施加",
  "Differentiating D(theta) at theta=0 gives the generator J on each site, and the Hamiltonian must satisfy [J,H(k)]=0. Finite unitary and antiunitary operations are then checked separately." -> "在 theta=0 处对 D(theta) 求导可得到每个格点上的生成元 J，Hamiltonian 必须满足 [J,H(k)]=0；随后分别检查有限幺正和反幺正操作。",
  "DirectProduct and Induced modes" -> "DirectProduct 与 Induced 模式",
  "This minimal example supplies the complete two-state matrix directly. For a multi-site induced construction, use RepresentationMode -> \"Induced\" with a reference-site stabilizer and its matrices, as shown by the executable option example on the initfromrep page." -> "这个最小例子直接给出完整二态矩阵。多格点诱导构造应使用 RepresentationMode -> \"Induced\"，并提供参考格点稳定子及其矩阵；initfromrep 帮助页的可执行选项例子给出了完整写法。",
  "14.1.2.3.L.1 at the 000 point" -> "14.1.2.3.L.1 的 000 点",
  "The optional SpinLayerCorepresentations paclet supplies the ordered eight finite representatives and their exact double-valued matrices through public functions. At the 000 seed the finite spatial action generates {0,0,0} and {1/2,1/2,0}; the site stabilizer is read from the same records rather than typed by hand." -> "可选的 SpinLayerCorepresentations paclet 通过公开函数提供有序的八个有限代表元及其精确双值矩阵。有限空间作用把 000 种子生成 {0,0,0} 和 {1/2,1/2,0}；位置稳定子也从同一批记录计算，不手工录入。",
  "Four-band DirectProduct branch" -> "四带 DirectProduct 分支",
  "Each of the two sites carries the two-component m=plus/minus one-half local space. The complete finite matrices and both continuous site matrices are supplied explicitly, so the full Hamiltonian has four bands." -> "两个格点各带一个 m=正负二分之一的二分量局域空间。完整有限矩阵和两个格点的连续矩阵均显式输入，因此完整 Hamiltonian 有四条能带。",
  "Two-band Induced branch" -> "两带 Induced 分支",
  "Only the one-dimensional reference-site representation on the site stabilizer is supplied. The antiunitary and unitary coset operations induce the second site. Its continuous weight is the opposite sign, so the two-band Hamiltonian is diagonal and the first inter-site shell is forbidden by C-infinity." -> "这里只输入位置稳定子上的一维参考点表示，幺正和反幺正陪集操作诱导出第二个格点。第二个格点的连续权符号相反，因此两带 Hamiltonian 为对角形式，第一层格点间跃迁被 C-infinity 禁止。",
  "Validation scope" -> "验证范围",

  "General MagneticTB Examples" -> "MagneticTB 综合范例",
  "This tutorial migrates every physical workflow from Examples/GeneralExamples.nb to the MagneticTB 2.0 API. Each model is initialized explicitly, every Hamiltonian is computed by the current implementation, and the displayed matrices, files, and plots are generated from the same inputs." -> "本教程把 Examples/GeneralExamples.nb 中的全部物理流程迁移到 MagneticTB 2.0 API。每个模型都显式初始化，每个 Hamiltonian 都由当前实现计算，显示的矩阵、文件和图形均来自同一组输入。",
  "Graphene" -> "石墨烯",
  StringJoin[
    "Start with graphene, retaining one pz orbital on each carbon atom. The complete input ",
    "below gives the two-band Hamiltonian with onsite, nearest-neighbor, and second-neighbor ",
    "terms:"
  ] -> StringJoin[
    "先考虑石墨烯，每个碳原子保留一个 pz 轨道。下面给出完整输入，构造包含 ",
    "onsite、最近邻和第二近邻项的两带 Hamiltonian："
  ],
  "Three-band MoS2 model" -> "三带 MoS2 模型",
  StringJoin[
    "For monolayer MoS2, keep the three local orbitals dz2, dxy, and dx2-y2. The ",
    "substitution after symham converts the parameters to the convention used in the ",
    "MagneticTB 1.0 example. Assign moS2Parameters only after this conversion, so each value ",
    "corresponds to the intended matrix element."
  ] -> StringJoin[
    "对于单层 MoS2，保留 dz2、dxy 和 dx2-y2 三个局域轨道。symham ",
    "后的替换把参数换成 MagneticTB 1.0 例子使用的约定，再代入 ",
    "moS2Parameters。这里要注意参数与矩阵元的对应关系，不能只按相同的参数名称直接赋值。"
  ],
  StringJoin[
    "Reference: G.-B. Liu, W.-Y. Shan, Y. Yao, W. Yao, and D. Xiao, \"Three-band ",
    "tight-binding model for monolayers of group-VIB transition metal dichalcogenides,\" ",
    "Physical Review B 88, 085433 (2013), https://doi.org/10.1103/PhysRevB.88.085433. The ",
    "example below uses the three-band model and parameter values from the MagneticTB 1.0 ",
    "example; the parameter conversion does not change its numerical Hamiltonian."
  ] -> StringJoin[
    "三带模型可参考 G.-B. Liu 等，\"Three-band tight-binding model ",
    "for monolayers of group-VIB transition metal ",
    "dichalcogenides,\" Physical Review B 88, 085433 ",
    "(2013)，https://doi.org/10.1103/PhysRevB.88.085433。下面",
    "沿用 MagneticTB 1.0 例子的模型和参数值，参数换写不改变其数值 Hamiltonian。"
  ],
  "The CartesianCoordinates option changes only the momentum convention. Here is the actual Cartesian form, not a dimension summary:" -> "CartesianCoordinates 只改变动量约定。下面直接显示实际笛卡尔形式，而不是维数摘要：",
  "banddata writes the eigenvalues used for plotting. Use a portable temporary path, then inspect the exported numerical rows directly:" -> "banddata 写出用于绘图的本征值。这里使用可移植临时路径，并直接检查导出的数值行：",
  "Magnetic C3 Weyl point and Chern diagnostic" -> "磁性 C3 Weyl 点与 Chern 诊断",
  StringJoin[
    "For magnetic space group 143.3, the local spinor pair gives a four-band Hamiltonian. As ",
    "in the MagneticTB 1.0 example, select the two-band block with row and column indices ",
    "{1,4}. The replacement c3OriginalParameterCoordinates puts its parameters in the 1.0 ",
    "convention before numerical values are assigned."
  ] -> StringJoin[
    "磁空间群 143.3 配合下面的局域旋量基，可以得到四带 Hamiltonian。与 ",
    "MagneticTB 1.0 例子一样，这里取行列编号为 {1,4} 的两带子块。代入数值前，先用 ",
    "c3OriginalParameterCoordinates 把参数换成 1.0 例子的约定。"
  ],
  "To reproduce the original Chern diagnostic, first keep the total cubic Taylor expansion around Gamma. The occupied eigenvector is then parallel-transported around 10-step latitude circles at 1001 polar angles, and the endpoint phase is reduced modulo 2 Pi. The final ListPlot is generated by the same discretization used in the original notebook." -> "为了重现原来的 Chern 诊断，先保留 Gamma 附近的总三次 Taylor 展开。随后在 1001 个极角上，沿每条分成 10 步的纬线环平行输运占据态本征矢，并把终点相位对 2 Pi 取模。最后的 ListPlot 使用与原 notebook 完全相同的离散方式生成。",
  "Magnetic cubic nodal line" -> "磁性三次节点线",
  StringJoin[
    "Next use magnetic space group 184.196 with the same local spinor basis. The operations ",
    "numbered {2,7,13} in the MagneticTB 1.0 example generate this 24-operation group; init ",
    "below receives the full group. Again select the {1,4} two-band block, and use ",
    "nodalOriginalParameterCoordinates to put its parameters in the 1.0 convention."
  ] -> StringJoin[
    "下面改用磁空间群 184.196，局域旋量基保持不变。MagneticTB 1.0 例子中编号为 ",
    "{2,7,13} 的操作可以生成这个含 24 个操作的群，下面直接把完整群传给 init。同样取 ",
    "{1,4} 两带子块，并用 nodalOriginalParameterCoordinates 换成 ",
    "1.0 例子的参数约定。"
  ],
  "C4 topological-insulator model" -> "C4 拓扑绝缘体模型",
  StringJoin[
    "This model has two inequivalent Wyckoff positions and two local orbitals. Keep the ",
    "indicated neighbor hoppings and make the parameter substitutions shown below, as in the ",
    "MagneticTB 1.0 example:"
  ] -> StringJoin[
    "这个模型包含两个不等价 Wyckoff 位置和两个局域轨道。按照 MagneticTB 1.0 ",
    "例子，保留下面列出的近邻跃迁，并按所给的关系替换参数："
  ],
  "Wyckoff, layer-group, and rod-group databases" -> "Wyckoff、层群和杆群数据库",
  "The database display functions return the actual ordered records used for model construction. showMSGWyckoff displays positions and allowed magnetic moments; mlgop and mrgop return magnetic layer- and rod-group operations." -> "数据库显示函数返回建模实际使用的有序记录。showMSGWyckoff 显示位置和允许磁矩，mlgop 与 mrgop 返回磁层群和磁杆群操作。",
  "API changes from the historical notebook" -> "相对历史 notebook 的 API 变化",
  "The 2.0 migration removes symmetryset from symham. A model is initialized once with its complete physical symmetry, and every later function stays within that compiled model. To construct a deliberately symmetry-broken model, call brokenSymmetryInitRules, inspect the returned complete init rules, and then evaluate init@@rules explicitly." -> "2.0 迁移从 symham 删除 symmetryset。模型只用完整物理对称性初始化一次，后续函数都限定在该已编译模型内。若要构造主动破缺对称性的模型，应调用 brokenSymmetryInitRules，检查返回的完整 init 规则，再显式计算 init@@rules。",

  "Getting Started with MagneticTB" -> "MagneticTB 入门",
  "MagneticTB starts from a lattice, Wyckoff seed, symmetry operations, and local orbitals. It expands symmetry-related sites, constructs the full representation, searches complete bond shells, solves the allowed hopping matrices, and assembles the Bloch Hamiltonian. This tutorial follows that physical order with graphene and shows every main result directly." -> "MagneticTB 从晶格、Wyckoff 种子、对称操作和局域轨道出发，展开对称相关格点，构造完整表示，搜索完整键层，求解允许的跃迁矩阵并组装 Bloch Hamiltonian。本教程按这一物理顺序贯穿石墨烯模型，直接显示每一步主要结果。",
  "Install, load, and open the documentation" -> "安装、加载与打开帮助",
  "Install, update, uninstall, and load" -> "安装、更新、卸载与加载",
  StringJoin[
    "Obtain a local MagneticTB .paclet file, either supplied to you or built from the ",
    "source. The installation commands below use this file; they do not download a release ",
    "from the Internet."
  ] -> StringJoin[
    "安装前需要取得 MagneticTB 的 .paclet ",
    "文件，也可以从源码构建。下面使用本地文件安装，不会从网上自动下载。"
  ],
  "Install" -> "安装",
  StringJoin[
    "Install the .paclet file with the following command. Replace the path and filename with ",
    "those of your file:"
  ] -> StringJoin[
    "取得 .paclet 文件后，在 Mathematica ",
    "中运行下面的命令，其中路径和文件名要换成实际文件："
  ],
  StringJoin[
    "After installation, quit Mathematica completely and open it again before loading ",
    "MagneticTB or opening its help. This also removes any definitions left by the ",
    "previously loaded version."
  ] -> StringJoin[
    "安装后，请完全退出 Mathematica 再重新打开，然后加载 MagneticTB ",
    "或打开帮助文档。这样也可以清除先前版本在内核中留下的定义。"
  ],
  "Update" -> "更新",
  StringJoin[
    "To update, install the newer .paclet file in the same way, then quit and reopen ",
    "Mathematica. If PacletInstall reports samevers, that version is already installed. To ",
    "replace it with another file of the same version, uninstall that version first and then ",
    "install the new file."
  ] -> StringJoin[
    "更新时，按同样的方法安装新的 .paclet 文件，再退出并重新打开 Mathematica。如果出现 ",
    "PacletInstall::samevers，说明相同版本已经安装；要换成同版本的另一个文件，需要先卸",
    "载该版本，再重新安装。"
  ],
  "Uninstall or reinstall the same version" -> "卸载或重装同一版本",
  "To uninstall MagneticTB, run the following command, then quit and reopen Mathematica:" -> "如果不再使用 MagneticTB，可以运行下面的命令卸载，再退出并重新打开 Mathematica：",
  "Load and open the documentation" -> "加载并打开帮助文档",
  StringJoin[
    "After reopening Mathematica, load the package with Needs. To open the help, search for ",
    "MagneticTB in the Documentation Center, or run the command below to open its home page. ",
    "The functions from version 1.0 are loaded with MagneticTBOld`; use a separate kernel ",
    "for that version."
  ] -> StringJoin[
    "重新打开 Mathematica 后，使用 Needs 加载程序包。在帮助中心搜索 ",
    "MagneticTB，或运行下面的命令，都可以打开帮助主页。1.0 版函数使用 ",
    "MagneticTBOld` 加载，需要与 2.0 版分别在不同内核中使用。"
  ],
  "Install a built MagneticTB 2.x paclet archive with PacletInstall. Installing a newer archive with the same paclet name updates the installed version. Quit the kernel before loading a newly installed build." -> "用 PacletInstall 安装构建好的 MagneticTB 2.x paclet。安装同名的更新版本会更新已安装版本；加载新构建前应退出当前内核。",
  "Search for MagneticTB in the Documentation Center, or open the guide URI below. MagneticTBOld` is included for compatibility but is fully isolated; use a second fresh kernel for old/new comparisons." -> "可在帮助中心搜索 MagneticTB，也可直接打开下列指南 URI。兼容后端 MagneticTBOld` 与新版完全隔离；比较新旧结果时请使用另一个全新内核。",
  "MagneticTBOld` is shipped as a fully isolated compatibility backend for running original notebooks and comparing old and new results. Load it in a separate kernel." -> "MagneticTBOld` 是与新版完全隔离的兼容后端，用于继续运行原有 notebook 并比较新旧结果。请在单独的内核中加载它。",
  "The five physical inputs" -> "五类物理输入",
  "A model is determined by five pieces of physical data. lattice contains primitive vectors as rows; lattpar assigns their symbolic parameters; wyckoffposition contains one seed position and magnetic moment for each Wyckoff orbit; symminformation supplies ordered finite operations; basisFunctions lists the local orbitals attached to each orbit." -> "模型由五类物理数据确定：lattice 以行为原胞基矢；lattpar 给符号晶格参数赋值；wyckoffposition 为每个 Wyckoff 轨道给出一个种子位置和磁矩；symminformation 提供有序有限操作；basisFunctions 列出附着在各轨道上的局域轨道。",
  "For graphene, use the database's hexagonal convention. The carbon seed at {1/3,2/3,0} expands to two equivalent sites. A large c separates repeated layers, and one pz orbital is placed on every generated carbon." -> "石墨烯使用数据库的六方约定。碳原子种子 {1/3,2/3,0} 展开为两个等价格点；较大的 c 分隔周期重复层，每个生成的碳原子放置一个 pz 轨道。",
  "msgop[gray[191]] returns the complete gray group. Therefore GenerateSymmetryGroup -> False means that init uses those records directly; it does not mean that every record is a generator. The completed multiplication table has the following compact generating set:" -> "msgop[gray[191]] 返回完整灰群。因此 GenerateSymmetryGroup -> False 表示 init 直接使用这些记录，并不表示每条记录都是生成元。完整乘法表的紧凑生成集如下：",
  "Check the orbital order" -> "检查轨道顺序",
  "orbitalTable[] gives the row and column order shared by every Hamiltonian and full representation matrix. FractionalPosition uses primitive-lattice coordinates; CartesianPosition is evaluated after lattpar. Spinless states use SpinState -> None so the table remains extensible to spin groups." -> "orbitalTable[] 给出所有 Hamiltonian 与完整表示矩阵共同使用的行列顺序。FractionalPosition 使用原胞分数坐标，CartesianPosition 在代入 lattpar 后计算；无自旋态写成 SpinState -> None，以便将来扩展到自旋群。",
  "Inspect prepared bond shells" -> "检查已准备键层",
  "Shell 1 is the zero-distance onsite shell. Shell 2 is the shortest nonzero shell and contains the three nearest-neighbour carbon-carbon displacements in both directions. showbonds reads data prepared by init; it never searches again or calls symham." -> "第 1 键层是零距离 onsite 层。第 2 键层是最短非零键层，包含三条碳—碳最近邻位移及其反向。showbonds 只读取 init 准备的数据，不重新搜索，也不调用 symham。",
  "Solve the Hamiltonian shell by shell" -> "逐键层求解 Hamiltonian",
  StringJoin[
    "First calculate the onsite term with symham[1]. The two carbon atoms are equivalent, so ",
    "both diagonal entries contain the same real parameter e1:"
  ] -> "先用 symham[1] 求 onsite 项。两个碳原子是等价的，所以对角元含有相同的实参数 e1：",
  StringJoin[
    "Next calculate the nearest-neighbor hopping with symham[2]. The three phase factors in ",
    "the upper-right entry come from the three nearest-neighbor bonds. For real t1, the ",
    "lower-left entry is its complex conjugate:"
  ] -> StringJoin[
    "再用 symham[2] 求最近邻跃迁。右上角的三个相位因子来自三条最近邻键；t1 ",
    "为实数时，左下角是它的复共轭："
  ],
  StringJoin[
    "symham[3] gives the second-neighbor hopping. These bonds connect carbon atoms on the ",
    "same sublattice, so their six phase factors appear in the diagonal entries:"
  ] -> "symham[3] 给出第二近邻跃迁。这些键连接同一子晶格上的碳原子，所以六个相位因子出现在对角元中：",
  "The working Hamiltonian is the ordinary sum of the prepared shell contributions. init is not called again, and repeated symham calls with identical options read the solved-shell cache." -> "工作 Hamiltonian 是已准备键层贡献的普通求和。这里不会再次调用 init；使用相同选项重复调用 symham 时直接读取已求解键层缓存。",
  "Explore the bands" -> "查看能带",
  "A path is a list of reciprocal-fractional line segments and endpoint labels. bandManipulate multiplies the coordinates by 2 Pi internally. Every Hamiltonian symbol other than kx, ky, and kz becomes a slider; all sliders start at zero and no hidden energy shift is inserted." -> "path 是倒空间分数坐标线段及端点标签的列表。bandManipulate 在内部把坐标乘以 2 Pi。Hamiltonian 中除 kx、ky、kz 外的每个符号都成为滑块；所有滑块从零开始，不插入隐藏能量平移。",
  StringJoin[
    "To draw a band plot with fixed parameters, assign their values as below. The ",
    "second-neighbor term r1 breaks particle-hole symmetry, but the Dirac crossing at K ",
    "remains:"
  ] -> StringJoin[
    "也可以给参数赋值，直接画出能带。第二近邻参数 r1 破坏粒子空穴对称性，但不会打开 K 点处的 ",
    "Dirac 交叉："
  ],
  "See how symmetry acts" -> "查看对称性如何作用",
  "showSymmetryRepresentations[] displays the matrices of the compact generators. The antiunitary flag already belongs to the operation, so only the matrix is displayed. Unitary operations obey U H(k) U^dagger = H(g k); antiunitary operations obey U Conjugate[H(k)] U^dagger = H(g k)." -> "showSymmetryRepresentations[] 显示紧凑生成元矩阵。反幺正标志已属于操作本身，所以这里只显示矩阵。幺正操作满足 U H(k) U^dagger = H(g k)；反幺正操作满足 U Conjugate[H(k)] U^dagger = H(g k)。",
  "Evaluation boundaries and next models" -> "计算边界与后续模型",
  "InitialBondShells defaults to 10; this tutorial requests 3. symham[4] must therefore fail and ask for a new init with a larger value. It must not extend the session or switch algorithms silently. A failed init clears the old session instead of leaving stale model data." -> "InitialBondShells 默认为 10，本教程显式请求 3。因而 symham[4] 必须失败，并要求用更大值重新 init；它不能静默扩展会话或切换算法。init 失败时会清空旧会话，不留下过期模型数据。",

  "Validated MagneticTB Models" -> "经验证的 MagneticTB 模型",
  "This page turns the main regression models into user-facing applications. Each section starts from a complete physical input and displays the Hamiltonian, orbital or hopping interpretation, symmetry matrices, and bands. Parameter counts and residuals are checked by the documentation regression but are not used as substitutes for the physical result." -> "本页把主要回归模型转化为面向用户的应用。每节都从完整物理输入开始，显示 Hamiltonian、轨道或跃迁解释、对称矩阵和能带。参数数目与残差在文档回归中检查，但不代替物理结果。",
  "CsCl: a 1 by 3 rectangular s-p hopping block" -> "CsCl：1×3 矩形 s-p 跃迁块",
  "The corner site carries one s orbital and the body-centred site carries px, py, pz. The off-diagonal row below is therefore the actual rectangular hopping block; it is not required to be a square or Hermitian matrix by itself." -> "角点格位带一个 s 轨道，体心格位带 px、py、pz。所以下面的非对角行就是实际矩形跃迁块；它本身不必是方阵，也不必厄米。",
  "FeSe: two Wyckoff orbits with d and p orbitals" -> "FeSe：含 d、p 轨道的两个 Wyckoff 轨道",
  "This model places dxz,dyz on the first orbit and px,py,pz on the second. The orbital table fixes the row/column order; the displayed 3 by 3 corner makes the symbolic onsite and hopping structure readable without hiding the full ten-band calculation." -> "该模型在第一轨道放置 dxz、dyz，在第二轨道放置 px、py、pz。轨道表固定行列顺序；显示 3×3 角块既保留完整十带计算，又使符号 onsite 和跃迁结构可读。",
  "Spinful square px/py model" -> "有自旋方格 px/py 模型",
  "The local ordering is px-up, px-down, py-up, py-down. The onsite Hamiltonian and full generator matrices show how orbital and spin transformations are combined before the bond constraints are solved." -> "局域顺序为 px-up、px-down、py-up、py-down。onsite Hamiltonian 与完整生成元矩阵展示在求解键约束前如何组合轨道与自旋变换。",
  "Cubic p orbitals: DirectProduct and Induced representations" -> "立方 p 轨道：DirectProduct 与 Induced 表示",
  "Collinear discrete spin-space group" -> "共线离散自旋空间群",
  "C4spin is a pure internal spin rotation. tauHalfT combines a half translation with an antiunitary spin action. GenerateSymmetryGroup closes these two generators before the physical action and Hamiltonian constraints are compiled." -> "C4spin 是纯内部自旋旋转，tauHalfT 把半平移与反幺正自旋作用组合起来。GenerateSymmetryGroup 先闭合两个生成元，再编译物理作用与 Hamiltonian 约束。",
  "Noncollinear octahedral spin-space group" -> "非共线八面体自旋空间群",
  "Two noncommuting pure-spin rotations generate the finite noncollinear group. The representation output displays the spinor matrices actually used to constrain the onsite and hopping terms." -> "两个不对易的纯自旋旋转生成有限非共线群。表示输出显示实际用于约束 onsite 和跃迁项的旋量矩阵。",
  "E plus PT: antiunitary induction between two sites" -> "E+PT：两格点间的反幺正诱导",
  "E alone cannot move the seed at x=1/4 to its partner at x=3/4. PT is therefore the required antiunitary coset representative. Its antiunitary property comes from symminformation; the displayed full matrix contains no separate K marker. Hamiltonian covariance uses U Conjugate[H(k)] U-dagger." -> "E 单独不能把 x=1/4 的种子移到 x=3/4 的伙伴点，因此 PT 是必需的反幺正陪集代表。反幺正属性来自 symminformation；显示的完整矩阵不另带 K 标记。Hamiltonian 协变关系使用 U Conjugate[H(k)] U-dagger。",
  "C3 times Z2^T over C3: a complex local character" -> "C3×Z2^T/C3：复局域特征标",
  "The reference-site stabilizer is C3 with the one-dimensional character {1,omega,omega^2}. An antiunitary coset maps the reference site to its partner. Because the target coset representative is antiunitary, the local character is complex conjugated during induction; the two diagonal C3 entries are omega and omega^2." -> "参考点稳定子为 C3，其一维特征标是 {1,omega,omega^2}。反幺正陪集把参考点映到伙伴点。由于目标陪集代表为反幺正，诱导时局域特征标取复共轭，所以 C3 的两个对角元分别为 omega 与 omega^2。",
  "Further validated applications" -> "其他经验证应用",
  "Graphene, the three-band MoS2 model, magnetic Weyl and nodal-line models, and the C4 topological-insulator example are developed in General MagneticTB Examples. The ten database-driven DirectProduct/Induced pairs are developed in Wyckoff Database Models in Seven Crystal Systems." -> "石墨烯、三带 MoS2、磁性 Weyl 与节点线以及 C4 拓扑绝缘体范例见《MagneticTB 综合范例》。十组数据库 DirectProduct/Induced 对照见《七晶系 Wyckoff 数据库模型》。",
  "Legacy backend" -> "旧版后端",
  "The old implementation must be evaluated in a separate fresh kernel. See Using the Legacy Backend for a complete initold, symhamold, and bandManipulateold example with outputs." -> "旧实现必须在独立全新内核中计算。完整且带实际输出的 initold、symhamold、bandManipulateold 范例见《使用旧版后端》。",

  "Wyckoff Database Models in Seven Crystal Systems" -> "七晶系 Wyckoff 数据库模型",
  "These ten examples cover all seven crystal systems and the exact BNS/Wyckoff entries used by the DirectProduct/Induced regression. Every seed is read from the paclet's formal wyckoffMSG.mx database. The helper below contains the complete init, symham, representation, bandManipulate, and fixed-parameter bandplot workflow used by every case." -> "这十个范例覆盖七大晶系，并对应 DirectProduct/Induced 回归使用的精确 BNS/Wyckoff 条目。每个种子都从 paclet 正式 wyckoffMSG.mx 数据库读取。下面的辅助代码包含每个案例完整的 init、symham、表示、bandManipulate 和固定参数 bandplot 流程。",
  "Database reader and complete model runner" -> "数据库读取器与完整模型运行函数",
  "The seed below is read from the formal wyckoffMSG.mx entry at evaluation time. The magnetic group requires an antiunitary coset to complete the two-site orbit." -> "下面的种子在计算时从正式 wyckoffMSG.mx 条目读取。该磁群必须使用反幺正陪集才能补全两格点轨道。",
  "DirectProduct representation" -> "DirectProduct 表示",
  "Induced representation" -> "Induced 表示",
  "Triclinic \[LongDash] BNS 1.3-a" -> "三斜晶系 \[LongDash] BNS 1.3-a",
  "Triclinic \[LongDash] BNS 2.6-i" -> "三斜晶系 \[LongDash] BNS 2.6-i",
  "Monoclinic \[LongDash] BNS 3.3-e" -> "单斜晶系 \[LongDash] BNS 3.3-e",
  "Monoclinic \[LongDash] BNS 3.4-a" -> "单斜晶系 \[LongDash] BNS 3.4-a",
  "Orthorhombic \[LongDash] BNS 16.3-q" -> "正交晶系 \[LongDash] BNS 16.3-q",
  "Orthorhombic \[LongDash] BNS 16.3-r" -> "正交晶系 \[LongDash] BNS 16.3-r",
  "Tetragonal \[LongDash] BNS 75.3-c" -> "四方晶系 \[LongDash] BNS 75.3-c",
  "Trigonal \[LongDash] BNS 143.3-a" -> "三方晶系 \[LongDash] BNS 143.3-a",
  "Hexagonal \[LongDash] BNS 168.111-b" -> "六方晶系 \[LongDash] BNS 168.111-b",
  "Cubic \[LongDash] BNS 195.3-a" -> "立方晶系 \[LongDash] BNS 195.3-a",
  "What was certified" -> "已认证内容",
  "For every pair, the documentation regression checks the real Hamiltonian parameter subspaces in both directions, aligned bands at five momenta, strict representation/corepresentation multiplication, six generic-k covariance reports, and all non-Gamma little groups on the {0,Pi/2,Pi}^3 grid with the reciprocal-lattice sewing matrix. Gamma-only checks are not used as evidence of full space-group covariance." -> "对每一对模型，文档回归双向检查实 Hamiltonian 参数子空间、五个动量点的对齐能带、严格表示/共表示乘法、六个一般 k 点的协变报告，以及 {0,Pi/2,Pi}^3 网格上全部非 Gamma 小群和倒格矢 sewing 矩阵。只检查 Gamma 点不能作为完整空间群协变性的证据。",

  "Read a compact physical summary from the current session. This inspection does not compile or solve another shell:" -> "从当前会话读取紧凑物理摘要；该检查不会编译或求解其他键层：",
  "The user-supplied labels determine the displayed orbital order; no basis-function transformation is inferred:" -> "用户提供的标签决定显示的轨道顺序；程序不会推断任何基函数变换：",
  "Initialize the two-dimensional px/py representation of the tetragonal {E,C4z,C2z,C4z^-1} group. The outer list enumerates Wyckoff orbits and the inner list enumerates the complete ordered finite group:" -> "用四方群 {E,C4z,C2z,C4z^-1} 的二维 px/py 显式表示初始化模型。repinformation 的外层列表枚举 Wyckoff 轨道，内层列表按照完整有限群的固定顺序给出矩阵：",
  "The one finite operation is represented by the supplied two-dimensional identity matrix:" -> "唯一有限操作由用户提供的二维单位矩阵表示：",
  "The corresponding onsite Hamiltonian is the most general Hermitian 2 by 2 matrix allowed by this representation:" -> "相应 onsite Hamiltonian 是该表示允许的最一般 2×2 厄米矩阵：",
  "InitialBondShells controls the exact range prepared by initialization. This explicit model prepares only its onsite shell and displays that Hamiltonian:" -> "InitialBondShells 精确控制初始化准备的范围。该显式模型只准备 onsite 键层，并显示对应 Hamiltonian：",
  "Initialize the archived graphene model and display its onsite Hamiltonian. This example is evaluated by a separate legacy-only kernel:" -> "在独立的纯旧版内核中初始化归档石墨烯模型，并显示其 onsite Hamiltonian：",
  "Initialize graphene from a complete ordered gray-group list. GenerateSymmetryGroup -> False uses those records directly; it does not mean that every operation is a generator:" -> "用完整有序灰群列表初始化石墨烯。GenerateSymmetryGroup -> False 表示直接使用这些记录，不表示每个操作都是生成元：",
  "The completed 48-element multiplication table has the following compact generating set:" -> "闭合后的 48 元乘法表具有如下紧凑生成集：",
  "Inspect the actual Hamiltonian ordering. The two rows are the symmetry-related carbon sites generated from the single Wyckoff seed:" -> "检查实际 Hamiltonian 顺序。两行对应由单个 Wyckoff 种子生成的两个对称相关碳格点：",
  "GenerateSymmetryGroup -> True accepts a generating set. Here the identity is expanded trivially, and the actual onsite Hamiltonian is shown:" -> "GenerateSymmetryGroup -> True 接受生成集。这里单位元只作平凡闭合，并直接显示实际 onsite Hamiltonian：",
  "Construct onsite plus nearest-neighbour hopping and display the actual archived Bloch Hamiltonian:" -> "构造 onsite 与最近邻跃迁之和，并显示实际归档 Bloch Hamiltonian：",
  "Initialize the complete graphene model used below. symham never performs initialization implicitly:" -> "先初始化下文使用的完整石墨烯模型；symham 绝不隐式执行初始化：",
  "Shell 1 is the onsite contribution. The same onsite parameter appears on the two equivalent carbon sites:" -> "第 1 键层是 onsite 贡献；两个等价碳格点上出现同一个 onsite 参数：",
  "Shell 2 is the nearest-neighbour hopping. Its off-diagonal entries contain the three symmetry-related Bloch phases:" -> "第 2 键层是最近邻跃迁；非对角元包含三个对称相关 Bloch 相位：",
  "Add onsite, nearest-neighbour, and next-nearest-neighbour shells to obtain the working Bloch Hamiltonian:" -> "把 onsite、最近邻和次近邻键层相加，得到工作 Bloch Hamiltonian：",
  "CartesianCoordinates -> True rewrites the same shell in Cartesian momentum coordinates defined by the reciprocal lattice:" -> "CartesianCoordinates -> True 把同一键层改写为倒格子定义的笛卡尔动量坐标：",
  StringJoin[
    "For graphene, give the hexagonal lattice, one representative carbon position, and its ",
    "pz orbital. The symmetry operations generate the second carbon atom. In version 1.0, ",
    "the function and option names below all end in old:"
  ] -> StringJoin[
    "构造石墨烯模型时，给出六方晶格、一个碳原子的代表位置和 pz 轨道，对称操作会生成另一个碳原子。1.0 ",
    "版的函数和选项名称末尾都带 old，完整输入如下："
  ],

  " — close concrete finite-group generators under an explicit product." -> " — 在显式乘法规则下闭合具体的有限群生成元。",
  " — compile canonical model data without installing a session." -> " — 编译规范模型数据，但不安装会话。",
  " — construct exact local point-operation matrices." -> " — 构造精确的局域点操作矩阵。",
  " — inspect an unconstrained directed-bond Hamiltonian." -> " — 查看未施加约束的有向键 Hamiltonian。",
  " — prepare Wannier centers and symmetry matrices for HR symmetrization." -> " — 为 HR 对称化准备 Wannier 中心与对称矩阵。",
  " — use the deprecated 2.x generator-search compatibility wrapper." -> " — 使用 2.x 中已弃用的生成元搜索兼容封装。",
  "All displays the complete numerical spectrum; {emin,emax} restricts only the visible energy interval." -> "All 显示完整数值能谱；{emin,emax} 只限制可见的能量区间。",
  "plotRange changes only the displayed vertical energy interval; the Hamiltonian and sampled eigenvalues are unchanged:" -> "plotRange 只改变显示的纵向能量范围；Hamiltonian 与采样得到的本征值不变：",
  "Compile a complete simple-cubic input and display its physical lattice, site orbit, exact representation matrices, and first bond class:" -> "编译一份完整的简单立方输入，并显示其物理晶格、格点轨道、精确表示矩阵和第一类键：",
  "compiles a canonical model association without installing it as the current session." -> "编译规范模型 Association，但不把它安装为当前会话。",
  "It prepares sites, the full exact representation, bond classes, and static directed-bond orbits; it does not solve a Hamiltonian shell and does not mutate CurrentModelSession." -> "它准备格点、完整精确表示、键分类和静态有向键轨道；它不求解 Hamiltonian 键层，也不修改 CurrentModelSession。",
  "Missing canonical fields are rejected; no partial model is returned:" -> "缺少规范字段的输入会被拒绝；函数不会返回不完整模型：",
  "This advanced programmatic entry point accepts the same physical information used by init, but with canonical string keys." -> "这一高级程序化入口接受与 init 相同的物理信息，但使用规范字符串键。",
  "Use this entry point when another program generates exact model data and needs a validated compiled value before deciding whether to call the public session API:" -> "当其他程序生成精确模型数据，并且需要先得到经过验证的编译结果再决定是否调用公开会话接口时，可使用该入口：",
  "Binary predicate used to decide whether two generated concrete elements are equal. It must return True or False for every pair." -> "用于判断两个已生成具体元素是否相等的二元谓词；对任意元素对都必须返回 True 或 False。",
  "closes concrete generator elements under a user-supplied multiplication law." -> "按照用户给定的乘法规则闭合具体生成元。",
  "Generate the additive cyclic group C4 from the element 1:" -> "由元素 1 生成加法循环群 C4：",
  "SameTest can provide a domain-specific exact equality predicate; here it reproduces the same C4 ordering:" -> "SameTest 可以提供与具体领域相符的精确相等判据；这里它得到相同的 C4 排序：",
  "The identity is always the first returned element, including for the trivial group." -> "单位元始终是返回列表中的第一个元素，平凡群也不例外。",
  "Use the operation-specific product and equality functions required by the physical symmetry type." -> "应使用与物理对称类型相符的专用乘法和相等判断函数。",
  "An explicit exact SameTest leaves the mathematical result unchanged:" -> "显式给出精确 SameTest 不会改变数学结果：",
  "Binary predicate used for exact equality of concrete group elements." -> "用于精确判断具体群元素是否相等的二元谓词。",
  "getGenerator is a deprecated compatibility wrapper retained for 2.x notebooks." -> "getGenerator 是为兼容 2.x notebook 而保留的弃用封装。",
  "Recover one generator of the complete additive C4 group:" -> "从完整的加法群 C4 中找回一个生成元：",
  "returns a greedy concrete generating subset of a complete finite group." -> "返回完整有限群的一个贪心具体生成子集。",
  "The canonical implementation works from the complete finite set and never guesses a different multiplication law." -> "规范实现直接处理完整有限集合，绝不猜测另一套乘法规则。",
  "For a supplied H(k), None uses zero centers, Automatic uses the current initialized orbital centers when available, and an explicit list supplies one fractional center per orbital." -> "对于给定的 H(k)，None 使用零中心；Automatic 在可用时采用当前已初始化的轨道中心；显式列表则为每条轨道给出一个分数坐标中心。",
  "For shell selection, choose the cached Hermitian or directed non-Hermitian symham solution. It is not used when extracting a supplied matrix H(k)." -> "按键层导出时，可选择已缓存的厄米或有向非厄米 symham 解；从给定矩阵 H(k) 提取时不使用此选项。",
  "For shell selection, choose the cached \"Iterative\" or \"Stacked\" null-space solution." -> "按键层导出时，选择已缓存的 \"Iterative\" 或 \"Stacked\" 零空间解。",
  "For shell selection, choose the cached \"None\", \"Basic\", or \"Full\" validation result." -> "按键层导出时，选择与 \"None\"、\"Basic\" 或 \"Full\" 对应的缓存验证结果。",
  "Integer at least 6 giving the number of decimal digits written for real and imaginary hopping components." -> "不小于 6 的整数，指定写出 hopping 实部和虚部时的小数位数。",
  "None returns and prints Wannier90 text; a file or directory path writes wannier90_hr.dat." -> "None 返回并打印 Wannier90 文本；文件或目录路径会写出 wannier90_hr.dat。",
  "Positive numerical tolerance used only when extracting integer real-space translations from exponential phases in a supplied H(k)." -> "正数容差，只在从给定 H(k) 的指数相位提取整数实空间平移时使用。",
  "The following complete export selects the cached Hermitian iterative Basic solution, uses current orbital centers, writes ten decimal digits, and applies the stated translation tolerance when phases need extraction:" -> "下面的完整导出选择已缓存的厄米、迭代、Basic 解，使用当前轨道中心，写出十位小数，并在需要从相位提取平移时采用给定容差：",
  "A complete, already closed, ordered finite magnetic or spin-space group. One additional pure-internal C_infty Association is allowed." -> "完整、已经闭合且有序的有限磁群或自旋空间群；另外允许一个纯内部的 C_infinity Association。",
  "Boolean compatibility flag. It does not change the representation, solver, or failure semantics." -> "布尔兼容选项；它不改变表示、求解器或失败语义。",
  "Display labels for the local representation coordinates. Automatic generates orb1, orb2, and so on." -> "局域表示坐标的显示标签；Automatic 依次生成 orb1、orb2 等标签。",
  "Exact representation matrices in exactly the finite-group order. Use the documented Discrete/Continuous Association for C_infty." -> "严格按照有限群顺序排列的精确表示矩阵；C_infinity 使用文档规定的 Discrete/Continuous Association。",
  "initfromrep deliberately rejects basisFunctions because exact ordered representation matrices, rather than transformed functions, define this input route:" -> "initfromrep 有意拒绝 basisFunctions，因为该输入途径由有序的精确表示矩阵而非变换后的函数定义：",
  "One fractional seed and magnetic moment per Wyckoff orbit; the orbit order must match repinformation and orbitalLabels." -> "每个 Wyckoff 轨道给出一个分数坐标种子和磁矩；轨道顺序必须与 repinformation 和 orbitalLabels 一致。",
  "Positive integer fixing the complete shell range prepared during initialization. symham never extends it." -> "正整数，用于固定初始化期间准备的完整键层范围；symham 绝不扩展该范围。",
  "Primitive direct-lattice vectors stored by rows. The numerical lattice obtained after lattpar must be real and nonsingular." -> "按行存储的原胞正格矢；代入 lattpar 后得到的数值晶格必须为实数且非奇异。",
  "Replacement rules assigning every symbolic lattice parameter used by lattice." -> "为 lattice 中出现的每个符号晶格参数赋值的替换规则。",
  "Use \"DirectProduct\" for full ordered matrices on every orbit, or \"Induced\" for reference-site stabilizer data." -> "每条轨道都有完整有序矩阵时使用 \"DirectProduct\"；从参考格点稳定子数据出发时使用 \"Induced\"。",
  "All six historical options remain accepted together. This complete explicit call displays the resulting onsite Hamiltonian inside the isolated old kernel:" -> "六个历史选项仍可同时使用。下面的完整显式调用在隔离的旧版内核中显示得到的 onsite Hamiltonian：",
  "Boolean legacy diagnostic flag retained unchanged for notebook compatibility." -> "为兼容旧 notebook 而原样保留的布尔诊断选项。",
  "Complete ordered legacy magnetic-group operation records." -> "完整有序的旧版磁群操作记录。",
  "Legacy list of fractional Wyckoff seeds and magnetic moments." -> "旧版分数 Wyckoff 种子与磁矩列表。",
  "Legacy local function basis, ordered exactly as wyckoffpositionold." -> "旧版局域函数基，其顺序与 wyckoffpositionold 完全一致。",
  "Legacy numerical substitutions for all symbolic lattice parameters." -> "旧版中对所有符号晶格参数的数值替换。",
  "Legacy primitive direct-lattice vectors stored by rows." -> "按行存储的旧版原胞正格矢。",
  "debugQ is a Boolean compatibility option; other values are rejected instead of being silently ignored:" -> "debugQ 是布尔兼容选项；其他取值会被拒绝，而不是被静默忽略：",
  "Unknown options fail before model compilation and leave no usable previous session:" -> "未知选项会在模型编译前失败，并且不会留下可继续使用的旧会话：",
  "A singular lattice cannot define the Cartesian point action:" -> "奇异晶格无法定义笛卡尔点操作：",
  "constructs exact unitary point-operation matrices on a local function basis." -> "在局域函数基上构造精确幺正点操作矩阵。",
  "operations use the traditional {label,R,t,flag} records; basis may contain scalar orbitals or two-component spinors." -> "operations 使用传统的 {label,R,t,flag} 记录；basis 可以包含标量轨道或二分量旋量。",
  "The identity and a twofold rotation about z act on the ordered {px,py} basis as follows:" -> "单位操作与绕 z 轴二重旋转在有序 {px,py} 基上的作用如下：",
  "The returned matrices follow the operation order and act on the supplied basis order." -> "返回矩阵遵循操作顺序，并作用于用户给定的基顺序。",
  "All keeps every translation; a nonnegative integer n keeps cells whose largest absolute translation component is at most n." -> "All 保留全部平移；非负整数 n 只保留最大平移分量绝对值不超过 n 的晶胞。",
  "Apply a zero-cell cutoff and a 10^-12 cleanup threshold after strictly parsing the complete file:" -> "严格解析完整文件后，只保留零晶胞并用 10^-12 阈值清除数值噪声：",
  "Nonnegative numerical threshold used to Chop small real and imaginary matrix entries after parsing." -> "非负数值阈值，用于在解析后 Chop 掉矩阵实部和虚部中的小量。",
  "KernelMethod and ValidationLevel select an already solved cache entry; they never trigger an implicit alternative solve. This Full iterative entry is solved explicitly before inspection:" -> "KernelMethod 与 ValidationLevel 选择已经求解的缓存条目；它们绝不会隐式触发另一种求解。下面先显式求解 Full 迭代条目，再查看其来源：",
  "Select provenance from the cached \"Iterative\" or \"Stacked\" solution; this function never runs a missing solver." -> "从已缓存的 \"Iterative\" 或 \"Stacked\" 解中选择来源；该函数绝不会运行缺失的求解器。",
  "Select provenance from the Hermitian shell solution or from a previously solved directed non-Hermitian shell." -> "选择厄米键层解或此前已求解的有向非厄米键层解的来源。",
  "Select the cached \"None\", \"Basic\", or \"Full\" solution matching the earlier symham call." -> "选择与此前 symham 调用相匹配的 \"None\"、\"Basic\" 或 \"Full\" 缓存解。",
  "Select the separately cached directed non-Hermitian solution. The table now traces the independent forward and reverse amplitudes rather than the Hermitian parameter:" -> "选择单独缓存的有向非厄米解。此时表格追踪相互独立的正向与反向振幅，而不是厄米参数：",
  "CartesianCoordinates -> True rewrites the archived nearest-neighbour shell with the legacy reciprocal lattice. symmetrysetold is retained only for historical notebooks; the full prepared set is used here:" -> "CartesianCoordinates -> True 使用旧版倒格子改写归档的最近邻键层。symmetrysetold 只为历史 notebook 保留；这里使用完整的已准备操作集合：",
  "False uses the legacy reciprocal-fractional momenta; True rewrites the result using the legacy reciprocal lattice." -> "False 使用旧版倒格分数动量；True 使用旧版倒格子改写结果。",
  "Legacy subset of symmetry-operation indices. This historical option is kept only in the isolated old backend." -> "旧版对称操作编号子集；该历史选项只保留在隔离的旧版后端中。",
  "A shell outside the range prepared by init fails; symham never extends the bond search implicitly:" -> "超出 init 已准备范围的键层会失败；symham 绝不隐式扩展键搜索：",
  "False uses reciprocal-fractional momentum components kx, ky, kz; True rewrites the result in Cartesian momentum coordinates." -> "False 使用倒格分数动量分量 kx、ky、kz；True 将结果改写为笛卡尔动量坐标。",
  "\"Iterative\" restricts the candidate null space constraint by constraint; \"Stacked\" solves one vertically assembled matrix. Neither silently falls back to the other." -> "\"Iterative\" 逐条约束候选零空间；\"Stacked\" 求解一个纵向拼接的矩阵。两者都不会静默回退到另一算法。",
  "\"None\" performs no mathematical certification, \"Basic\" checks residuals, and \"Full\" also checks independence and rank-nullity." -> "\"None\" 不进行数学认证；\"Basic\" 检查残差；\"Full\" 还检查独立性与秩－零化度关系。",
  "The Hermitian option is strictly Boolean:" -> "Hermitian 选项严格要求布尔值：",
  "True relates every directed hopping to its reverse and returns a Hermitian shell. False keeps reverse amplitudes independent." -> "True 把每条有向 hopping 与其反向项联系起来并返回厄米键层；False 保持反向振幅独立。",
  "All six options are shown explicitly here: the lattice, its numerical parameters, the Wyckoff seed, ordered symmetry records, local basis, and VASP ordering together determine the returned centers and matrices:" -> "这里显式给出全部六个选项：晶格、数值晶格参数、Wyckoff 种子、有序对称记录、局域基和 VASP 排序共同决定返回的中心与矩阵：",
  "Complete ordered symmetry records." -> "完整有序的对称操作记录。",
  "Local basis used to build Wannier centers and point matrices." -> "用于构造 Wannier 中心和点操作矩阵的局域基。",
  "Numerical replacement rules for the symbolic lattice." -> "符号晶格的数值替换规则。",
  "Only the documented VASP ordering convention is currently accepted." -> "目前只接受文档规定的 VASP 排序约定。",
  "Orbital ordering convention; VASP is the only supported value." -> "轨道排序约定；目前唯一支持的取值是 VASP。",
  "prepares Wannier centers, full symmetry matrices, and ordered symmetry records for HR symmetrization." -> "为 HR 对称化准备 Wannier 中心、完整对称矩阵和有序对称记录。",
  "Prepare the exact one-orbital symmetrization data for a simple-cubic model:" -> "为简单立方单轨道模型准备精确的对称化数据：",
  "Symbolic direct-lattice vectors stored by rows." -> "按行存储的符号正格矢。",
  "This compatibility entry point returns the native association keys wcc, DR, and symmetry and never loads the archived backend." -> "这一兼容入口返回原生 Association 键 wcc、DR 和 symmetry，并且绝不加载归档后端。",
  "Unsupported software conventions fail explicitly:" -> "不支持的软件约定会明确失败：",
  "Wyckoff seeds and magnetic moments in the same nesting used by init." -> "采用与 init 相同嵌套结构的 Wyckoff 种子与磁矩。",
  "constructs the unconstrained symbolic Hamiltonian for prepared bond shell n." -> "为已准备的第 n 个键层构造未约束的符号 Hamiltonian。",
  "For the same prepared shell, symham reduces these independent entries to the symmetry-allowed Hermitian subspace:" -> "对于同一个已准备键层，symham 把这些独立矩阵元约化到对称性允许的厄米子空间：",
  "Initialize a simple-cubic s-orbital model, then display the independent complex amplitude carried by every directed nearest-neighbour bond:" -> "初始化简单立方 s 轨道模型，然后显示每条有向最近邻键携带的独立复振幅：",
  "It is an advanced diagnostic; use symham for a physical symmetry-constrained model." -> "它是高级诊断函数；物理的对称约束模型应使用 symham。",
  "The shell must already have been prepared by init; the function never extends the bond search:" -> "该键层必须已经由 init 准备；函数绝不扩展键搜索：",
  "unsymham uses the bond geometry prepared by init but deliberately applies no symmetry or Hermiticity constraints." -> "unsymham 使用 init 准备的键几何，但有意不施加任何对称性或厄米性约束。",
  "Convert the finite records into the exact ordered symminformation association, generate the two sites, find the reference-site stabilizer, and append the one pure-internal C_infinity operation. These are preparation steps only; no Hamiltonian has yet been solved." -> "把有限操作记录转换成精确有序的 symminformation Association，生成两个格点，找出参考格点稳定子，并附加唯一的纯内部 C_infinity 操作。这些都只是准备步骤；此时尚未求解任何 Hamiltonian。",
  "The following immutable snapshot contains the ordered eight finite representatives and exact double-valued matrices generated from SpinLayerCorepresentations. It is bundled here so the tutorial has no optional runtime dependency. At the 000 seed the finite spatial action generates {0,0,0} and {1/2,1/2,0}; the site stabilizer is derived from the same records." -> "下面的不可变快照包含由 SpinLayerCorepresentations 生成的八个有序有限代表元和精确双值矩阵。它随教程一同提供，因此教程没有可选运行时依赖。在 000 种子处，有限空间作用生成 {0,0,0} 与 {1/2,1/2,0}；格点稳定子也由同一组记录导出。",
  "Replacing pz by the ordered spinor pair pzup,pzdn gives four bands. The full gray group supplies the antiunitary operation, while the symbolic Hamiltonian retains every symmetry-allowed spin-dependent term. This section reproduces the spinful graphene, bandManipulate, static-band, and hop workflows of the historical notebook." -> "把 pz 替换为有序旋量对 pzup、pzdn 后得到四条能带。完整灰群提供反幺正操作，而符号 Hamiltonian 保留所有对称性允许的自旋相关项。本节复现历史 notebook 中含自旋石墨烯、bandManipulate、静态能带和 hop 工作流。",
  "Spinful graphene, interactive SOC bands, and Wannier90 export" -> "含自旋石墨烯、交互式 SOC 能带与 Wannier90 导出",
  "Spinful MoS2 and two explicit basis orders" -> "含自旋 MoS2 与两种显式基顺序",
  "The historical notebook evaluated the spinful onsite model in two local basis orders. Both workflows are retained here. The first groups all spin-up orbitals before all spin-down orbitals; the subgroup is installed explicitly before symham is called." -> "历史 notebook 使用两种局域基顺序计算含自旋 onsite 模型。这里保留两套工作流。第一种先排列全部自旋向上轨道，再排列全部自旋向下轨道；调用 symham 前会显式安装子群。",
  "The numerical choices in this Weyl section and in the following nodal-line and C4 sections are illustrative symmetry-compatible parameters migrated from the historical example notebook; they are not presented as a fit to a named material." -> "本 Weyl 节以及后续 nodal-line 和 C4 各节中的数值，都是从历史示例 notebook 迁移来的、与对称性相容的示意参数，并不表示对某种具体材料的拟合。",
  StringJoin[
    "The full two-band c3Hamiltonian is used for the band plot and Chern-number calculation ",
    "below. To see its matrix elements more clearly, first display it on Gamma-A by setting ",
    "kx=ky=0:"
  ] -> StringJoin[
    "下面的能带和 Chern 数都使用完整的两带 c3Hamiltonian。为了看清矩阵元，先令 ",
    "kx=ky=0，显示 Gamma-A 线上的形式："
  ],
  "The original notebook also inspected the lowest-band constant-energy contours across several Brillouin zones. The polygon drawn in thick lines is the reciprocal primitive cell computed by init, not a hand-entered shape:" -> "原 notebook 还查看了跨越多个 Brillouin 区的最低能带等能轮廓。粗线绘制的多边形是 init 计算出的倒空间原胞，而不是手工输入的形状：",
  "The removed symmetryset option used to solve selected group operations inside symham. In 2.0 the symmetry is changed explicitly at initialization: close the desired subgroup, inspect the complete returned init rules, install them, and only then solve the new Hamiltonian." -> "已删除的 symmetryset 选项过去用于在 symham 内只求解选定群操作。2.0 中必须在初始化时显式改变对称性：闭合所需子群，查看返回的完整 init 规则，安装这些规则，然后再求解新的 Hamiltonian。",
  StringJoin[
    "Alternatively, place the spin-up and spin-down states next to each other for each ",
    "orbital. The order shown by orbitalTable is also the row and column order of the ",
    "Hamiltonian:"
  ] -> StringJoin[
    "也可以把每个轨道的自旋向上和自旋向下分量排在一起。orbitalTable 中的顺序就是 ",
    "Hamiltonian 的行列顺序："
  ],
  "The solved-shell overload of hop exports the exact cached real-space terms, so no Fourier coefficient has to be reconstructed from the displayed Bloch phases. The generated file is then inspected directly." -> "hop 的已求解键层形式直接导出精确缓存的实空间项，因此不必从显示的 Bloch 相位重建 Fourier 系数。随后直接查看生成的文件。",
  "A Gamma-point commutator is only a little-group check at k=0. It cannot test the momentum transformation, nonsymmorphic reciprocal sewing, or a non-Gamma little group. The release regression therefore also uses six generic momenta and reciprocal-sewn non-Gamma points for the seven-crystal-system models." -> "Gamma 点对易子只是 k=0 处的小群检查，无法检验动量变换、非点式倒格 sewing 或非 Gamma 小群。因此发布回归还会对七晶系模型使用六个一般动量点和经过倒格 sewing 的非 Gamma 点。",
  "A spinful antiunitary model" -> "含自旋的反幺正模型",
  "Calling a model-dependent function after that failed init exposes the missing session immediately." -> "在该 init 失败后调用依赖模型的函数，会立即暴露会话缺失。",
  "CsCl has two inequivalent Wyckoff orbits. The corner site carries one s orbital, while the body-centred site carries px, py, pz. The outer orders of wyckoffposition and basisFunctions match, so the inter-orbit block is physically a 1 by 3 rectangular hopping matrix." -> "CsCl 有两条不等价的 Wyckoff 轨道。角点格点携带一个 s 轨道，体心格点携带 px、py、pz。wyckoffposition 与 basisFunctions 的外层顺序一致，因此轨道间分块在物理上是一个 1×3 矩形 hopping 矩阵。",
  "Each magnetic-group record is {label,R,t,\"F\"|\"T\"}. The last field records whether the operation is unitary or antiunitary." -> "每条磁群记录都是 {label,R,t,\"F\"|\"T\"}；最后一个字段记录该操作是幺正还是反幺正。",
  "Evaluation boundaries and common errors" -> "计算边界与常见错误",
  "From graphene to unequal local dimensions" -> "从石墨烯到不等的局域维数",
  "initfromrep requires exact unitary matrices. The following upper-triangular matrix is not unitary and is rejected instead of being orthogonalized automatically." -> "initfromrep 要求精确幺正矩阵。下面的上三角矩阵不是幺正矩阵，因此会被拒绝，而不会被自动正交化。",
  "Numerical assignments for every symbolic lattice constant. MagneticTB does not impose an angstrom unit; all lengths must use one consistent unit." -> "为每个符号晶格常数给出数值。MagneticTB 不强制使用埃作为单位；所有长度必须采用一致单位。",
  "One nonempty local basis list per Wyckoff orbit. Different entries may have different dimensions, which produces rectangular inter-orbit hopping blocks." -> "每条 Wyckoff 轨道对应一个非空局域基列表。不同条目可以具有不同维数，从而产生矩形的轨道间 hopping 分块。",
  "Operation 25 is the antiunitary time-reversal operation in this ordered gray group. The output is its matrix only; the antiunitary flag remains in the operation record and tells the covariance check to complex-conjugate H(k)." -> "操作 25 是该有序灰群中的反幺正时间反演操作。输出只显示它的矩阵；反幺正标志仍保留在操作记录中，并指示协变检查对 H(k) 取复共轭。",
  "Restore the three-shell graphene session, then request an unprepared shell. The call stops with $Failed; it does not search another shell or change the session." -> "先恢复准备了三个键层的石墨烯会话，再请求尚未准备的键层。调用以 $Failed 停止；它不会搜索其他键层，也不会改变会话。",
  "Set pacletArchive to a built MagneticTB 2.x archive. PacletInstall installs a newer version alongside the normal Wolfram paclet registry and PacletFind reports the version that will be found. PacletInstall::samevers means that exactly the same version is already installed; either increase the archive version or explicitly uninstall that installed version before reinstalling. Quit the kernel before loading a newly installed build." -> "把 pacletArchive 设为已经构建的 MagneticTB 2.x 归档文件。PacletInstall 通过标准 Wolfram paclet 注册表安装新版本，PacletFind 会报告当前可找到的版本。PacletInstall::samevers 表示完全相同的版本已经安装；应提高归档版本号，或先显式卸载该已安装版本再重新安装。加载新安装的构建前请退出内核。",
  StringJoin[
    "A spin-space group allows the spin rotation to be specified separately from the spatial ",
    "operation. The following example has a pure spin C4 rotation and a half translation ",
    "combined with time reversal. In each dictionary entry, spin -> {S,0} denotes a unitary ",
    "operation and spin -> {S,1} an antiunitary one; S is the spin-rotation matrix."
  ] -> StringJoin[
    "自旋空间群可以分别指定空间操作和自旋旋转。下面的例子包含纯自旋 C4 ",
    "旋转，以及半平移与时间反演的组合。在字典中，spin -> {S,0} 表示幺正操作，spin -> ",
    "{S,1} 表示反幺正操作，其中 S 是自旋旋转矩阵。"
  ],
  "The following report evaluates every prepared operation at the generic momentum {Pi/7,Pi/11,Pi/13}. For a unitary operation it compares U H(k) U-dagger with H(p^(-T) k); for an antiunitary operation it compares U Conjugate[H(k)] U-dagger with H(-p^(-T) k). The per-operation residuals are computed from every independent Hamiltonian parameter basis, not from one fitted parameter choice." -> "下面的报告在一般动量 {Pi/7,Pi/11,Pi/13} 处计算每个已准备操作。对于幺正操作，它比较 U H(k) U-dagger 与 H(p^(-T) k)；对于反幺正操作，它比较 U Conjugate[H(k)] U-dagger 与 H(-p^(-T) k)。逐操作残差由每个独立 Hamiltonian 参数基计算，而不是只检查一组拟合参数。",
  "The outer basisFunctions list must have one entry per Wyckoff seed. This input has one seed but two basis lists, so initialization fails and clears the preceding graphene session." -> "basisFunctions 外层列表必须为每个 Wyckoff 种子提供一个条目。该输入只有一个种子却给出两个基列表，因此初始化失败并清除此前的石墨烯会话。",
  "The outer list enumerates inequivalent Wyckoff orbits. Each entry is {seed position, magnetic moment}; the same outer order is used by basisFunctions." -> "外层列表枚举不等价 Wyckoff 轨道。每个条目为 {种子位置, 磁矩}；basisFunctions 使用相同的外层顺序。",
  "Three primitive vectors stored as rows. Coordinates in wyckoffposition are fractional coordinates in this basis." -> "按行存储三条原胞基矢。wyckoffposition 中的坐标是相对于该基的分数坐标。",
  "Plot the archived bands" -> "绘制归档版本的能带",
  "MoS2: an eleven-band d-p orbital model" -> "MoS2：十一带 d-p 轨道模型",
  "SimpleCubicS: one s orbital on a simple-cubic site" -> "SimpleCubicS：简单立方格点上的单个 s 轨道",
  "The first Wyckoff orbit carries the five metal d orbitals. The second seed expands to two chalcogen sites, each with px, py, pz, giving eleven orbitals in total. The parameters below are illustrative symmetry-allowed amplitudes; this example demonstrates the complete orbital hierarchy and does not claim a fitted reproduction of a particular material data set." -> "第一条 Wyckoff 轨道携带五个金属 d 轨道。第二个种子展开成两个硫族元素格点，每个格点带有 px、py、pz，总计十一条轨道。下面的参数是对称性允许振幅的示意值；该例展示完整轨道层次，并不声称拟合复现某个特定材料数据集。",
  "This is the smallest nontrivial three-dimensional workflow. The gray Pm-3m operations keep the origin fixed, while shell 2 supplies the six nearest-neighbour translations. The result is shown as an actual one-band Bloch Hamiltonian rather than only a parameter count." -> "这是最小的非平凡三维工作流。灰 Pm-3m 操作保持原点不动，而第 2 键层提供六个最近邻平移。结果直接显示为实际的单带 Bloch Hamiltonian，而不只是参数数目。",
  "Wyckoff Database Models" -> "Wyckoff 数据库模型",
  "Magnetic Wyckoff conventions follow the BNS magnetic-group tables described by S. V. Gallego et al., \"Magnetic symmetry in the Bilbao Crystallographic Server,\" Journal of Applied Crystallography 45, 1236-1247 (2012), https://doi.org/10.1107/S0021889812042185. The bundled database is the executable source used here; the citation documents the corresponding magnetic Wyckoff data convention." -> "磁 Wyckoff 约定遵循 S. V. Gallego 等人在 \"Magnetic symmetry in the Bilbao Crystallographic Server,\" Journal of Applied Crystallography 45, 1236-1247 (2012)，https://doi.org/10.1107/S0021889812042185 中说明的 BNS 磁群表。随包数据库是这里实际执行的数据源；该文献说明了相应的磁 Wyckoff 数据约定。",
  "The magnetic group needs an antiunitary coset to complete this two-site orbit. Both branches below read the formal database seed and then initialize the complete physical model explicitly." -> "该磁群需要一个反幺正陪集来补全这条两格点轨道。下面两个分支都读取正式数据库种子，然后显式初始化完整物理模型。",
  "These ten examples cover all seven crystal systems and the exact BNS/Wyckoff entries used by the DirectProduct/Induced regression. Every branch visibly reads its seed from the paclet's formal wyckoffMSG.mx database and then gives its complete init, symham, representation-matrix, and bandManipulate workflow." -> "这十个范例覆盖七大晶系以及 DirectProduct/Induced 回归使用的精确 BNS/Wyckoff 条目。每个分支都明确从 paclet 正式 wyckoffMSG.mx 数据库读取种子，然后给出完整的 init、symham、表示矩阵与 bandManipulate 工作流。",
  " — overlay tight-binding bands with supplied reference eigenvalues." -> " — 将紧束缚能带与输入的参考本征值叠加比较。",
  "All displays the complete overlay; {emin,emax} restricts the visible energy interval." -> "All 显示完整叠加图；{emin,emax} 仅限制可见能量区间。",
  "A reference list whose length is not a multiple of npoint cannot be partitioned into path segments and is rejected by the underlying list operation:" -> "长度不是 npoint 整数倍的参考列表无法按路径段分组，底层列表操作会拒绝该输入：",
  "Construct a complete one-band tight-binding comparison and overlay four supplied reference samples:" -> "构造一个完整的单带紧束缚对比，并叠加四个给定的参考采样点：",
  "overlays tight-binding eigenvalues of H with supplied reference-band eigenvalues along the same path." -> "沿同一路径叠加 H 的紧束缚本征值与输入的参考能带本征值。",
  "Reference values normally come from a DFT or experimental workflow. Keep their path order and npoint grouping identical to comparisonPath before overlaying them:" -> "参考值通常来自 DFT 或实验工作流。叠加前必须保持其路径顺序和 npoint 分组与 comparisonPath 完全一致：",
  "Restrict only the visible vertical range of the same tight-binding/reference overlay:" -> "仅限制同一紧束缚/参考叠加图的可见纵轴范围：",
  "The plotRange option defaults to All and changes only the displayed vertical energy interval." -> "plotRange 选项默认为 All，只改变显示的纵向能量区间。",
  "The referenceBands argument is a list of {momentum,eigenvalue-list} records, grouped into npoint records per path segment." -> "referenceBands 是 {momentum,eigenvalue-list} 记录列表，每个路径段按 npoint 条记录分组。",
  "This example deliberately generates the reference values from the same analytic one-band model so the alignment can be inspected without an external DFT file." -> "本例有意从同一解析单带模型生成参考值，以便不依赖外部 DFT 文件就能检查对齐。",
  "None uses the orbital centers prepared by init. An explicit list supplies one fractional center per Hamiltonian orbital in the same order as the rows and columns of H." -> "None 使用 init 准备的轨道中心；显式列表必须按照 H 的行列顺序，为每个 Hamiltonian 轨道提供一个分数坐标中心。"
  ,"Triclinic [LongDash] BNS 1.3-a" -> "三斜晶系 [LongDash] BNS 1.3-a"
  ,"Triclinic [LongDash] BNS 2.6-i" -> "三斜晶系 [LongDash] BNS 2.6-i"
  ,"Monoclinic [LongDash] BNS 3.3-e" -> "单斜晶系 [LongDash] BNS 3.3-e"
  ,"Monoclinic [LongDash] BNS 3.4-a" -> "单斜晶系 [LongDash] BNS 3.4-a"
  ,"Orthorhombic [LongDash] BNS 16.3-q" -> "正交晶系 [LongDash] BNS 16.3-q"
  ,"Orthorhombic [LongDash] BNS 16.3-r" -> "正交晶系 [LongDash] BNS 16.3-r"
  ,"Tetragonal [LongDash] BNS 75.3-c" -> "四方晶系 [LongDash] BNS 75.3-c"
  ,"Trigonal [LongDash] BNS 143.3-a" -> "三方晶系 [LongDash] BNS 143.3-a"
  ,"Hexagonal [LongDash] BNS 168.111-b" -> "六方晶系 [LongDash] BNS 168.111-b"
  ,"Cubic [LongDash] BNS 195.3-a" -> "立方晶系 [LongDash] BNS 195.3-a"
  ,"Advanced programmatic interfaces" -> "供程序调用的高级接口"
  ,"MagneticTB constructs analytical tight-binding Hamiltonians from crystal symmetry and local orbital information. It supports ordinary space groups, magnetic space groups, and finite discrete spin-space groups." -> "MagneticTB 根据晶体对称性和局域轨道构造解析的紧束缚 Hamiltonian，可处理普通空间群、磁空间群以及有限离散自旋空间群。"
  ," — check program-generated model input without making it the current model." -> " — 检查其他程序生成的模型输入，但不把它设为当前模型。"
  ," — construct a Hamiltonian with the original implementation." -> " — 用原来的实现构造 Hamiltonian。"
  ," — explore bands obtained from the original implementation." -> " — 交互查看原实现得到的能带。"
  ," — find a generating subset of a complete finite group (2.x compatibility function)." -> " — 从完整有限群中寻找一组生成元（2.x 兼容函数）。"
  ," — initialize a model with the original implementation." -> " — 用原来的实现初始化模型。"
  ," — inspect the present model, orbital order, symmetry representation, and prepared bond shells." -> " — 查看当前模型的轨道顺序、对称表示和已经准备好的键层。"

  ,"writes the numerical band energies of H along path to a data file." -> "把 H 沿 path 的数值能量写入数据文件。"
  ,"After applying rules, the function samples the same reciprocal-fractional path convention used by bandplot. The output contains one row for each band and one column for each sampled momentum." -> "代入 rules 后，函数按 bandplot 相同的倒空间分数坐标路径取样。输出中每一行对应一条能带，每一列对应一个取样动量。"
  ,"Use FileNameJoin rather than a machine-specific literal path in notebooks that will be shared." -> "需要共享 notebook 时，请用 FileNameJoin 构造路径，不要写死某台机器上的绝对路径。"

  ,"draws the bands of an old-backend Hamiltonian and provides sliders for its parameters." -> "绘制原实现 Hamiltonian 的能带，并为其中的参数提供滑块。"
  ,"This function belongs to MagneticTBOld` and uses kxold, kyold, and kzold. Construct H with symhamold in the same fresh kernel before plotting it." -> "这个函数属于 MagneticTBOld`，动量变量为 kxold、kyold 和 kzold。请在同一个只加载旧版的新内核中先用 symhamold 构造 H，再画能带。"
  ,"Define the graphene path and open the interactive band panel of the original implementation:" -> "定义石墨烯的动量路径，并打开原实现的交互式能带面板："
  ,"Fix the original parameters to obtain a static band plot:" -> "固定旧版参数，得到一幅静态能带图："

  ,"draws the bands of H along path and provides sliders for the tight-binding parameters." -> "沿 path 绘制 H 的能带，并为紧束缚参数提供滑块。"
  ,"path is a list of reciprocal-fractional line segments and endpoint labels. npoint controls the sampling of each segment." -> "path 由若干倒空间分数坐标线段及其端点标签组成，npoint 决定每一段的取样密度。"
  ,"Every symbol in H other than kx, ky, and kz becomes an interactive parameter. The sliders start at zero, so the initial plot may show coincident flat bands when all Hamiltonian parameters vanish." -> "H 中除 kx、ky、kz 以外的符号都会成为可调参数。滑块的初值为零，因此刚打开时所有 Hamiltonian 参数都为零，若干条能带可能重合成一条水平线。"
  ,"No energy shift is added. For a genuinely non-Hermitian Hamiltonian, choose explicitly whether to plot the real part, imaginary part, or complex spectrum before calling this function." -> "函数不会自动平移能量。若 Hamiltonian 确实是非厄米的，应先明确要画本征值的实部、虚部还是完整复能谱。"
  ,StringJoin[
    "Use standardKPath[] as the first argument when MagneticTB should choose the built-in conventional ",
    "path from the current initialized lattice. bandManipulate receives the explicit path returned by ",
    "standardKPath[]."
  ] -> "若希望 MagneticTB 根据当前已初始化的晶格自动选择内置约定路径，可把 standardKPath[] 作为第一个参数。bandManipulate 使用 standardKPath[] 返回的显式路径。"

  ,"plots the numerical bands of H along a reciprocal-fractional momentum path after the model parameters are fixed." -> "固定模型参数后，沿倒空间分数坐标路径绘制 H 的数值能带。"
  ,"Each path segment is given by two reciprocal-fractional endpoints and their labels. MagneticTB converts these coordinates to the phase convention of H." -> "每段路径由两个倒空间分数坐标端点及其标签给出。MagneticTB 会把这些坐标换算到 H 所采用的相位约定。"
  ,"rules must assign numerical values to every tight-binding parameter. plotRange changes only the displayed energy window." -> "rules 必须为全部紧束缚参数指定数值。plotRange 只改变图中显示的能量范围。"

  ,"builds a new set of init rules for the subgroup generated by selected symmetry operations." -> "根据选定的对称操作所生成的子群，构造一套新的 init 输入规则。"
  ,"Use this function when studying a deliberate symmetry-breaking pattern. The selected operations generate the retained subgroup; all atoms of the present model are kept." -> "研究人为指定的对称性破缺时可使用这个函数。所选操作生成要保留的子群，当前模型中的原子全部保留。"
  ,"If the subgroup no longer relates all atoms of a Wyckoff orbit, that orbit is split into separate entries of the new wyckoffposition input." -> "若保留的子群不能再把某个 Wyckoff 轨道中的全部原子互相联系，该轨道会拆成新的 wyckoffposition 条目。"
  ,"The current model is not changed. Inspect the returned rules and evaluate init@@rules explicitly. Use {} to retain only the identity." -> "函数不会改动当前模型。请先检查返回的规则，再显式执行 init@@rules。若只保留单位元，输入 {}。"

  ,"plots tight-binding bands from H together with VASP eigenvalues read from EIGENVAL by vaspEig." -> "把 H 的紧束缚能带与 vaspEig 从 EIGENVAL 读取的 VASP 本征值画在同一张图中。"
  ,"Use vaspEig[file,efermi,spin,startBand,endBand] to read the required bands from a VASP EIGENVAL file before calling compareBand." -> "调用 compareBand 前，先用 vaspEig[file,efermi,spin,startBand,endBand] 从 VASP EIGENVAL 文件读取需要比较的能带。"
  ,"The EIGENVAL k points must follow the same path-segment order as path. npoint is the number of VASP samples in each segment." -> "EIGENVAL 中的 k 点必须与 path 采用相同的路径段顺序；npoint 是每一段中的 VASP 取样点数。"
  ,"plotRange changes only the visible energy interval." -> "plotRange 只改变图中显示的能量区间。"
  ,"Read bands 21-23 from the included VASP EIGENVAL file, construct a three-orbital hexagonal model, and compare its tight-binding bands with the VASP bands:" -> "从随页提供的 VASP EIGENVAL 文件读取第 21-23 条能带，构造一个六方晶格三轨道模型，再把紧束缚能带与 VASP 能带画在一起比较："
  ,"Using the same VASP and tight-binding bands, plotRange -> {-1,2} restricts only the visible energy window:" -> "继续使用上面的 VASP 与紧束缚能带；plotRange -> {-1,2} 只限制图中显示的能量范围："

  ,"checks and prepares program-generated MagneticTB input without making it the current model." -> "检查并整理其他程序生成的 MagneticTB 输入，但不把它设为当前模型。"
  ,"This is an advanced interface for another program that generates MagneticTB models. Ordinary notebook users should call init or initfromrep instead." -> "这是供其他建模程序调用的高级接口。一般在 notebook 中建模时，请直接使用 init 或 initfromrep。"
  ,"input contains the same lattice, Wyckoff, symmetry, and orbital information used by init, written with the documented string keys. The result contains the checked sites, representation matrices, and bond geometry, but no hopping shell is solved and the current model is unchanged." -> "input 包含与 init 相同的晶格、Wyckoff 位置、对称性和轨道信息，只是采用文档规定的字符串键。返回结果给出检查后的格点、表示矩阵和键几何；它不会求解任何 hopping，也不会改变当前模型。"
  ,"A model-generating program can inspect the exact representation returned here before deciding whether to install the model through the public initialization interface:" -> "生成模型的程序可以先检查这里返回的精确表示，再决定是否通过公开初始化接口把它设为当前模型："
  ,"All physical inputs must be present. For example, a lattice by itself does not define a model:" -> "物理输入必须完整。例如，只有晶格还不足以定义一个模型："

  ,"returns the read-only data of the model most recently prepared by init or initfromrep." -> "返回最近一次由 init 或 initfromrep 建立的模型数据，供只读查看。"
  ,"The result contains the physical model, full symmetry representation, prepared bond shells, and any shell solutions already obtained by symham. It is intended for inspection and advanced programmatic use." -> "返回结果包括物理模型、完整对称表示、已经准备好的键层，以及 symham 已求出的各层结果，主要用于检查模型和高级程序调用。"
  ,"Calling CurrentModelSession before a successful initialization returns $Failed." -> "在成功初始化之前调用 CurrentModelSession 会返回 $Failed。"
  ,"Reading this value does not change the model, initialize another model, or solve a missing shell." -> "读取这些数据不会改变模型，不会重新初始化，也不会求解尚未计算的键层。"
  ,"Initialize a complete model before inspecting it:" -> "先初始化一个完整模型，再查看其中的数据："
  ,"Read the representation mode, generated site orbits, local bases, and prepared shell numbers of the graphene model. This inspection performs no additional Hamiltonian calculation:" -> "查看石墨烯模型所采用的表示模式、生成的格点轨道、局域基底和已经准备好的键层编号。这个操作不会增加任何 Hamiltonian 计算："

  ,"generates a complete finite group from its generators and an explicitly supplied multiplication law." -> "按显式给定的乘法规则，由生成元生成完整有限群。"
  ,"The returned list begins with the identity. For the trivial group, the result therefore contains the identity rather than an empty list." -> "返回列表以单位元开头。因此，平凡群返回 {E}，而不是空列表。"
  ,"The caller supplies both multiplication and equality because magnetic-space-group and spin-space-group elements need not use the same concrete product." -> "调用者需要同时给出乘法和相等判据，因为磁空间群元与自旋空间群元的具体乘法形式并不相同。"
  ,"finds a small generating subset of a complete finite group." -> "从完整有限群中寻找一组较小的生成元。"
  ,"getGenerator is retained for notebooks that already call it. Supply the complete group and its multiplication law explicitly." -> "为兼容已有 notebook，getGenerator 仍可继续使用。调用时必须显式给出完整群及其乘法规则。"
  ,"The complete group, identity, and multiplication law must be supplied. The function does not infer a product from the appearance of the elements." -> "必须给出完整群、单位元和乘法规则。函数不会仅凭群元的写法猜测它们如何相乘。"

  ,"writes the real-space hoppings obtained from solved shells 1 through n in Wannier90 HR format." -> "把已经求出的第 1 到第 n 键层的实空间 hopping 写成 Wannier90 HR 格式。"
  ,"writes only the selected solved shells in Wannier90 HR format." -> "只把列表中选定且已经求出的键层写成 Wannier90 HR 格式。"
  ,"Fourier-decomposes a supplied finite exponential Bloch Hamiltonian and writes the resulting real-space matrices in Wannier90 HR format." -> "对输入的有限指数形式 Bloch Hamiltonian 作 Fourier 分解，并把得到的实空间矩阵写成 Wannier90 HR 格式。"
  ,"For hop[n,rules] or hop[{n1,n2,...},rules], first evaluate symham for every shell that will be exported. The integer form selects shells 1 through n, whereas the list form selects exactly the listed shells. hop reads those exact real-space matrices and does not solve a missing shell." -> "使用 hop[n,rules] 或 hop[{n1,n2,...},rules] 前，应先对要导出的每个键层计算 symham。整数形式选择第 1 到第 n 层，列表形式只选择明确列出的键层。hop 直接读取这些精确实空间矩阵，不会替你求解缺少的键层。"
  ,"For hop[H,rules], H must be written as a finite sum of exponential Bloch phases so that integer lattice translations can be identified." -> "使用 hop[H,rules] 时，H 必须写成有限个 Bloch 指数相位之和，这样才能识别对应的整数晶格平移。"
  ,"Set \"hrExport\" to a file or directory path to write wannier90_hr.dat. With None, the generated text is returned in the notebook. Use \"wcc\" -> Automatic when the orbital centers should be taken from init." -> "把 \"hrExport\" 设为文件或目录路径即可写出 wannier90_hr.dat；设为 None 时，生成的文本直接返回 notebook。若轨道中心取自 init，请用 \"wcc\" -> Automatic。"
  ,"Set \"hrExport\" to a path to write wannier90_hr.dat; with None, the text is shown in the notebook. Use \"wcc\" -> Automatic when the orbital centers should come from init." -> "把 \"hrExport\" 设为路径即可写出 wannier90_hr.dat；设为 None 时，文本直接显示在 notebook 中。若轨道中心取自 init，请使用 \"wcc\" -> Automatic。"
  ,"For shell export, select the solution previously obtained with the same Hermitian setting. This option is not used for a supplied matrix H(k)." -> "按键层导出时，选择先前用相同 Hermitian 设置求出的结果。若直接输入矩阵 H(k)，这个选项不起作用。"
  ,"For shell export, select the solution previously obtained with the same \"Iterative\" or \"Stacked\" method." -> "按键层导出时，选择先前用同一 \"Iterative\" 或 \"Stacked\" 方法求出的结果。"
  ,"For shell export, select the solution previously obtained with the same \"None\", \"Basic\", or \"Full\" validation level." -> "按键层导出时，选择先前用相同 \"None\"、\"Basic\" 或 \"Full\" 检查等级求出的结果。"
  ,"Initialize graphene, solve three shells, write a Wannier90 HR file, and inspect its header and first hopping entries:" -> "初始化石墨烯，求出三个键层，写入 Wannier90 HR 文件，并查看文件头和前几条 hopping："
  ,"Initialize graphene and construct its Hamiltonian from the first three neighbour contributions:" -> "初始化石墨烯，并由前三个近邻项构造 Hamiltonian："
  ,"Convert the Hamiltonian directly to a Wannier90 HR file:" -> "可以直接把 Hamiltonian 转换成 Wannier90 HR 文件："
  ,"Alternatively, export only a selected neighbour contribution:" -> "也可以只把指定的近邻项转换成 Wannier90 HR 文件："
  ,"This export uses the Hermitian shell obtained with the iterative method and Basic checks. It takes the orbital centers from the present model and writes ten decimal digits:" -> "下面导出用迭代方法和 Basic 检查得到的厄米键层，从当前模型读取轨道中心，并保留十位小数："

  ,"sets up a tight-binding model when the exact symmetry matrices are supplied directly instead of being derived from basis functions." -> "直接输入精确对称表示矩阵来建立紧束缚模型，而不是由基函数推导这些矩阵。"
  ,"Use this interface when another calculation has already produced the symmetry matrices. It is deliberately separate from init: basisFunctions, SiteLocalData, and GenerateSymmetryGroup are not accepted." -> "如果其他计算已经给出了对称表示矩阵，就使用这个接口。它与 init 分开，不能输入 basisFunctions、SiteLocalData 或 GenerateSymmetryGroup。"
  ,"symminformation must contain the complete finite group in a fixed order, and repinformation must contain one exact unitary matrix for every operation in exactly that order." -> "symminformation 必须按固定顺序给出完整有限群；repinformation 必须以完全相同的顺序，为每个操作给出一个精确幺正矩阵。"
  ,"For the supported continuous case, one pure-internal C_infinity generator may be supplied together with the finite matrices through the documented Discrete and Continuous form. The continuous operation constrains the local states but does not move atoms or generate bond orbits." -> "在目前支持的连续情形中，可以用文档规定的 Discrete/Continuous 形式，在有限矩阵之外再给出一个纯内部 C_infinity 生成元。这个连续操作只约束局域态，不移动原子，也不生成新的键轨道。"
  ,"By default initfromrep finds the first 10 complete bond shells. If a model needs a farther neighbour hopping, choose a larger InitialBondShells value and run initfromrep again; symham does not search for new bonds." -> "initfromrep 默认准备前 10 个完整近邻 hopping 层。若模型需要更远的 hopping，请设置更大的 InitialBondShells 并重新运行 initfromrep；symham 不会自行搜索新的近邻。"
  ,"Every call to initfromrep starts a new model. If the ordered group and representation data are inconsistent, initialization fails and MagneticTB does not continue with the preceding model." -> "每次运行 initfromrep 都是在建立新模型。若有序群与表示数据不一致，初始化会失败，MagneticTB 不会沿用上一个模型。"
  ,"Names used to display the coordinates of the supplied representation. Automatic uses orb1, orb2, and so on; these labels do not change the matrices." -> "为输入表示的各个坐标指定显示名称。Automatic 使用 orb1、orb2 等名称；这些标签不会改变矩阵。"
  ,"Number of complete bond shells prepared during initialization. Shell 1 contains onsite terms; shell 2 is the first nonzero-distance neighbour shell." -> "初始化时准备的完整近邻 hopping 层数。第 1 层是 onsite 项；第 2 层是第一个非零距离近邻层。"
  ,"Use \"DirectProduct\" when repinformation already gives the full local action, or \"Induced\" when it gives the representation of a reference site's symmetry group." -> "若 repinformation 已给出完整局域作用，使用 \"DirectProduct\"；若它给出的是参考格点位置稳定子群的表示，使用 \"Induced\"。"

  ,"sets up a model with the original MagneticTB implementation shipped as an isolated compatibility backend." -> "使用随包提供且完全隔离的原 MagneticTB 实现建立模型。"
  ,"Use initold for notebooks written for the original implementation, or whenever the original and current constructions need to be compared." -> "原版 notebook 可以继续使用 initold；需要比较原实现与当前实现时也可使用它。"
  ,"Its physical inputs are the original lattice, Wyckoff, symmetry, and basis-function options, renamed with the suffix old so that they cannot be confused with the new interface." -> "物理输入仍是原来的晶格、Wyckoff 位置、对称性和基函数选项，只是在名称后加 old，以免与新版接口混淆。"
  ,"MagneticTBOld` and MagneticTB` do not load or call one another. Evaluate them in separate fresh Wolfram kernels when comparing results." -> "MagneticTBOld` 与 MagneticTB` 不会互相加载或调用。比较结果时，请分别在两个全新的 Wolfram 内核中运行。"
  ,"Initialize graphene with the original implementation and display its onsite Hamiltonian. Run this example in a kernel that loads only MagneticTBOld`:" -> "用原来的实现初始化石墨烯并显示 onsite Hamiltonian。请在只加载 MagneticTBOld` 的内核中运行这个例子："
  ,"All six original options remain accepted together. This complete call displays the resulting onsite Hamiltonian in a kernel that loads only MagneticTBOld`:" -> "原来的六个选项仍可同时使用。下面给出完整调用，并在只加载 MagneticTBOld` 的内核中显示 onsite Hamiltonian："

  ,"sets up a tight-binding model from its lattice, Wyckoff positions, symmetry, and local orbitals." -> "根据晶格、Wyckoff 位置、对称性和局域轨道建立紧束缚模型。"
  ,"The five main inputs have a direct physical meaning: lattice and lattpar define the Bravais lattice, wyckoffposition specifies the inequivalent sites and their moments, symminformation gives the symmetry operations, and basisFunctions specifies the local orbital space on each Wyckoff orbit." -> "五项主要输入都对应明确的物理信息：lattice 与 lattpar 定义 Bravais 晶格；wyckoffposition 给出不等价格点及其磁矩；symminformation 给出对称操作；basisFunctions 给出每个 Wyckoff 轨道上的局域轨道空间。"
  ,"After a successful call, MagneticTB has determined the symmetry-related atoms, the common orbital order of all matrices, and the requested bond shells. The hopping parameters themselves are obtained later, one shell at a time, by symham." -> "初始化成功后，MagneticTB 已经确定对称性联系的原子、所有矩阵共同采用的轨道顺序，以及所需的键层。hopping 参数随后由 symham 按键层逐一求出。"
  ,"InitialBondShells is 10 by default. This is the exact range prepared by init; symham does not search a farther shell on its own." -> "InitialBondShells 默认为 10，这就是 init 准备好的完整范围；symham 不会自行搜索更远的键层。"
  ,"RepresentationMode selects how the full orbital representation is built. Use \"DirectProduct\" when the action on the complete local basis is known, and \"Induced\" when the input starts from a representation of a site-symmetry group." -> "RepresentationMode 决定如何构造完整轨道表示。已知对称性在完整局域基底上的作用时使用 \"DirectProduct\"；从位置稳定子群表示出发时使用 \"Induced\"。"
  ,"Calling init starts a new model. If the new input is inconsistent, the call stops and the preceding model is not kept as an implicit fallback." -> "每次调用 init 都是在建立一个新模型。若新输入不自洽，计算会停止，不会暗中继续使用上一个模型。"
  ,"Boolean compatibility option retained for existing notebooks. It does not alter the Hamiltonian or select another algorithm." -> "为兼容已有 notebook 而保留的布尔选项。它不改变 Hamiltonian，也不选择其他算法。"
  ,"Number of complete bond shells prepared by init. Shell 1 contains the zero-distance onsite terms; shell 2 is the first nonzero-distance neighbour shell." -> "init 要准备的完整键层数。第 1 层是零距离 onsite 项，第 2 层是第一个非零距离近邻层。"
  ,"\"DirectProduct\" uses the local action on the full basis together with the permutation of equivalent sites. \"Induced\" starts from the representation on a reference site's symmetry group and transports it to the other sites." -> "\"DirectProduct\" 把完整局域基底上的作用与等价格点置换组合起来。\"Induced\" 从参考格点位置稳定子群的表示出发，把局域态搬运到其他等价点。"
  ,"\"DirectProduct\" uses the local action on the full basis together with the permutation of equivalent sites. \"Induced\" derives the reference-site representation from basisFunctions and transports it to the other sites. Use initfromrep when exact site-symmetry matrices are already known." -> "\"DirectProduct\" 把完整局域基上的作用与等价格点置换组合起来。\"Induced\" 根据 basisFunctions 推导参考点表示，再把它搬运到其他等价点。若已经知道精确的格点对称矩阵，请使用 initfromrep。"
  ,"Data used only in Induced mode. Automatic derives the reference-site representation from basisFunctions. Program-generated input may instead give the reference site, its site-symmetry operations and matrices, and optional coset representatives explicitly." -> "仅用于 Induced 模式。Automatic 从 basisFunctions 推导参考点表示；程序生成的输入也可以显式给出参考点、位置稳定子群操作及其矩阵，并可选择给出陪集代表元。"
  ,"Initialize graphene with the complete ordered list returned by msgop. GenerateSymmetryGroup -> False means that this list is already the whole group; it does not mean that every operation is a generator:" -> "用 msgop 返回的完整有序群元列表初始化石墨烯。GenerateSymmetryGroup -> False 表示这已经是完整群，并不表示每个操作都是生成元："
  ,"These four operations generate all 48 operations of the gray group used above:" -> "下面四个操作可以生成上述灰群的全部 48 个操作："
  ,"Inspect the Hamiltonian ordering. The two rows are the symmetry-related carbon sites generated from the single Wyckoff seed:" -> "检查 Hamiltonian 的轨道顺序。下面两行是由一个 Wyckoff 种子生成的两个对称等价碳原子："
  ,"InitialBondShells fixes the available shell range. This graphene input prepares only the onsite shell, whose Hamiltonian is displayed:" -> "InitialBondShells 固定可用键层的范围。下面的石墨烯输入只准备 onsite 层，并显示该层的 Hamiltonian："
  ,"GenerateSymmetryGroup -> True accepts generators instead of a complete operation list. Here a pure spin rotation and a half translation combined with time reversal generate an eight-element collinear spin-space group. The displayed matrices show their action on the local spinor:" -> "GenerateSymmetryGroup -> True 允许只输入生成元。这里一个纯自旋转动和一个与时间反演结合的半平移生成八元共线自旋空间群；所显示的矩阵给出它们在局域旋量上的作用："
  ,"An unknown option stops initialization before a model is created, and an earlier model is not used implicitly:" -> "未知选项会使初始化在模型建立之前停止，而且不会暗中沿用之前的模型："

  ,"reads the symmetry operations of a magnetic layer group from the bundled database." -> "从随包数据库读取磁层群的对称操作。"
  ,"reads the symmetry operations of a magnetic rod group from the bundled database." -> "从随包数据库读取磁杆群的对称操作。"
  ,"As for msgop, each item contains the operation label, point matrix, fractional translation, and unitary or antiunitary flag." -> "与 msgop 相同，每一项包含操作标签、点操作矩阵、分数平移以及幺正或反幺正标志。"
  ,"reads the symmetry operations of a magnetic space group from the bundled database." -> "从随包数据库读取磁空间群的对称操作。"
  ,"Use gray[n] for a gray group, or typeI[n], typeIII[n], and typeIV[n] for the corresponding magnetic-group types." -> "灰群使用 gray[n]；I、III、IV 型磁群分别使用 typeI[n]、typeIII[n] 和 typeIV[n]。"
  ,"Each returned item is {label,R,t,flag}: R and t are the point and translation parts in fractional coordinates, and flag states whether the operation is unitary or antiunitary." -> "每个群元写成 {label,R,t,flag}：R 和 t 分别是点操作与分数坐标平移，flag 表示它是幺正还是反幺正操作。"
  ,"The list may be passed directly to symminformation in init. Its order is also the order used when displaying the corresponding representation matrices." -> "这个列表可以直接作为 init 的 symminformation。列表顺序同时决定表示矩阵的显示顺序。"

  ,"lists the orbitals in the row and column order used by the Hamiltonian and symmetry matrices." -> "按 Hamiltonian 和对称矩阵实际采用的行列顺序列出全部轨道。"
  ,"Each row identifies the Wyckoff orbit, the equivalent atom, its fractional and Cartesian coordinates, the local basis state, and any spin label." -> "每一行给出 Wyckoff 轨道、等价原子、分数与笛卡尔坐标、局域基态以及自旋标签。"
  ,"OrbitalID is the matrix index. The same order is used by symham, showHamiltonianBasis, and every full representation matrix." -> "OrbitalID 就是矩阵指标。symham、showHamiltonianBasis 和所有完整表示矩阵都采用同一顺序。"
  ,"orbitalTable only displays the model prepared by init; it does not call symham." -> "orbitalTable 只显示 init 建立的模型，不会调用 symham。"

  ,"gives the exact matrix of each point operation on a chosen local orbital or spinor basis." -> "给出各点操作在所选局域轨道或旋量基底上的精确矩阵。"
  ,"Each operation is written as {label,R,t,flag}. Only its point action enters the local orbital matrix; the translation is retained so the same operation data can be used by init." -> "每个操作写成 {label,R,t,flag}。局域轨道矩阵只由点操作决定；平移仍保留在群元中，因此同一组操作数据可以直接交给 init。"
  ,"basis may contain scalar orbitals such as px and py or two-component spinors. The matrices are returned in the same operation order and act on the basis in the order supplied." -> "basis 可以包含 px、py 等标量轨道，也可以包含二分量旋量。返回矩阵与输入群元顺序一致，并按给定的基底顺序作用。"

  ,"reads a Wannier90 HR file and returns its lattice translations and real-space hopping matrices." -> "读取 Wannier90 HR 文件，返回晶格平移和实空间 hopping 矩阵。"
  ,"The complete header, degeneracies, translation vectors, orbital indices, and complex matrix elements are checked before a result is returned. A malformed or incomplete file is rejected." -> "返回结果前会检查完整文件头、简并度、平移向量、轨道指标和复矩阵元。格式错误或不完整的文件会被拒绝。"
  ,"The ncell option may then keep only translations inside a chosen real-space box." -> "随后可用 ncell 只保留指定实空间范围内的平移。"
  ,"Create a graphene HR file with hop, read it back, and display its lattice translations and real-space hopping matrices:" -> "先用 hop 写出石墨烯 HR 文件，再读回并显示其中的晶格平移和实空间 hopping 矩阵："
  ,"Keep only the onsite cell and set matrix elements smaller than 10^-12 to zero after reading the file:" -> "读取文件后只保留 onsite 晶胞，并把小于 10^-12 的矩阵元置零："

  ,"shows the bond length and all directed bonds in prepared shell n." -> "显示已经准备好的第 n 键层的键长和全部有向键。"
  ,"Shell 1 is the zero-distance onsite shell; shell 2 is the shortest nonzero-distance neighbour shell." -> "第 1 层是零距离 onsite 层；第 2 层是距离最短的非零近邻层。"
  ,"The output identifies the starting site, the translated ending site, and the displacement of each directed bond." -> "输出给出每条有向键的起点、平移后的终点以及位移向量。"
  ,"The shell must have been prepared by init. showbonds does not search farther and does not evaluate symham." -> "该键层必须已经由 init 准备。showbonds 不会向更远距离搜索，也不会计算 symham。"

  ,"shows the basis states in the exact row and column order of the Hamiltonian." -> "按 Hamiltonian 的精确行列顺序显示基态。"
  ,"identifies the final bra state and initial ket state of matrix element H[[i,j]]." -> "指出矩阵元 H[[i,j]] 的末态左矢和初态右矢。"
  ,"For H[[i,j]], column j is the initial ket and row i is the final bra. In Induced mode the table also shows how a reference-site state was transported to an equivalent site." -> "对 H[[i,j]]，第 j 列是初态右矢，第 i 行是末态左矢。在 Induced 模式下，表中还会显示参考点基态如何搬运到等价点。"

  ,"shows the representative real-space hopping associated with each independent parameter in a solved shell." -> "显示已求键层中每个独立参数所对应的代表性实空间 hopping。"
  ,"shows every real-space matrix element generated from parameter p by symmetry and Hermiticity." -> "显示参数 p 通过对称性和厄米关系生成的全部实空间矩阵元。"
  ,"Evaluate symham[n] first. The table then connects a symbol such as t1 to a definite initial orbital, final orbital, and lattice translation." -> "请先计算 symham[n]。随后表格会把 t1 这样的参数明确对应到初态轨道、末态轨道和晶格平移。"
  ,"The two-argument form follows the same parameter through all symmetry-related bonds, making it possible to see where every term of the Bloch Hamiltonian comes from." -> "双参数形式沿全部对称等价键追踪同一个参数，从而看清 Bloch Hamiltonian 中每一项的实空间来源。"
  ,"This is an inspection function; it does not solve a missing shell." -> "这是查看结果的函数，不会求解尚未计算的键层。"
  ,"Choose the earlier Hermitian or non-Hermitian shell calculation whose hopping origins are to be displayed." -> "选择先前算出的厄米或非厄米键层，并显示其中 hopping 的来源。"
  ,"Choose the earlier shell calculation made with the \"Iterative\" or \"Stacked\" null-space method. A missing calculation is not run automatically." -> "选择先前用 \"Iterative\" 或 \"Stacked\" 零空间方法得到的键层。若该结果不存在，函数不会自动补算。"
  ,"Choose the earlier shell calculation made with the corresponding \"None\", \"Basic\", or \"Full\" consistency check." -> "选择先前采用相应 \"None\"、\"Basic\" 或 \"Full\" 一致性检查得到的键层。"
  ,"Initialize graphene, solve its nearest-neighbour shell, and display which real-space hopping is represented by t1:" -> "初始化石墨烯并求出最近邻键层，再查看 t1 代表哪一条实空间 hopping："
  ,"After solving the non-Hermitian shell, display the independent forward and reverse hoppings instead of identifying them by Hermiticity:" -> "求出非厄米键层后，分别显示正向和反向的独立 hopping，不再用厄米关系把它们等同："
  ,"Use the same null-space method and checking level in both calls. The shell is solved explicitly before its hopping origins are displayed:" -> "两次调用应使用同一个零空间方法和检查等级。下面先显式求解键层，再显示 hopping 的来源："

  ,"shows the Wyckoff positions and symmetry-allowed magnetic moments of a magnetic space group." -> "显示磁空间群的 Wyckoff 位置以及对称性允许的磁矩方向。"
  ,"group may be a database identifier such as gray[n] or the corresponding magnetic-group number." -> "group 可以是 gray[n] 这样的数据库标识，也可以是相应的磁群编号。"
  ,"Each row gives a representative position. Symbols x, y, and z are free Wyckoff coordinates; the moment entry gives the directions allowed by the site symmetry." -> "每行给出一个代表位置。x、y、z 表示自由 Wyckoff 坐标；磁矩一栏给出点群允许的方向。"
  ,"Display the Wyckoff table of gray magnetic space group 193:" -> "显示第 193 号灰磁空间群的 Wyckoff 表："

  ,"shows how selected symmetry operations act on the complete Hamiltonian basis." -> "显示所选对称操作在完整 Hamiltonian 基底上的作用。"
  ,"Automatic displays a generating set; All displays every operation. An integer or list of integers selects operations by their prepared order." -> "Automatic 显示一组生成元，All 显示全部操作。整数或整数列表按准备好的群元顺序选择操作。"
  ,"The table places the spatial operation, spin action, and unitary or antiunitary type beside its matrix. Since antiunitarity is already part of the operation data, the matrix is shown without a separate formal K." -> "表格把空间操作、自旋作用以及幺正或反幺正类型与矩阵并列显示。反幺正性已经写在操作数据中，因此矩阵后不再另加形式上的 K。"
  ,"The function only reads the representation prepared by init or initfromrep." -> "函数只读取 init 或 initfromrep 已建立的表示。"

  ,"changes the Bloch basis convention of H by including the orbital-center phases." -> "把轨道中心相位吸收到 Bloch 基底中，从而改变 H 的基底约定。"
  ,"In convention I, the Bloch phase is attached to the full displacement between orbitals. Convention II moves the intracell orbital-center phases into the basis states." -> "在约定 I 中，Bloch 相位对应轨道之间的完整位移；约定 II 把晶胞内轨道中心的相位移入基态。"
  ,StringJoin[
    "Initialize the model first so that the fractional orbital centers are known, or supply them with ",
    "\"wcc\". An explicit list must follow the row and column order of H, and lets symhamII convert a ",
    "user-supplied Hamiltonian without relying on the current init session. The transformation changes ",
    "the matrix phases but not the energy spectrum."
  ] -> "应先初始化模型以确定分数坐标轨道中心，也可以用 \"wcc\" 显式给出。显式列表必须遵循 H 的行列顺序，使 symhamII 可以在不依赖当前 init 会话的情况下转换用户自己的 Hamiltonian。这个变换改变矩阵中的相位位置，但不改变能谱。"

  ,"returns bond shell n of the symmetry-allowed Hamiltonian obtained with the original MagneticTB algorithm." -> "返回原 MagneticTB 算法得到的第 n 键层对称性允许 Hamiltonian。"
  ,"Call initold first. The shell numbering and Bloch Hamiltonian follow the original package, while every symbol carries the old suffix where needed for isolation." -> "请先调用 initold。键层编号和 Bloch Hamiltonian 沿用原包约定；需要隔离的符号都带 old 后缀。"
  ,"symhamold is kept for existing notebooks and for independent old/new comparisons. It does not call the new linear-algebra implementation." -> "保留 symhamold 是为了运行已有 notebook 和独立比较新旧结果；它不会调用新版线性代数实现。"
  ,"Subset of symmetry-operation indices used by the original algorithm. This option exists only in MagneticTBOld`." -> "原算法使用的对称操作编号子集。这个选项只存在于 MagneticTBOld`。"
  ,"Add the onsite and nearest-neighbour terms to obtain the original Bloch Hamiltonian:" -> "把 onsite 和最近邻项相加，得到原实现的 Bloch Hamiltonian："
  ,"CartesianCoordinates -> True writes the nearest-neighbour term in Cartesian momentum coordinates. symmetrysetold is an option of the original implementation; here All uses every prepared operation:" -> "CartesianCoordinates -> True 用笛卡尔动量坐标写出最近邻项。symmetrysetold 是原实现的选项；这里用 All 选择全部已准备操作："

  ,"returns the symmetry-allowed contribution of bond shell n to the Bloch Hamiltonian." -> "返回第 n 键层对 Bloch Hamiltonian 的对称性允许贡献。"
  ,"constructs the same bond shell without imposing the Hermitian relation between a hopping and its reverse." -> "构造同一个键层，但不施加 hopping 与反向 hopping 之间的厄米关系。"
  ,"Call init first. For the selected shell, symham relates all real-space hopping matrices by the prepared symmetry operations and solves for the independent amplitudes. It then places the corresponding phase factors in the Bloch Hamiltonian." -> "请先调用 init。对选定键层，symham 用对称操作联系全部实空间 hopping 矩阵，求出独立振幅，再把对应相位因子写入 Bloch Hamiltonian。"
  ,"Shell 1 contains onsite terms. Shell 2 is the first nonzero-distance neighbour shell, shell 3 is the next one, and so forth. A complete model is obtained by adding the required shell contributions." -> "第 1 层包含 onsite 项，第 2 层是第一个非零距离近邻层，第 3 层是再下一层，依此类推。把物理模型需要的各层相加即可得到总 Hamiltonian。"
  ,"The default result is Hermitian and uses the iterative exact null-space solver with Basic validation. CartesianCoordinates -> False keeps the reciprocal-fractional momentum convention used by the input lattice." -> "默认结果是厄米的，采用迭代式精确零空间算法并进行 Basic 检查。CartesianCoordinates -> False 保留输入晶格所对应的倒空间分数坐标动量约定。"
  ,"symham does not repeat init or the bond search. If shell n was not prepared, rerun init with InitialBondShells -> n." -> "symham 不会重复 init 或重新搜索键。若第 n 层尚未准备，请用 InitialBondShells -> n 重新运行 init。"
  ,"True imposes H_ij(R)=Conjugate[H_ji(-R)] and gives a Hermitian shell contribution. False treats the two directions as independent amplitudes." -> "True 施加 H_ij(R)=Conjugate[H_ji(-R)]，得到厄米键层；False 把两个方向视为独立振幅。"
  ,"Exact method used to find the common null space of the symmetry constraints. \"Iterative\" applies the constraints successively; \"Stacked\" solves their vertically joined matrix. A failed method is reported rather than replaced silently." -> "求全部对称约束公共零空间的精确方法。\"Iterative\" 逐个加入约束；\"Stacked\" 求纵向拼接约束矩阵的零空间。若所选方法失败，程序会直接报告，不会暗中换用另一算法。"
  ,"Amount of consistency checking after the null space is found. \"None\" skips these checks, \"Basic\" checks the constraint residuals, and \"Full\" also checks independence and rank-nullity." -> "找到零空间后的一致性检查等级。\"None\" 不做这些检查；\"Basic\" 检查约束余量；\"Full\" 还检查基底独立性和秩—零度关系。"
  ,"Hermitian -> False keeps the two hopping directions independent, as shown by the non-Hermitian matrix in Basic Examples. KernelMethod and ValidationLevel change only how the same symmetry equations are solved and checked; they do not change the physical model." -> "Hermitian -> False 保持两个 hopping 方向相互独立，基本范例给出了相应非厄米矩阵。KernelMethod 和 ValidationLevel 只改变同一组对称方程的求解与检查方式，不改变物理模型。"

  ,"prepares the orbital centers and symmetry matrices needed to symmetrize a Wannier90 HR model." -> "准备对 Wannier90 HR 模型作对称化所需的轨道中心和对称矩阵。"
  ,"The lattice, Wyckoff positions, symmetry operations, and local orbitals determine the Wannier centers wcc and the full matrices DR returned by this compatibility interface." -> "晶格、Wyckoff 位置、对称操作和局域轨道共同决定这个兼容接口返回的 Wannier 中心 wcc 与完整矩阵 DR。"
  ,"The orbital order must match the external electronic-structure convention. At present only the documented VASP order is supported. The old MagneticTB backend is not loaded." -> "轨道顺序必须与外部电子结构程序的约定一致。目前只支持文档规定的 VASP 顺序，而且不会加载旧版 MagneticTB。"

  ,"returns the translation-periodic Hamiltonian of shell n before point-group, magnetic, or Hermitian constraints are imposed." -> "返回第 n 键层只满足平移周期性、尚未施加点群、磁对称或厄米约束的 Hamiltonian。"
  ,"Every directed hopping matrix element is independent, although the bond geometry and Bloch phases still come from the lattice prepared by init." -> "每个有向 hopping 矩阵元都相互独立，但键几何和 Bloch 相位仍取自 init 建立的晶格。"
  ,"This is mainly useful for understanding how symmetry reduces the parameter space. Use symham to construct the physical symmetry-constrained model." -> "这个函数主要用于理解对称性如何缩小参数空间。真正构造满足物理对称性的模型时，请使用 symham。"

  ,"For an ordered basis {f1,...,fn}, column j contains the coefficients of the transformed state g fj in that same basis: g fj = Sum_i fi D_ij(g)." -> "对有序基底 {f1,...,fn}，第 j 列给出变换后状态 g fj 在同一基底中的展开系数，即 g fj = Sum_i fi D_ij(g)。"
  ,"The identity and a twofold rotation about z act on the ordered {px,py} basis as follows. C2z sends both px and py to their negatives, so its matrix is -IdentityMatrix[2]:" -> "单位元和绕 z 轴的二重转动在有序基底 {px,py} 上的作用如下。C2z 使 px、py 都变号，因此其矩阵为 -IdentityMatrix[2]："
  ,"The full unitary tetragonal C4 point group from gray group 75 acts on the ordered {px,py} basis as follows:" -> "灰群 75 的完整幺正四方 C4 点群在有序基底 {px,py} 上的作用如下："
  ,"In the returned Association, wcc contains the fractional Wannier centers, DR contains the full symmetry matrices in the chosen orbital order, and symmetry lists the corresponding operations in the same order." -> "返回的 Association 中，wcc 是分数坐标 Wannier 中心，DR 是所选轨道顺序下的完整对称矩阵，symmetry 则以同样顺序列出与这些矩阵对应的操作。"
  ,"Every matrix acts in the orbital order displayed by orbitalTable[]. The function only reads the representation prepared by init or initfromrep." -> "每个矩阵都按 orbitalTable[] 显示的轨道顺序作用。函数只读取 init 或 initfromrep 已建立的表示。"
  ,"Initialize graphene and display the matrices of four operations that generate its complete symmetry group:" -> "初始化石墨烯，并显示能够生成其完整对称群的四个操作矩阵："
  ,"Each translation R is paired with a real-space hopping matrix H(R). Its row and column order is the Wannier-orbital order used when the HR file was written." -> "每个平移 R 对应一个实空间 hopping 矩阵 H(R)。矩阵行列顺序就是写出 HR 文件时采用的 Wannier 轨道顺序。"
  ,"Each row gives a representative position. Symbols x, y, and z are free Wyckoff coordinates; the moment entry gives the directions allowed by the site symmetry. The multiplicity is the number of equivalent positions in the conventional cell." -> "每行给出一个代表位置。x、y、z 是自由 Wyckoff 坐标，磁矩一栏给出点群允许的方向；multiplicity 表示常规晶胞中等价位置的数目。"
  ,StringJoin[
    "The two bands meet at K. plotRange only changes the displayed energy range; it does not ",
    "change the energies or the position of the crossing."
  ] -> "图中两条能带在 K 点相交。plotRange 只改变显示的能量范围，不会改变能量本征值或交叉位置。"
  ,StringJoin[
    "The header gives the number of Wannier orbitals and translation vectors. Each hopping ",
    "entry then gives a translation vector, two orbital indices, and the real and imaginary ",
    "parts of the matrix element."
  ] -> StringJoin[
    "文件头给出 Wannier 轨道数和平移矢量数。后面的每条跃迁数据依次给出平移矢量、两个轨道编号，以及矩",
    "阵元的实部和虚部。"
  ]
  ,"This tutorial uses graphene, MoS2, and several magnetic topological models to show the full MagneticTB workflow. Each calculation starts from the lattice, Wyckoff positions, symmetry, and local orbitals, constructs the Hamiltonian, and then evaluates bands, constant-energy contours, a Chern diagnostic, or Wannier90 data. The models reproduce the physical tasks collected in the original GeneralExamples.nb with the MagneticTB 2.0 interface." -> "本教程用石墨烯、MoS2 和几个磁性拓扑模型说明 MagneticTB 的完整用法。每个计算都从晶格、Wyckoff 位置、对称性和局域轨道出发，先构造 Hamiltonian，再计算能带、等能线、Chern 诊断或 Wannier90 数据。这些模型用 MagneticTB 2.0 接口重现原 GeneralExamples.nb 中的物理任务。"
  ,"MagneticTBOld` never loads or calls MagneticTB`, and MagneticTB` never loads or calls MagneticTBOld`. They do not share the current model, parameter symbols, or any intermediate result. To compare them, run the same physical input in two separate fresh kernels and compare the resulting Hamiltonians and bands." -> "MagneticTBOld` 与 MagneticTB` 不会互相加载或调用，也不共享当前模型、参数符号或任何中间结果。比较新旧实现时，应在两个全新的内核中分别运行同一物理输入，再比较所得 Hamiltonian 和能带。"

  ,"After applying rules, the function samples the same reciprocal-fractional path convention used by bandplot. Each sampled momentum occupies one row of the file, and each band occupies one column." -> "代入 rules 后，函数按 bandplot 相同的倒空间分数坐标路径取样。文件中每个采样动量占一行，每条能带占一列。"
  ,"Initialize graphene, construct its Hamiltonian and reciprocal-space path, and choose a portable temporary CSV file:" -> "初始化石墨烯，构造 Hamiltonian 和倒空间路径，并选择一个可移植的临时 CSV 文件："
  ,"Export the sampled band energies. Each row of the CSV file is one momentum point and each column is one band:" -> "导出取样得到的能带能量。CSV 文件中每行对应一个动量点，每列对应一条能带："
  ,"Import the CSV file and plot the two graphene bands:" -> "载入 CSV 文件，画出石墨烯的两条能带："
  ,"The two rows of orbitalTable are the two sites generated from the database seed. They are also the two basis states used by the Hamiltonian and by every displayed symmetry matrix." -> "orbitalTable 的两行是由数据库种子生成的两个格点，也就是 Hamiltonian 和所有对称矩阵共同采用的两个基态。"
  ,"Induced starts from the site-symmetry representation of the first row and constructs the partner state with a coset representative. Although its intermediate construction differs from DirectProduct, the two Hamiltonians span the same allowed matrix space and their band degeneracies agree." -> "Induced 从第一行格点的位置稳定子群表示出发，用陪集代表元构造伙伴格点的状态。虽然中间构造与 DirectProduct 不同，两者所得 Hamiltonian 张成相同的允许矩阵空间，能带简并也一致。"
  ,"Local states may carry different SO(2) weights m. A matrix element can be nonzero only when its total weight is conserved, so the continuous constraint can eliminate forbidden onsite and hopping terms before the remaining finite symmetries are imposed." -> "局域态可以带有不同的 SO(2) 权 m。只有总权守恒的矩阵元才可能非零，因此连续对称约束可以先消去禁戒的 onsite 和 hopping 项，再施加其余有限对称性。"
  ,"Each of the two sites carries two local states with m=+1/2 and m=-1/2. The complete finite matrices and both continuous site matrices are supplied explicitly, so the full Hamiltonian has four bands." -> "两个格点各自带有 m=+1/2 和 m=-1/2 两个局域态。显式输入完整有限群矩阵及两个格点上的连续表示矩阵后，总 Hamiltonian 有四条能带。"
  ,"Continuous and finite symmetry conditions" -> "连续与有限对称条件"
  ,"Representative MagneticTB Models" -> "MagneticTB 典型模型"
  ," — generate a finite group from an explicit product." -> " — 按显式给定的乘法规则生成有限群。"
  ,"length" -> "长度"
  ,"fractional position and magnetic moment" -> "分数坐标与磁矩"
  ,"symmetry operations" -> "对称操作"
  ,"local orbital basis" -> "局域轨道基底"
  ,"The following models illustrate unequal local dimensions, several Wyckoff orbits, spinful magnetic symmetry, induced representations, and antiunitary transport between sites. Every section starts from complete physical input and displays the Hamiltonian, orbital order or hopping interpretation, symmetry matrices, and bands." -> "下列模型分别展示不同局域维数、多个 Wyckoff 轨道、含自旋的磁对称性、诱导表示，以及反幺正操作在格点之间的搬运。每一节都从完整物理输入出发，直接显示 Hamiltonian、轨道顺序或 hopping 含义、对称矩阵和能带。"

  ,"A continuous internal C_infinity symmetry is specified by a one-parameter matrix D(theta) acting on the local states. MagneticTB obtains its Hermitian generator from the derivative at theta=0 and requires each hopping matrix to satisfy the corresponding continuous constraint. The finite space-group operations are imposed afterwards. This continuous operation acts only on the local states: it does not move atoms, generate a finite group, or create new bonds." -> "连续内部 C_infinity 对称性由作用在局域态上的单参数矩阵 D(theta) 指定。MagneticTB 在 theta=0 处求导得到厄米生成元，并要求每个 hopping 矩阵满足相应连续约束，随后再施加有限空间群操作。这个连续操作只作用于局域态，不移动原子，不生成有限群，也不产生新的键。"
  ,"For 14.1.2.3.L.1, the eight finite operations and their exact double-valued matrices are listed explicitly below. They were generated with SpinLayerCorepresentations and are included in the example so that it can run without loading another paclet. Starting from 000, the finite spatial operations generate the two sites {0,0,0} and {1/2,1/2,0}; the operations that leave 000 fixed form its site-symmetry group." -> "对 14.1.2.3.L.1，下面显式列出八个有限操作及其精确双值矩阵。这些数据由 SpinLayerCorepresentations 生成，并直接放入范例，因此运行时无需加载其他 paclet。从 000 出发，有限空间操作生成 {0,0,0} 与 {1/2,1/2,0} 两个格点；保持 000 不变的操作组成 000 位置的稳定子群。"
  ,"The next cell puts the finite operations in the order required by symminformation, generates the two sites, finds the site-symmetry group of 000, and adds the pure internal C_infinity rotation. It only prepares the symmetry data; the Hamiltonian is calculated in the following sections." -> "下一段代码按 symminformation 所需顺序排列有限操作，生成两个格点，找出 000 位置的稳定子群，并加入纯内部 C_infinity 转动。这里仅准备对称性数据，Hamiltonian 在后续小节中计算。"
  ,"For the continuous symmetry, differentiating D(theta) at theta=0 gives the generator J on each site. The Hamiltonian must satisfy [J,H(k)]=0. The finite operations are checked separately: unitary operations act by U H(k) U-dagger, while antiunitary operations act by U Conjugate[H(k)] U-dagger. The example also checks the exact group multiplication and Hermiticity." -> "对连续对称性，在 theta=0 处对 D(theta) 求导可得到每个格点上的生成元 J，Hamiltonian 必须满足 [J,H(k)]=0。有限操作另行处理：幺正操作按 U H(k) U-dagger 作用，反幺正操作按 U Conjugate[H(k)] U-dagger 作用。本例还检查精确群乘法关系和厄米性。"

  ,StringJoin[
    "After symham has obtained the hopping terms, use hop to write them to a Wannier90 HR ",
    "file. The beginning of the file is shown below:"
  ] -> StringJoin[
    "symham 求出跃迁项后，可以用 hop 写成 Wannier90 的 HR ",
    "文件。下面显示文件开头的内容："
  ]
  ,"CartesianCoordinates changes only the definition of the momentum components. The Hamiltonian in Cartesian momentum coordinates is:" -> "CartesianCoordinates 只改变动量分量的定义。用笛卡尔动量坐标写出的 Hamiltonian 为："
  ,"Changing the symmetry in version 2.0" -> "在 2.0 中改变模型对称性"
  ,StringJoin[
    "To lower the symmetry, change the model inputs and run init again before calculating a ",
    "new Hamiltonian. brokenSymmetryInitRules gives the retained subgroup and the resulting ",
    "Wyckoff positions. Inspect these rules and pass them to init@@rules. The symmetry is ",
    "set by init, not by a symmetryset option in symham."
  ] -> StringJoin[
    "如果要降低对称性，需要修改模型输入，重新 init 后再求 ",
    "Hamiltonian。brokenSymmetryInitRules 可以给出保留的子群及新的 ",
    "Wyckoff 位置，查看这些输入后，用 init@@rules 重新初始化。对称性由 init ",
    "设置，不能在 symham 中用 symmetryset 改变。"
  ]
  ,"The BNS 184.196 example uses the same local spinor basis but a different magnetic group. The displayed block contains the momentum-dependent matrix elements, and the short path around the crossing resolves the symmetry-protected nodal connection." -> "BNS 184.196 使用相同的局域旋量基底，但磁群不同。显示的矩阵块给出含动量的矩阵元，穿过交叉点附近的短路径展示对称性保护的节点连接。"
  ,StringJoin[
    "The parameters in this Weyl model and the following nodal-line and C4 models illustrate ",
    "symmetry-allowed band structures. They are the values used in the MagneticTB 1.0 ",
    "examples, not fitted parameters for a particular material."
  ] -> StringJoin[
    "这个 Weyl 模型，以及后面的节点线和 C4 模型，使用 MagneticTB 1.0 ",
    "例子的参数展示对称性允许的能带。这些参数不是对某个具体材料的拟合结果。"
  ]
  ,StringJoin[
    "The symmetry databases can also be used on their own. showMSGWyckoff lists the ",
    "representative positions and allowed moments of a magnetic space group. mlgop and mrgop ",
    "give the operations of magnetic layer groups and magnetic rod groups, respectively."
  ] -> StringJoin[
    "对称性数据库也可以单独查询。showMSGWyckoff ",
    "列出磁空间群的代表位置和允许的磁矩，mlgop 和 mrgop 分别给出磁层群和磁杆群的对称操作。"
  ]
  ,StringJoin[
    "For the spinful onsite Hamiltonian, first place all spin-up orbitals before the ",
    "spin-down orbitals. Specify the reduced symmetry in init before calculating the ",
    "Hamiltonian:"
  ] -> StringJoin[
    "考虑自旋后，先把所有自旋向上的轨道排在前面，自旋向下的轨道排在后面。所用的较低对称性在 init ",
    "中给出，然后求 onsite Hamiltonian："
  ]

  ,StringJoin[
    "To construct a tight-binding model, first specify the lattice, atomic positions, ",
    "symmetry operations, and local orbitals. With these inputs, MagneticTB finds the ",
    "equivalent atoms and the symmetry-allowed hopping terms. We use graphene to show how to ",
    "obtain the Hamiltonian and plot its bands."
  ] -> StringJoin[
    "要构造紧束缚模型，需要给出晶格、原子位置、对称操作和局域轨道。有了这些输入，MagneticTB ",
    "就可以生成等价原子，求出对称性允许的跃迁。下面以石墨烯为例，介绍如何构造 Hamiltonian ",
    "并画出能带。"
  ]
  ,StringJoin[
    "init sets up the crystal and its local orbitals. Here lattice gives the three lattice ",
    "vectors as rows, and lattpar supplies their parameter values. wyckoffposition gives one ",
    "representative position and its magnetic moment for each inequivalent Wyckoff orbit. ",
    "symminformation gives the symmetry operations, and basisFunctions lists the orbitals in ",
    "the same order as the positions."
  ] -> StringJoin[
    "有了对称操作，就可以使用 init 设置晶格、原子位置和基函数。其中 lattice ",
    "的三行分别是三个晶格基矢，lattpar 给出晶格参数；wyckoffposition 给出各个不等价 ",
    "Wyckoff 位置的代表原子及磁矩；symminformation ",
    "是对称操作，basisFunctions 按相同顺序列出各位置的局域轨道。"
  ]
  ,StringJoin[
    "The two carbon atoms in graphene are symmetry equivalent. Enter only {1/3,2/3,0}; the ",
    "other atom at {2/3,1/3,0} is generated by symmetry. Keep one pz orbital on each carbon. ",
    "The lattice parameter c separates the periodically repeated graphene layers."
  ] -> StringJoin[
    "石墨烯的两个碳原子是对称等价的，只需输入 {1/3,2/3,0}，程序就会生成另一个位置 ",
    "{2/3,1/3,0}。每个碳原子保留一个 pz 轨道，晶格参数 c 用来分隔周期重复的石墨烯层。"
  ]
  ,"msgop[gray[191]] already returns all operations of the gray group, so GenerateSymmetryGroup -> False tells init to use this complete list as given. The following four operations are sufficient to generate the same group; this information is useful when one wants to provide generators instead of the full list:" -> "msgop[gray[191]] 已返回灰群的全部操作，因此 GenerateSymmetryGroup -> False 表示 init 直接使用这个完整列表。下面四个操作足以生成同一个群；只想输入生成元时可以参考这组结果："
  ,"Before calculating a Hamiltonian, first check the atoms and orbital order produced by init. orbitalTable[] lists the basis used for every Hamiltonian and representation matrix. The two rows below are the two equivalent carbon atoms; their coordinates are fractional coordinates of the primitive lattice, and each carries one pz orbital. Spin is None because this model is spinless." -> "计算 Hamiltonian 前，先检查 init 生成的原子和轨道顺序。orbitalTable[] 列出所有 Hamiltonian 与表示矩阵共同采用的基底。下面两行对应两个等价碳原子，坐标是原胞分数坐标，每个原子带一个 pz 轨道；模型无自旋，所以 Spin 显示为 None。"
  ,StringJoin[
    "Hoppings are numbered by distance. Number 1 is the onsite term, number 2 is the nearest ",
    "neighbor, and number 3 is the second neighbor. For graphene, showbonds[2] displays the ",
    "three nearest-neighbor carbon-carbon bonds and their reverse directions. No hopping ",
    "parameters need to be specified to draw these bonds."
  ] -> StringJoin[
    "跃迁按距离编号，1 表示 onsite 项，2 表示最近邻，3 ",
    "表示第二近邻。对于石墨烯，showbonds[2] 画出三条最近邻碳碳键及其反向。查看这些键时还不需要给",
    "跃迁参数赋值。"
  ]
  ,StringJoin[
    "Add the contributions to obtain the Hamiltonian. Here we keep the onsite, ",
    "nearest-neighbor, and second-neighbor terms:"
  ] -> "把所需的项相加，就得到 Hamiltonian。这里保留 onsite、最近邻和第二近邻跃迁："
  ,StringJoin[
    "The path is given in reciprocal fractional coordinates. Each entry contains a segment ",
    "and its endpoint labels. bandManipulate draws the bands along this path and provides a ",
    "slider for each parameter. The sliders initially equal zero, and no extra energy shift ",
    "is added."
  ] -> StringJoin[
    "能带路径使用倒空间分数坐标，每一项给出一段路径及其端点名称。bandManipulate ",
    "沿这条路径画出能带，并为每个参数提供一个滑块。滑块初值为零，程序不会额外移动能量零点。"
  ]
  ,"The symmetry operations act both on momentum and on the orbital basis. showSymmetryRepresentations[] displays the matrices for a generating set. Whether an operation is unitary or antiunitary is already recorded in symminformation, so the table displays only its matrix U. A unitary operation satisfies U H(k) U^dagger = H(g k), whereas an antiunitary operation satisfies U Conjugate[H(k)] U^dagger = H(g k)." -> "对称操作同时作用于动量和轨道基底。showSymmetryRepresentations[] 显示一组生成元的矩阵。操作是幺正还是反幺正已经记录在 symminformation 中，因此表格只显示矩阵 U。幺正操作满足 U H(k) U^dagger = H(g k)，反幺正操作满足 U Conjugate[H(k)] U^dagger = H(g k)。"
  ,"In the ordering returned by msgop[gray[191]], operation 25 is time reversal. Its displayed matrix is only the linear matrix U; the operation data already state that it is antiunitary, so the action on the Hamiltonian includes complex conjugation." -> "在 msgop[gray[191]] 返回的顺序中，第 25 个操作是时间反演。显示内容只有线性矩阵 U；操作数据已经说明它是反幺正的，因此作用于 Hamiltonian 时还要包含复共轭。"
  ,"The following calculation checks every symmetry operation at the generic momentum {Pi/7,Pi/11,Pi/13}. The check is carried out for every independent term of the symbolic Hamiltonian, so it does not depend on a particular choice of hopping parameters. A zero residual means that the corresponding unitary or antiunitary covariance relation is satisfied." -> "下面在一般动量 {Pi/7,Pi/11,Pi/13} 检查全部对称操作。检查逐个作用于符号 Hamiltonian 的独立项，因此不依赖某一组 hopping 参数。余量为零表示相应的幺正或反幺正协变关系成立。"
  ,"Checking only the Gamma-point is not sufficient: at Gamma the momentum transformation is invisible. A complete check must also use generic momenta and non-Gamma points, including the reciprocal-lattice phase needed when a transformed momentum differs by a reciprocal vector. The examples for the seven crystal systems are tested in this stronger way." -> "只检查 Gamma 点是不够的，因为在 Gamma 点看不出动量如何变换。完整检查还必须取一般动量和非 Gamma 高对称点；若变换后的动量相差一个倒格矢，还要包含相应的轨道位置相位。七晶系范例都按这种方式检查。"
  ,StringJoin[
    "DirectProduct uses the symmetry action on the orbitals supplied at each site. Induced ",
    "starts from local states at one representative site and generates the symmetry-related ",
    "states. For example, in the cubic p-orbital model, DirectProduct requires px, py, and ",
    "pz explicitly, whereas Induced can start from px at the representative site. The ",
    "operation relating two sites may also include time reversal."
  ] -> StringJoin[
    "DirectProduct 使用各位置所给轨道的对称变换，Induced ",
    "则从代表位置的局域态出发，生成对称相关的态。例如，立方晶系的 p 轨道模型用 ",
    "DirectProduct 时需要给出 px、py、pz，而用 Induced 时可以从代表位置的 ",
    "px 出发。连接两个位置的操作也可以含有时间反演。"
  ]
  ,"InitialBondShells is 10 by default; the graphene calculation above deliberately requested only 3. Asking for shell 4 must therefore stop and request a new init with a larger value. MagneticTB does not search farther or change the null-space method without being asked. If a new initialization fails, no earlier model is used in its place." -> "InitialBondShells 默认为 10；上面的石墨烯计算特意只请求 3 层。因此请求第 4 层时必须停止，并要求用更大的 InitialBondShells 重新 init。MagneticTB 不会擅自搜索更远键层，也不会自行更换零空间方法。若新的初始化失败，也不会改用先前模型。"
  ,"Initialize the three-shell graphene model again, then request an unprepared shell. The call stops with $Failed and leaves the requested model boundary explicit." -> "重新初始化三键层石墨烯模型，再请求一个未准备的键层。调用会返回 $Failed，明确保留用户指定的计算范围。"
  ,"The outer basisFunctions list must have one entry per Wyckoff seed. This input has one seed but two basis lists, so initialization stops and the preceding graphene model is not reused." -> "basisFunctions 外层列表必须与 Wyckoff 种子一一对应。下面只有一个种子却给出两组基函数，因此初始化停止，也不会继续使用此前的石墨烯模型。"
  ,"A model-dependent function called after this failed init reports immediately that no model is available." -> "在这次 init 失败后调用依赖模型的函数，会立即报告当前没有可用模型。"
  ,"The remaining tutorials change one physical ingredient at a time: CsCl introduces rectangular s-p hopping, the cubic p-orbital example compares DirectProduct and Induced, E+PT illustrates antiunitary induction, and the spin-space-group examples distinguish collinear and noncollinear spin actions." -> "后续教程每次只改变一个物理要素：CsCl 引入矩形 s-p hopping；立方 p 轨道范例比较 DirectProduct 与 Induced；E+PT 说明反幺正诱导；自旋空间群范例则区分共线和非共线自旋作用。"
  ,"After restarting Mathematica, search for MagneticTB in the Documentation Center, or open the guide directly with the following URI. The package also contains the old implementation under MagneticTBOld`; use a separate fresh kernel when comparing old and new calculations." -> "重启 Mathematica 后，可在帮助中心搜索 MagneticTB，也可用下面的 URI 直接打开指南。包内还包含 MagneticTBOld` 形式的原实现；比较新旧计算时，请分别使用两个全新的内核。"
  ,"Set pacletArchive to the downloaded MagneticTB paclet and evaluate the following cell. PacletFind shows the installed version. If PacletInstall reports samevers, the same version is already present; install a newer archive, or uninstall that version explicitly before installing it again. Restart the kernel after replacing an installed version." -> "把 pacletArchive 设为下载的 MagneticTB paclet 路径，再运行下面代码。PacletFind 会显示已安装版本。若 PacletInstall 报告 samevers，说明同一版本已经安装；可以安装更新版本，或先显式卸载再重装。替换已安装版本后请重启内核。"

  ,"MagneticTB 2.0 includes the original implementation so that established notebooks remain usable and old/new results can be compared directly. The old code is exposed through names ending in old and is loaded from MagneticTBOld`, completely separately from the new package." -> "MagneticTB 2.0 同时包含原来的实现，使已有 notebook 可以继续运行，也便于直接比较新旧结果。旧版函数名以 old 结尾，并从 MagneticTBOld` 加载，与新版完全隔离。"
  ,"Plot the bands with the original implementation" -> "用原来的实现绘制能带"

  ,"The cubic rotations of Pm-3m mix px, py, and pz. DirectProduct must therefore receive the complete three-orbital basis on the representative site. Induced starts with only px at {1/2, 0, 0}; the site-symmetry representation and coset transport generate the corresponding py and pz states on the other two symmetry-equivalent sites. The two orbital tables make this reduction explicit." -> "Pm-3m 的立方旋转会混合 px、py、pz，因此 DirectProduct 必须在代表点输入完整的三个轨道。Induced 只从 {1/2, 0, 0} 上的 px 出发；位置稳定子群表示与陪集搬运会在另外两个对称等价格点生成相应的 py、pz 态。两个轨道表直接展示了输入规模的缩减。"
  ,"C4spin rotates spin by 90 degrees without rotating real space. tauHalfT combines a half translation with an antiunitary spin action. With GenerateSymmetryGroup -> True, these two operations generate the complete finite spin-space group used to constrain the Hamiltonian." -> "C4spin 只把自旋旋转 90 度而不转动实空间；tauHalfT 把半平移与反幺正自旋作用组合起来。设 GenerateSymmetryGroup -> True 后，这两个操作生成约束 Hamiltonian 所需的完整有限自旋空间群。"
  ,"Place one s orbital at the origin of a simple-cubic lattice. The gray Pm-3m group keeps this site fixed, and the six bonds in shell 2 connect it to the nearest cells along plus or minus x, y, and z. The resulting one-band Hamiltonian contains one onsite energy and the usual sum of three cosine dispersions." -> "在简单立方晶格原点放置一个 s 轨道。灰 Pm-3m 群保持该点不变，第 2 键层的六条键沿正负 x、y、z 方向连接最近晶胞。所得单带 Hamiltonian 包含一个 onsite 能量和三个余弦色散之和。"
  ,"Related applications" -> "相关应用"
  ,"The first Wyckoff orbit carries dxz and dyz, while the second carries px, py, and pz. Because both seeds generate equivalent partners, the complete basis has ten states. The orbital table fixes their order; the displayed 3 by 3 corner is a readable part of the same ten-band Hamiltonian used for the band plot." -> "第一个 Wyckoff 轨道带 dxz、dyz，第二个带 px、py、pz。两个种子都会生成一个等价伙伴，因此完整基底共有十个态。轨道表固定它们的顺序；显示的 3×3 左上角只是同一个十带 Hamiltonian 中便于阅读的一块，能带图使用完整矩阵。"
  ,"The General Examples tutorial develops graphene, the three-band MoS2 model, magnetic Weyl and nodal-line models, and a C4 topological-insulator model. The Wyckoff tutorial gives ten DirectProduct/Induced pairs drawn from all seven crystal systems." -> "“通用范例”教程包含石墨烯、三带 MoS2、磁性 Weyl 与节点线模型以及 C4 拓扑绝缘体。Wyckoff 教程给出覆盖七个晶系的十组 DirectProduct/Induced 对照。"
  ,"The identity E cannot move the site at x=1/4 to its partner at x=3/4. In this model the required coset representative is PT, so the second site is reached by an antiunitary operation. symminformation already records this fact; the representation table therefore shows only its matrix U. Its action on the Hamiltonian is U Conjugate[H(k)] U-dagger." -> "单位元 E 不能把 x=1/4 的格点送到 x=3/4 的伙伴点。在这个模型中，所需陪集代表元是 PT，因此第二个格点由反幺正操作得到。symminformation 已记录反幺正性，所以表示表只显示矩阵 U；它在 Hamiltonian 上的作用为 U Conjugate[H(k)] U-dagger。"
  ,"The reference site has C3 symmetry and carries the one-dimensional character {1,omega,omega^2}. An antiunitary coset maps it to the partner site. Antiunitarity complex-conjugates the local character during this transport, so C3 is represented by omega on one site and omega^2 on the other." -> "参考点具有 C3 对称性，并带一维特征标 {1,omega,omega^2}。一个反幺正陪集把它映到伙伴点；搬运过程中反幺正性会复共轭局域特征标，因此 C3 在一个格点上表示为 omega，在另一个格点上表示为 omega^2。"

  ,"The ten examples below cover all seven crystal systems. In each case the representative atomic position is read directly from the bundled magnetic Wyckoff database, and the same physical orbit is used to construct both DirectProduct and Induced models. Each branch shows its complete init input, Hamiltonian, representation matrices, and bands." -> "下面十个范例覆盖七个晶系。每个范例都从随包磁性 Wyckoff 数据库直接读取代表原子位置，并用同一个物理轨道分别构造 DirectProduct 和 Induced 模型。两条分支都给出完整 init 输入、Hamiltonian、表示矩阵和能带。"
  ,"For this Wyckoff position, the unitary subgroup leaves the representative site in its own orbit, while an antiunitary operation generates the second site. The position shown below is read from the magnetic Wyckoff database. We use it first in DirectProduct mode and then in Induced mode." -> "对这个 Wyckoff 位置，幺正子群只生成代表点本身，第二个格点由一个反幺正操作生成。下面的位置直接读取自磁性 Wyckoff 数据库，先用于 DirectProduct 模式，再用于 Induced 模式。"
  ,"How the two constructions are compared" -> "如何比较两种构造"
  ,"For each DirectProduct/Induced pair, the comparison asks whether both methods span the same real space of allowed Hamiltonians and give the same bands; the names of their parameters need not agree. Their representation or corepresentation matrices must obey the group multiplication law. The Hamiltonian is also checked at six generic momenta and at every non-Gamma little group on the {0,Pi/2,Pi}^3 grid. When two momenta differ by a reciprocal vector, the required orbital-position phase is included. A commutator at Gamma alone cannot establish the full space-group symmetry." -> "对每一组 DirectProduct/Induced，首先比较两种方法是否张成同一个实 Hamiltonian 空间并给出相同能带；它们的参数名称不必相同。表示或共表示矩阵还必须满足群乘法关系。Hamiltonian 的对称性在六个一般动量以及 {0,Pi/2,Pi}^3 网格中全部非 Gamma 小群点上检查；当两个动量相差倒格矢时，还要包含轨道位置产生的相位。只在 Gamma 点计算对易子不能证明完整空间群对称性。"
  ,StringJoin[
    "To include spin-orbit coupling, replace pz by {pzup,pzdn} on each carbon atom. There ",
    "are now four basis states. The gray group includes time reversal, and symham gives the ",
    "spin-dependent terms allowed by this symmetry. We use this Hamiltonian below to plot ",
    "bands and export Wannier90 data."
  ] -> StringJoin[
    "考虑自旋轨道耦合时，把每个碳原子的 pz 换成 {pzup,pzdn}，两个碳原子共有四个基函数。这里使",
    "用的灰群含有时间反演，symham 会给出该对称性允许的自旋相关项。下面用这个 Hamiltonian ",
    "画能带，并导出 Wannier90 数据。"
  ]

  ,"Use this function when another program has already assembled the lattice, Wyckoff positions, symmetry operations, and local orbitals of a model. For an ordinary notebook calculation, init or initfromrep is simpler." -> "如果晶格、Wyckoff 位置、对称操作和局域轨道是由另一个程序生成的，可以用这个函数先整理并检查这些输入。平常在 notebook 中从头建模时，直接使用 init 或 initfromrep 更方便。"
  ,"The string keys of input describe the same physical quantities used by init. The result lets the calling program inspect the generated sites, representation matrices, and bond geometry. It does not calculate any hopping parameter or replace the model currently used by symham." -> "input 中的字符串键仍对应 init 的各项物理输入。调用程序可以从结果中查看生成的格点、表示矩阵和键的几何。这里不计算 hopping 参数，也不改变 symham 当前使用的模型。"
  ,"After init finishes, first inspect the generated atoms, the orbital order, and the neighbour shells. The hopping parameters are not needed at this stage; symham determines them later, one shell at a time." -> "init 完成后，先检查生成的原子、轨道顺序和近邻键层。此时还没有必要求 hopping 参数；随后用 symham 按键层逐一求出即可。"
  ,"By default init finds the first six complete bond shells. If a model needs a farther shell, choose a larger InitialBondShells value and run init again; symham itself does not search for new bonds." -> "init 默认找出前六个完整键层。若模型需要更远的 hopping，应增大 InitialBondShells 后重新运行 init；symham 本身不再向外搜索新键。"
  ,"Every call to init starts a new model. If the input is inconsistent, the calculation stops and asks the user to correct it; MagneticTB does not silently continue with the preceding model." -> "每次运行 init 都是在建立一个新模型。若输入彼此不相容，计算会停下来提示用户检查；MagneticTB 不会沿用上一个模型。"
  ,"Number of complete bond shells found during initialization. Shell 1 contains the zero-distance onsite terms; shell 2 is the first nonzero-distance neighbour shell." -> "初始化时要找出的完整键层数。第 1 层是零距离 onsite 项，第 2 层是第一个非零距离近邻层。"
  ,"Call init first. For bond shell n, symham starts from the real-space hoppings in that shell, relates symmetry-equivalent matrix elements, and keeps the independent amplitudes. Each hopping is then multiplied by its Bloch phase and placed in the Hamiltonian." -> "先运行 init。对第 n 个键层，symham 把对称等价的实空间 hopping 联系起来，只留下彼此独立的振幅；再给每条 hopping 乘上相应的 Bloch 相位，写入 Hamiltonian。"
  ,"The result is Hermitian by default. The KernelMethod and ValidationLevel options control how the common zero space is found and checked; they do not change the symmetry or orbitals of the model. CartesianCoordinates -> False keeps reciprocal-fractional momentum components." -> "默认得到厄米 Hamiltonian。KernelMethod 和 ValidationLevel 只决定怎样求共同零空间以及怎样检查结果，不改变模型的对称性和轨道。CartesianCoordinates -> False 保留倒空间分数坐标下的动量分量。"
  ,"symham uses the atoms and bonds already found by init; it does not initialize the model again. If shell n was not included, rerun init with InitialBondShells -> n." -> "symham 直接使用 init 已经找出的原子和键，不会再次初始化。若第 n 层不在原来的范围内，请设定 InitialBondShells -> n 后重新运行 init。"
  ,"Add onsite, nearest-neighbour, and next-nearest-neighbour shells to obtain the Bloch Hamiltonian used below:" -> "把 onsite、最近邻和次近邻三部分相加，得到下面使用的 Bloch Hamiltonian："
  ,"Most users can inspect a model with orbitalTable, showbonds, and showSymmetryRepresentations. CurrentModelSession is useful when a program needs all of that information together: the crystal model, the full symmetry matrices, the bond shells found by init, and the shell results already calculated by symham." -> "一般查看模型时，使用 orbitalTable、showbonds 和 showSymmetryRepresentations 就够了。若要在程序中一次取得更完整的信息，可以使用 CurrentModelSession：其中包括晶体模型、全部对称矩阵、init 找到的键层，以及 symham 已经算出的各层结果。"
  ,"Read the representation mode, generated site orbits, local bases, and available shell numbers of the graphene model. Merely viewing them does not calculate another Hamiltonian:" -> "下面读出石墨烯采用的表示方式、生成的格点轨道、局域基底和可用键层编号。查看这些信息不会另外计算 Hamiltonian："
  ,"There are two ways to use hop. In hop[n,rules] and hop[{n1,n2,...},rules], first calculate the required shells with symham. The integer form exports shells 1 through n, while the list form exports only the listed shells. The real-space hopping matrices already obtained by symham are written directly." -> "hop 有两种用法。使用 hop[n,rules] 或 hop[{n1,n2,...},rules] 时，要先用 symham 算出所需键层。整数 n 表示导出第 1 至第 n 层；列表只导出列出的键层。symham 已经求出的实空间 hopping 矩阵会被直接写出。"
  ,"The matrix form hop[H,rules] extracts the coefficients of the finite exponential phases in H. The examples below obtain H from symham, so the model construction remains visible and reproducible. Set \"hrExport\" to a path to write wannier90_hr.dat; with None, the text is shown in the notebook. Use \"wcc\" -> Automatic when the orbital centers should come from init." -> "矩阵形式 hop[H,rules] 提取 H 中有限个指数相位的系数。下面的例子都从 symham 得到 H，使模型构造过程保持可见且可复现。把 \"hrExport\" 设为路径即可写出 wannier90_hr.dat；设为 None 时文本显示在 notebook 中。若轨道中心取自 init，使用 \"wcc\" -> Automatic。"
  ,"Use this function when a term in the model preserves only part of the original symmetry. The selected operations generate the symmetry subgroup that is to remain, while all atoms are kept." -> "若模型中的某一项只保留原来的一部分对称性，可以使用这个函数。选中的操作生成要保留的子群，原模型中的原子则全部留下。"
  ,"The function only writes a new set of init rules; it does not change the model by itself. Check the new Wyckoff positions and symmetry operations, then evaluate init@@rules. Use {} when only the identity is to remain." -> "函数只给出一组新的 init 规则，本身不会改变当前模型。先检查新的 Wyckoff 位置和对称操作，再执行 init@@rules。若只保留单位元，输入 {}。"
  ," \[LongDash] read all information already available for the current model." -> " \[LongDash] 读取当前模型已有的全部信息。"
  ," \[LongDash] obtain the symmetry-allowed Hamiltonian of one bond shell." -> " \[LongDash] 求一个键层中对称性允许的 Hamiltonian。"
  ," \[LongDash] display the length and geometry of a bond shell." -> " \[LongDash] 显示一个键层的键长和几何关系。"
  ,"Functions from the original version" -> "原版函数"
  ,"The MagneticTB 1.0 functions are included under MagneticTBOld`, with old appended to their names. They can run earlier notebooks and provide a direct comparison with version 2.0. Load the 1.0 and 2.0 versions in separate kernels." -> "MagneticTB 1.0 版函数也随软件一同提供，位于 MagneticTBOld` 中，函数名末尾加 old。它们可继续运行 1.0 版 notebook，也可用来与 2.0 版比较。两个版本请分别在不同内核中载入。"
  ,"Check the neighbour shells" -> "检查近邻键层"
  ,"To see the momentum transformation explicitly, use the general point {Pi/7,Pi/11,Pi/13}. The following calculation applies every unitary and antiunitary operation to each independent term of the symbolic Hamiltonian. A zero residual means that the two sides of the corresponding covariance equation are the same." -> "为了真正看出动量怎样变换，下面取一般动量 {Pi/7,Pi/11,Pi/13}。程序把每个幺正和反幺正操作分别作用到符号 Hamiltonian 的各个独立项上。余量为零，表示相应协变关系的两边完全相同。"
  ,"The Gamma point alone cannot show how momentum is transformed. At a non-Gamma high-symmetry point, the transformed momentum may also differ by a reciprocal vector; the orbital positions then supply an additional phase. The tutorial for the seven crystal systems works through both general momenta and these non-Gamma little groups." -> "只看 Gamma 点看不出动量如何变化。在非 Gamma 高对称点，变换后的动量还可能相差一个倒格矢，此时必须补上轨道位置带来的相位。七晶系教程同时检查了一般动量和这些非 Gamma 小群。"
  ,StringJoin[
    "The following examples show how to construct models for graphene, MoS2, magnetic Weyl ",
    "points, nodal lines, and a C4 topological insulator. For each model, first specify the ",
    "lattice, atomic positions, symmetry, and orbitals. Then calculate the Hamiltonian, ",
    "assign its parameters, and plot the bands. Some examples also calculate a Chern number ",
    "or export Wannier90 hopping data."
  ] -> StringJoin[
    "下面介绍石墨烯、MoS2、磁性 Weyl 点、节点线和 C4 ",
    "拓扑绝缘体的模型。每个例子都先给出晶格、原子位置、对称操作和轨道，再求 ",
    "Hamiltonian、给参数赋值并画出能带。部分例子还介绍 Chern 数计算和 Wannier90 ",
    "跃迁数据的导出。"
  ]
  ,"The following ten examples cover the seven crystal systems. Each representative position is read directly from the magnetic Wyckoff database supplied with MagneticTB. The same orbit is then described in two ways: DirectProduct uses the action on a complete local basis, whereas Induced starts from the site-symmetry representation of one reference point. For both descriptions the complete init input, Hamiltonian, symmetry matrices, and bands are shown." -> "下面十个例子覆盖七个晶系。每个代表位置都直接从 MagneticTB 自带的磁性 Wyckoff 数据库读取。同一个轨道分别用两种方式描述：DirectProduct 使用完整局域基底上的作用；Induced 从一个参考点的位置稳定子群表示出发。两种方式都给出完整 init 输入、Hamiltonian、对称矩阵和能带。"
  ,"The parameter names in the two descriptions need not be the same. What matters is that the two sets of matrices describe the same allowed Hamiltonians and give the same bands. The comparison therefore checks both directions between the two matrix spaces. It also applies every symmetry at general momenta and at non-Gamma little-group points. When two momenta differ by a reciprocal vector, the phase from the orbital positions is included. A commutator evaluated only at Gamma would miss this information." -> "两种表示中的参数名不必相同，关键是它们给出同一个允许的 Hamiltonian 空间和相同的能带。因此比较时要确认两边的矩阵都能互相表示，还要在一般动量和非 Gamma 小群点逐个施加对称操作。若两个动量相差一个倒格矢，还必须计入轨道位置产生的相位；只在 Gamma 点算对易子会漏掉这些信息。"
  ,"Consider two local states that acquire different phases under a rotation by theta. A hopping between them is allowed only when their SO(2) weights are compatible. In initfromrep this continuous C_infinity symmetry is entered as a one-parameter matrix D(theta). MagneticTB differentiates D(theta) at theta=0, obtains the Hermitian generator, and first removes the forbidden matrix elements. The remaining finite symmetries are then applied in the usual way." -> "先考虑两个在转角 theta 下取得不同相位的局域态。只有 SO(2) 权相容时，它们之间的矩阵元才允许非零。在 initfromrep 中，这个连续 C_infinity 对称性写成单参数矩阵 D(theta)。MagneticTB 在 theta=0 处求导得到厄米生成元，先去掉连续对称性禁止的矩阵元，再照常加入其余有限对称约束。"
  ,"The continuous operation is internal: it acts on the local states but does not move atoms or create new bonds. The first two-state example makes the selection rule visible directly in the onsite Hamiltonian." -> "这个连续操作只作用于局域态，不移动原子，也不产生新键。下面先用最简单的二态例子，让这一选择定则直接显示在 onsite Hamiltonian 中。"
  ,"This page collects models that require more than the two-band graphene example. They include unequal orbital dimensions on two atoms, several Wyckoff orbits, spinful magnetic symmetry, induced representations, and an antiunitary operation that carries a state from one site to another. Each section gives the complete input and then shows the resulting Hamiltonian, orbital order or hopping, symmetry matrices, and bands." -> "这一页收集比两带石墨烯更复杂的模型，包括两个原子上轨道数不同、多个 Wyckoff 轨道、含自旋磁对称性、诱导表示，以及把一个格点上的态搬到另一个格点的反幺正操作。每一节都先给出完整输入，再显示 Hamiltonian、轨道或 hopping 的对应关系、对称矩阵和能带。"
  ,"sets up a model with the original MagneticTB implementation." -> "使用原版 MagneticTB 函数初始化模型。"
  ,"Primitive direct-lattice vectors used by the original version, stored by rows." -> "原版使用的原胞正格矢，按行排列。"
  ,"Numerical substitutions for the symbolic lattice parameters in latticeold." -> "为 latticeold 中的符号晶格参数给出数值。"
  ,"Fractional Wyckoff seeds and magnetic moments used by the original version." -> "原版使用的分数坐标 Wyckoff 代表点和磁矩。"
  ,"Complete ordered magnetic-group operation records used by the original version." -> "原版使用的完整、有序磁群操作。"
  ,"Local function basis, ordered exactly as wyckoffpositionold." -> "局域基函数，其顺序与 wyckoffpositionold 完全一致。"
  ,"Original Boolean diagnostic option, retained so that earlier notebooks still run." -> "原版的布尔诊断选项，为使以前的 notebook 继续运行而保留。"
  ,"Load MagneticTB 1.0" -> "载入 MagneticTB 1.0 版"
  ,StringJoin[
    "Start a fresh kernel and load MagneticTBOld` with the following command. Use separate ",
    "kernels for versions 1.0 and 2.0:"
  ] -> StringJoin[
    "先启动一个新内核，用下面的命令加载 MagneticTBOld`。1.0 版与 2.0 ",
    "版请分别在不同内核中使用："
  ]
  ,"Initialize the model with version 1.0" -> "用 1.0 版初始化模型"
  ,"Construct and display the Hamiltonian with version 1.0" -> "用 1.0 版求出并显示 Hamiltonian"
  ,"Compare results from versions 1.0 and 2.0" -> "比较 1.0 版与 2.0 版的结果"
  ,StringJoin[
    "To compare versions 1.0 and 2.0, enter the same lattice, positions, symmetry, and ",
    "orbitals in two separate kernels. Then compare the Hamiltonians and bands. Each version ",
    "needs its own initialization; neither uses the model set up in the other kernel."
  ] -> StringJoin[
    "如果要比较 1.0 版与 2.0 版，在两个独立内核中输入相同的晶格、位置、对称操作和轨道，再比较得到的",
    " Hamiltonian 和能带。两个版本都需要各自初始化，不会沿用另一个内核中的模型。"
  ]
  ,"Original version" -> "原版程序"
  ,"Run the original implementation in a separate fresh kernel. The tutorial Using the Original Version gives complete initold, symhamold, and bandManipulateold examples with their outputs." -> "原版程序要在单独的新内核中运行。《使用原版 MagneticTB》教程给出了完整的 initold、symhamold 和 bandManipulateold 例子及其输出。"
  ,"Here is a complete simple-cubic model in this input form. The output shows the lattice, generated sites, exact representation matrices, and first bond class:" -> "例如，下面用这种输入写一个完整的简单立方模型。结果依次列出晶格、生成的格点、精确表示矩阵和第一类键："
  ,"When another program writes the model input, it can first inspect the exact symmetry matrices returned here. The same physical input can then be passed to init or initfromrep for an interactive calculation:" -> "若输入是由其他程序产生的，可以先在这里检查精确的对称矩阵。确认无误后，再把同一组物理输入交给 init 或 initfromrep 继续计算："
  ,"For example, starting from one generator of C4, repeated multiplication produces its four elements. The returned list always begins with the identity; the trivial group is therefore {E}, not an empty list." -> "例如，从 C4 的一个生成元出发，反复相乘就得到四个群元。返回列表总以单位元开头，所以平凡群是 {E}，而不是空列表。"
  ,"When the complete finite group is already known, getGenerator looks for a smaller set whose products recover every group element." -> "已知完整有限群时，getGenerator 会寻找一组更小的生成元，使它们的乘积仍能给出全部群元。"
  ,"Supply the complete group, its identity, and its multiplication law. The function does not infer a product from the appearance of the elements. The name is retained so that existing notebooks continue to run." -> "输入必须包含完整群、单位元和乘法规则。函数不会只凭群元的写法猜测它们如何相乘。保留这个函数名是为了让已有 notebook 继续运行。"
  ,"To symmetrize a Wannier90 HR model, one needs the center of every Wannier orbital and the matrix of every symmetry operation in the same orbital order. This function constructs those data from the lattice, Wyckoff positions, symmetry operations, and local orbitals." -> "对称化 Wannier90 HR 模型时，需要知道每条 Wannier 轨道的中心，以及每个对称操作在同一轨道顺序下的矩阵。这个函数根据晶格、Wyckoff 位置、对称操作和局域轨道构造这些数据。"
  ,"The orbital order must agree with the external electronic-structure calculation. At present the function supports the documented VASP ordering." -> "轨道顺序必须与外部电子结构计算一致。目前该函数支持说明书中给出的 VASP 顺序。"
  ,"The result contains wcc for the fractional Wannier centers, DR for the full symmetry matrices, and symmetry for the corresponding ordered operations." -> "结果中，wcc 是分数坐标 Wannier 中心，DR 是完整对称矩阵，symmetry 则按相同顺序列出对应操作。"
  ,"The following complete call shows all six inputs. Together they fix the Wannier centers, orbital order, and symmetry matrices:" -> "下面的完整调用给出六项输入。它们共同确定 Wannier 中心、轨道顺序和对称矩阵："
  ,"Only the documented VASP orbital ordering is currently available:" -> "目前只支持说明书中给出的 VASP 轨道顺序："
  ,"debugQ accepts only True or False; any other value stops initialization:" -> "debugQ 只能取 True 或 False；其他取值会使初始化停止："
  ,"Retained for earlier calls and restricted to True or False. Its value does not change the supplied representation matrices or the resulting Hamiltonian." -> "为了兼容以前的调用而保留，只能取 True 或 False。它不改变用户给定的表示矩阵，也不改变最后的 Hamiltonian。"
  ,"Automatic displays a generating set, while All displays every operation. An integer or list of integers selects operations by their order in symminformation." -> "Automatic 显示一组生成元，All 显示全部操作。整数或整数列表则按 symminformation 中的顺序选择操作。"
  ,"Every matrix acts in the orbital order displayed by orbitalTable[]. Calling this function only displays the matrices already determined by init or initfromrep; it does not solve a Hamiltonian." -> "每个矩阵都按 orbitalTable[] 显示的轨道顺序作用。调用这个函数只会显示 init 或 initfromrep 已经确定的矩阵，不会另外求解 Hamiltonian。"
  ,"Each lattice translation R in a Wannier90 HR file is paired with a real-space hopping matrix H(R). Its row and column order is the Wannier-orbital order used when the file was written." -> "Wannier90 HR 文件中每个晶格平移 R 都对应一个实空间 hopping 矩阵 H(R)。它的行、列顺序就是写出该文件时采用的 Wannier 轨道顺序。"
  ,"Before returning these matrices, readHR checks the file header, degeneracies, translation vectors, orbital indices, and complex matrix elements. A malformed or incomplete file is rejected." -> "在返回这些矩阵之前，readHR 会检查文件头、简并度、平移向量、轨道编号和复矩阵元。格式错误或内容不完整的文件会被拒绝。"
  ,"Pass the returned Association to buildBlochHamiltonian to recover H(k). For repeated evaluation, define a function such as h[k_List] := buildBlochHamiltonian[hrData,k]." -> "把 readHR 返回的 Association 交给 buildBlochHamiltonian 即可恢复 H(k)。需要反复计算时，可以定义 h[k_List] := buildBlochHamiltonian[hrData,k]。"
  ,"Create a graphene HR file with hop, read it back, reconstruct H(k), and display the Hamiltonian at the Gamma point:" -> "用 hop 写出石墨烯 HR 文件，再用 readHR 读回、恢复 H(k)，并显示 Gamma 点的 Hamiltonian："
  ,"Use initold to run notebooks written for the original MagneticTB, or to compare the original and current constructions." -> "以前为原版 MagneticTB 编写的 notebook 可以继续使用 initold；比较原版与新版结果时也使用这个函数。"
  ,"Initialize graphene with the original version and display its onsite Hamiltonian. Run this example in a kernel that loads only MagneticTBOld`:" -> "用原版函数初始化石墨烯，并显示 onsite Hamiltonian。请在只载入 MagneticTBOld` 的内核中运行："
  ,"Add the onsite and nearest-neighbour terms to obtain the Bloch Hamiltonian of the original version:" -> "把 onsite 和最近邻项相加，得到原版的 Bloch Hamiltonian："
  ,"CartesianCoordinates -> True writes the nearest-neighbour term in Cartesian momentum coordinates. symmetrysetold belongs to the original version; here All uses every symmetry operation supplied to initold:" -> "CartesianCoordinates -> True 用笛卡尔动量坐标写出最近邻项。symmetrysetold 属于原版；这里用 All 选取交给 initold 的全部对称操作："
  ,"draws the bands of a Hamiltonian obtained with the original version and provides sliders for its parameters." -> "画出原版 Hamiltonian 的能带，并用滑块调节其中参数。"
  ,"Define the graphene path and open the interactive band panel of the original version:" -> "定义石墨烯的动量路径，然后打开原版的交互能带图："
  ," \[LongDash] initialize a model with version 1.0." -> " \[LongDash] 用 1.0 版初始化模型。"
  ," \[LongDash] construct a Hamiltonian with version 1.0." -> " \[LongDash] 用 1.0 版构造 Hamiltonian。"
  ," \[LongDash] explore bands obtained from version 1.0." -> " \[LongDash] 交互查看 1.0 版 Hamiltonian 的能带。"
  ,"Plot the bands with version 1.0" -> "用 1.0 版画能带"
  ,"After restarting Mathematica, search for MagneticTB in the Documentation Center, or open the guide directly with the following URI. The original MagneticTB functions are available from MagneticTBOld`; use a separate fresh kernel when comparing original and current calculations." -> "重启 Mathematica 后，可在帮助中心搜索 MagneticTB，也可用下面的 URI 直接打开指南。原版 MagneticTB 函数位于 MagneticTBOld` 中；比较原版和新版计算时，请分别使用两个新内核。"
  ,StringJoin[
    "MagneticTB 2.0 includes the functions from version 1.0. To continue using a 1.0 ",
    "notebook, load MagneticTBOld` and use the function names ending in old. The graphene ",
    "example below shows how to initialize a model, calculate its Hamiltonian, and plot the ",
    "bands with these functions."
  ] -> StringJoin[
    "MagneticTB 2.0 也提供了 1.0 版函数。要继续使用 1.0 版 ",
    "notebook，可以加载 MagneticTBOld`，调用名称末尾带 old ",
    "的函数。下面以石墨烯为例，介绍 1.0 版的初始化、Hamiltonian 和能带计算。"
  ]
  ,"Select an operation by its position in symminformation. Operation 25 is pure time reversal in this ordered gray group:" -> "按操作在 symminformation 中的位置进行选择。在这个有序灰群中，第 25 个操作是纯时间反演："
  ,"Select a short ordered list containing unitary operations and pure time reversal:" -> "选择一个同时包含幺正操作与纯时间反演的短有序列表："
  ,"Select one operation by its position in symminformation. Operation 25 is pure time reversal in this ordered gray group:" -> "按操作在 symminformation 中的位置选择一个操作。在这个有序灰群中，第 25 个操作是纯时间反演："
  ,"symhamold is kept for existing notebooks and for independent comparisons with the current version. Its calculation is entirely the one used by the original MagneticTB." -> "保留 symhamold 是为了继续运行已有 notebook，也便于与新版独立比较。它完全按照原版 MagneticTB 的方法计算。"
  ,"False uses the reciprocal-fractional momenta of the original version; True rewrites the result with its reciprocal lattice." -> "False 使用原版的倒空间分数动量；True 用原版倒格子改写结果。"
  ,"Initialize the complete original model used in this example. symhamold does not perform initialization:" -> "先完整初始化本例所用的原版模型。symhamold 不会执行初始化："
  ,"Initialize the original graphene model and construct its Hamiltonian before plotting:" -> "绘图前，先初始化原版石墨烯模型并构造 Hamiltonian："
  ,"Run the original MagneticTB functions in a separate fresh kernel. The tutorial Using the Original Version gives complete initold, symhamold, and bandManipulateold examples with their outputs." -> "原版 MagneticTB 函数要在单独的新内核中运行。《使用原版 MagneticTB》教程给出了完整的 initold、symhamold 和 bandManipulateold 例子及其输出。"

  ,"puts a collection of exact constraint matrices into one common cyclotomic field for exact row reduction and null-space calculations." -> "把一组精确约束矩阵编译到同一个分圆域中，以进行精确行化简与零空间计算。"
  ,"Use this advanced entry when a program already stores roots of unity in the documented root_sum form. For ordinary exact Wolfram Language matrices containing rationals, I, radicals, or rational-angle roots of unity, CyclotomicCommonNullSpace is the simpler entry." -> "当程序已经用文档规定的 root_sum 形式存储单位根时，可使用这个高级入口。对于含有有理数、I、根式或有理角单位根的普通精确 Wolfram Language 矩阵，使用 CyclotomicCommonNullSpace 更直接。"
  ,"input is an Association with \"coordinate_dimension\", \"constraints\", and optional \"max_cyclotomic_degree\" keys. Every constraint has the same number of columns. A constraint may be a rectangular nested list or an Association with \"rows\", \"columns\", and flat \"entries\"." -> "input 是包含 \"coordinate_dimension\"、\"constraints\" 以及可选 \"max_cyclotomic_degree\" 键的 Association。每个约束的列数必须相同。约束既可以是矩形嵌套列表，也可以是包含 \"rows\"、\"columns\" 和展平 \"entries\" 的 Association。"
  ,"A root_sum scalar has an exact rational constant and a list of terms with exact rational \"coefficient\", positive integer \"order\", and integer \"power\". The compiler takes the least common multiple of all root orders and stores every matrix in that one exact field." -> "root_sum 标量由一个精确有理常数和若干项组成；每项包含精确有理 \"coefficient\"、正整数 \"order\" 与整数 \"power\"。编译器取全部单位根阶数的最小公倍数，并把所有矩阵存入同一个精确域。"
  ,"The result contains exact internal matrix objects for CyclotomicRREF, CyclotomicNullSpace, and CyclotomicCommonKernel. Use RestoreCyclotomicExpression before presenting their entries as ordinary Wolfram Language expressions." -> "结果包含供 CyclotomicRREF、CyclotomicNullSpace 和 CyclotomicCommonKernel 使用的精确内部矩阵对象。若要把其中元素显示为普通 Wolfram Language 表达式，请先调用 RestoreCyclotomicExpression。"
  ,"Compile two three-coordinate constraints containing the fourth root of unity i. The summary shows that both blocks share conductor 4:" -> "编译两个含四次单位根 i 的三坐标约束。摘要表明两个矩阵块共用导数 4："
  ,"Restore the two compiled blocks to see the exact matrices represented by the shared field:" -> "还原两个已编译矩阵块，查看这个公共域所表示的精确矩阵："
  ,"Every matrix must have the declared coordinate dimension. A row with only two entries in a three-coordinate problem fails explicitly:" -> "每个矩阵都必须具有声明的坐标维数。在三坐标问题中，只有两个元素的行会明确失败："
  ,"The field degree may not exceed max_cyclotomic_degree. This boundary stops the calculation instead of changing to approximate arithmetic:" -> "域次数不得超过 max_cyclotomic_degree。超过上限时计算会停止，而不会改用近似算术："

  ,"computes a deterministic exact reduced row echelon form for one matrix returned by CompileCyclotomicMatrices." -> "对 CompileCyclotomicMatrices 返回的一个矩阵计算确定性的精确行最简阶梯形。"
  ,"matrix must be an exact matrix object from CompileCyclotomicMatrices. The function returns \"reduced_matrix\", \"pivot_columns\", \"free_columns\", \"rank\", and \"nullity\"." -> "matrix 必须是 CompileCyclotomicMatrices 产生的精确矩阵对象。函数返回 \"reduced_matrix\"、\"pivot_columns\"、\"free_columns\"、\"rank\" 和 \"nullity\"。"
  ,"pivot_columns and free_columns use 0-based indices so that the result has the same deterministic convention as the cross-language fixture format. Add one when using an index with Wolfram Language Part." -> "pivot_columns 和 free_columns 使用从 0 开始的索引，以便结果与跨语言 fixture 格式采用同一确定性约定。把索引用于 Wolfram Language 的 Part 时须加 1。"
  ,"The production calculation uses MagneticTB's exact cyclotomic row-reduction algorithm. It does not call the built-in RowReduce and does not switch to floating point, SVD, or another fallback." -> "生产计算使用 MagneticTB 自身的分圆域精确行化简算法。它不调用内置 RowReduce，也不会切换到浮点、SVD 或其他 fallback。"
  ,"Compile a rank-one rational matrix and display its exact RREF together with its 0-based pivot and free columns:" -> "编译一个秩为 1 的有理矩阵，并显示其精确 RREF 以及从 0 开始的主元列和自由列："
  ,"The same call works in a nontrivial cyclotomic field. Here the matrix entry i is kept exact throughout the reduction:" -> "同一调用也适用于非平凡分圆域。这里的矩阵元素 i 在整个化简过程中保持精确："
  ,"An ordinary nested list is not a compiled exact matrix object. Use CompileCyclotomicMatrices first, or use CyclotomicCommonNullSpace for ordinary exact matrices:" -> "普通嵌套列表不是已编译的精确矩阵对象。请先用 CompileCyclotomicMatrices 编译；对于普通精确矩阵，也可直接使用 CyclotomicCommonNullSpace："

  ,"finds a deterministic exact null-space basis for one matrix returned by CompileCyclotomicMatrices." -> "为 CompileCyclotomicMatrices 返回的一个矩阵求确定性的精确零空间基。"
  ,"The result contains both \"nullspace_rows\" and \"basis_matrix\". nullspace_rows has one basis vector per row; basis_matrix is its transpose, so matrix . basis_matrix is exactly zero." -> "结果同时包含 \"nullspace_rows\" 和 \"basis_matrix\"。nullspace_rows 每行存放一个基向量；basis_matrix 是它的转置，因此 matrix . basis_matrix 精确为零。"
  ,"The result also gives the exact reduced matrix, 0-based pivot_columns and free_columns, rank, nullity, and exact_residual_verified. The row basis is ordered by increasing free-column index." -> "结果还给出精确化简矩阵、从 0 开始的 pivot_columns 与 free_columns、rank、nullity 以及 exact_residual_verified。行基按自由列索引递增排列。"
  ,"The production algorithm does not call the built-in NullSpace or RowReduce. A failure is returned directly; the function never changes to floating point, SVD, or another null-space algorithm." -> "生产算法不调用内置 NullSpace 或 RowReduce。发生失败时直接返回 Failure；函数绝不会改用浮点、SVD 或其他零空间算法。"
  ,"For this rank-one matrix, the two row-basis vectors are {-2,1,0} and {-3,0,1}. Their transpose is the column-basis matrix annihilated by the input:" -> "对于这个秩为 1 的矩阵，两个行基向量为 {-2,1,0} 和 {-3,0,1}。它们的转置就是被输入矩阵湮灭的列基矩阵："
  ,"Full column rank gives a zero-dimensional null space. The basis matrix retains the exact 2 by 0 shape even though it has no entries:" -> "满列秩对应零维零空间。尽管没有元素，基矩阵仍保留精确的 2×0 形状："
  ,"CyclotomicNullSpace accepts a compiled exact matrix, not a raw nested list. The direct ordinary-matrix entry is CyclotomicCommonNullSpace:" -> "CyclotomicNullSpace 接受已编译的精确矩阵，而不是原始嵌套列表。普通矩阵的直接入口是 CyclotomicCommonNullSpace："

  ,"finds the vectors annihilated by every exact constraint block in a compiled cyclotomic problem." -> "求出被已编译分圆域问题中每个精确约束块共同湮灭的向量。"
  ,"compiled must be the complete Association returned by CompileCyclotomicMatrices. All blocks have the same coordinate dimension and exact field context." -> "compiled 必须是 CompileCyclotomicMatrices 返回的完整 Association。所有矩阵块具有相同坐标维数和精确域上下文。"
  ,"The result stores the common basis by columns in \"basis_matrix\" and by rows in \"nullspace_rows\". iteration_nullities begins with the coordinate dimension and records the remaining nullity after each input block." -> "结果在 \"basis_matrix\" 中按列存储公共基，在 \"nullspace_rows\" 中按行存储公共基。iteration_nullities 以坐标维数开始，随后记录处理每个输入块后的剩余零度。"
  ,"The method is the production iterative common-kernel algorithm. Every block is checked against the final exact basis. If any step or residual check fails, the function returns Failure; it never stacks all blocks as an automatic fallback and never uses floating point or SVD." -> "该方法是生产用迭代公共核算法。最终精确基会逐块接受检验。若任一步骤或残差检验失败，函数返回 Failure；它不会自动堆叠所有矩阵块作为 fallback，也不会使用浮点或 SVD。"
  ,"The first constraint removes the first coordinate and the second removes the second. Their common kernel is therefore spanned by {0,0,1}:" -> "第一个约束去掉第一坐标，第二个约束去掉第二坐标，因此它们的公共核由 {0,0,1} 张成："
  ,"With no constraint blocks, every coordinate is free. The row basis is the identity and iteration_nullities contains only the initial dimension:" -> "没有约束块时，每个坐标都是自由的。行基为单位矩阵，iteration_nullities 只含初始维数："
  ,"A problem with coordinate dimension zero has rank and nullity zero and retains a 0 by 0 basis matrix:" -> "坐标维数为零的问题，其秩与零度均为零，并保留一个 0×0 基矩阵："
  ,"The function accepts only a complete compiled problem. It does not infer missing context or dimensions from an arbitrary Association:" -> "函数只接受完整的已编译问题，不会从任意 Association 猜测缺失的上下文或维数："

  ,"accepts ordinary exact Wolfram Language matrices and returns a row basis for their common null space." -> "接受普通精确 Wolfram Language 矩阵，并返回其公共零空间的行基。"
  ,"This is the main entry for ordinary exact matrices. If rows is the returned list, Transpose[rows] is the column-basis matrix and every input block Ai satisfies Ai . Transpose[rows] == 0 exactly." -> "这是普通精确矩阵的主入口。若返回列表为 rows，则 Transpose[rows] 是列基矩阵，并且每个输入块 Ai 都精确满足 Ai . Transpose[rows] == 0。"
  ,"Supported scalar inputs include exact integers and rationals, I, square roots of exact rationals, roots of unity written with rational-angle Exp, Sin, or Cos, Conjugate, and exact sums, products, and supported powers of these expressions. Machine numbers, symbolic parameters, and unsupported algebraic numbers such as general cube roots return Failure." -> "支持的标量输入包括精确整数和有理数、I、精确有理数的平方根、用有理角 Exp、Sin 或 Cos 写出的单位根、Conjugate，以及这些表达式的精确和、积与受支持幂。机器数、符号参数以及一般立方根等不受支持的代数数会返回 Failure。"
  ,"All matrices must have the same column count. Empty constraints require CoordinateDimension because no matrix is available from which to infer that count." -> "所有矩阵必须具有相同列数。空约束必须指定 CoordinateDimension，因为此时没有矩阵可供推断列数。"
  ,"The function compiles one exact field and calls the production iterative common-kernel algorithm. It does not use the built-in NullSpace or RowReduce and never changes to floating point, SVD, stacked constraints, or another fallback." -> "函数编译一个公共精确域，并调用生产用迭代公共核算法。它不使用内置 NullSpace 或 RowReduce，也绝不会改用浮点、SVD、堆叠约束或其他 fallback。"
  ,"Number of vector coordinates. Automatic infers the common column count from nonempty input; an explicit nonnegative integer is required when the constraint list is empty." -> "向量坐标数。Automatic 从非空输入推断公共列数；约束列表为空时必须显式给出非负整数。"
  ,"Positive integer upper bound on the exact cyclotomic field degree. Exceeding the bound returns Failure instead of using approximate arithmetic." -> "精确分圆域次数的正整数上限。超过上限时返回 Failure，而不使用近似算术。"
  ,"Find the common null space of two ordinary exact constraint blocks containing three independent square roots. The two displayed residual matrices are exactly zero:" -> "求两个普通精确约束块的公共零空间，其中含三个相互独立的平方根。显示的两个残差矩阵都精确为零："
  ,"CoordinateDimension fixes the vector space when the constraint list is empty. With three coordinates, the complete row basis is the identity:" -> "约束列表为空时，CoordinateDimension 用来固定向量空间。三坐标空间的完整行基是单位矩阵："
  ,"The degree bound applies to the exact field degree, not to a numerical tolerance. The fourth-root field has degree two, so limit two succeeds and limit one fails:" -> "次数上限约束的是精确域次数，而不是数值容差。四次单位根域的次数为 2，因此上限 2 成功，上限 1 失败："
  ,"For a zero-dimensional vector space, empty constraints return the unique empty row basis:" -> "对于零维向量空间，空约束返回唯一的空行基："
  ,"Without a nonempty matrix or an explicit CoordinateDimension, the vector-space dimension is undefined:" -> "既没有非空矩阵也没有显式 CoordinateDimension 时，向量空间维数未定义："
  ,"All blocks must have the same number of columns:" -> "所有矩阵块必须具有相同列数："
  ,"Machine numbers and unsupported exact algebraic expressions fail for different reasons; neither is approximated or sent to a numerical fallback:" -> "机器数和不受支持的精确代数表达式会因不同原因失败；二者都不会被近似化，也不会转入数值 fallback："

  ,"turns an exact matrix, null-space result, or common-kernel result into ordinary exact Wolfram Language expressions." -> "把精确矩阵、零空间结果或公共核结果转换为普通精确 Wolfram Language 表达式。"
  ,"restores one exact cyclotomic element using its matching field context." -> "使用匹配的域上下文还原一个精确分圆域元素。"
  ,"The one-argument form accepts an exact_matrix, kernel_result, or common_kernel_result produced by the public cyclotomic functions. A restored matrix records \"rows\", \"columns\", and row-major flat \"entries\"." -> "单参数形式接受公共分圆域函数产生的 exact_matrix、kernel_result 或 common_kernel_result。还原后的矩阵记录 \"rows\"、\"columns\" 以及按行优先展平的 \"entries\"。"
  ,"For an RREF result, restore its \"reduced_matrix\" field. The complete rref_result Association is metadata and is not itself a restorable object." -> "对于 RREF 结果，应还原其中的 \"reduced_matrix\" 字段。完整 rref_result Association 属于元数据，其本身不是可还原对象。"
  ,"The two-argument form requires the exact context to which element belongs. A context mismatch or malformed object returns Failure; no numerical approximation is used." -> "双参数形式要求提供 element 所属的精确上下文。上下文不匹配或对象格式错误时返回 Failure；不会采用数值近似。"
  ,"Compile a matrix containing the fourth root of unity and restore it as the ordinary exact expression I:" -> "编译一个含四次单位根的矩阵，并把它还原为普通精确表达式 I："
  ,"Restore a complete kernel result, then display its row-basis and column-basis matrices:" -> "还原完整核结果，然后显示其中的行基与列基矩阵："
  ,"The two-argument form restores one element from its matching context. The second entry of this matrix is exactly I:" -> "双参数形式从匹配上下文中还原一个元素。这个矩阵的第二个元素精确为 I："
  ,"An arbitrary Association is not an exact cyclotomic object and fails explicitly:" -> "任意 Association 不是精确分圆域对象，因此会明确失败："

  ,"Exact cyclotomic linear algebra" -> "分圆域精确线性代数"
  ," \[LongDash] compile exact constraint blocks into one cyclotomic field." -> " \[LongDash] 把精确约束块编译到同一个分圆域。"
  ," \[LongDash] compute a deterministic exact reduced row echelon form." -> " \[LongDash] 计算确定性的精确行最简阶梯形。"
  ," \[LongDash] compute an exact null-space row and column basis for one compiled matrix." -> " \[LongDash] 计算一个已编译矩阵的精确零空间行基与列基。"
  ," \[LongDash] compute the iterative exact common kernel of compiled constraint blocks." -> " \[LongDash] 计算已编译约束块的迭代精确公共核。"
  ," \[LongDash] find a common null-space row basis directly from ordinary exact matrices." -> " \[LongDash] 直接从普通精确矩阵求公共零空间行基。"
  ," \[LongDash] restore compiled cyclotomic objects as exact Wolfram Language expressions." -> " \[LongDash] 把已编译分圆域对象还原为精确 Wolfram Language 表达式。"

  ,"RepresentationMode -> \"Induced\" accepts a reference-site stabilizer representation. For the identity group, the induced and direct constructions give the same exact onsite Hamiltonian:" -> "RepresentationMode -> \"Induced\" 接受参考点稳定子表示。对于单位群，诱导构造与直积构造给出相同的精确 onsite Hamiltonian："
  ,StringJoin[
    "A chiral P3 model makes the Hermitian option visible. With Hermitian -> True the reverse ",
    "hoppings are fixed by t1 and t2; Hermitian -> False introduces independent reverse ",
    "amplitudes t3 and t4:"
  ] -> "手性 P3 模型可以清楚显示 Hermitian 选项的作用：Hermitian -> True 时反向 hopping 由 t1、t2 固定；Hermitian -> False 时新增独立的反向振幅 t3、t4："
  ,"KernelMethod -> \"Stacked\" solves the same exact symmetry constraints as one joined matrix. It gives the same nearest-neighbour Hamiltonian as the default iterative method:" -> "KernelMethod -> \"Stacked\" 把相同的精确对称约束合并为一个矩阵求解，得到的最近邻 Hamiltonian 与默认迭代方法相同："
  ,"ValidationLevel -> \"Full\" additionally checks independence and rank-nullity after solving. These checks leave the exact nearest-neighbour Hamiltonian unchanged:" -> "ValidationLevel -> \"Full\" 在求解后进一步检查基的线性独立性和秩-零度关系；这些检查不改变精确最近邻 Hamiltonian："
  ,"By default init finds the first 10 complete bond shells. If a model needs a farther shell, choose a larger InitialBondShells value and run init again; symham itself does not search for new bonds." -> "init 默认查找前 10 个完整近邻 hopping 层。若模型需要更远的 hopping，应设置更大的 InitialBondShells 并重新运行 init；symham 不会自行继续搜索。"
  ,"InitialBondShells, the testing-line default is 10 complete shells, matching init and initfromrep. An explicit value in the input Association always takes precedence." -> "InitialBondShells 在当前 testing 版本线中默认为 10 个完整近邻 hopping 层，与 init 和 initfromrep 一致。输入 Association 中的显式设置始终优先。"
  ,"InitialBondShells is 10 by default; the graphene calculation above deliberately requested only 3. Asking for shell 4 must therefore stop and request a new init with a larger value. MagneticTB does not search farther or change the null-space method without being asked. If a new initialization fails, no earlier model is used in its place." -> "InitialBondShells 默认为 10；上面的石墨烯计算特意只请求了 3 个近邻 hopping 层。因此请求第 4 层时必须停止，并要求用更大值重新 init。MagneticTB 不会未经请求就扩大搜索范围或更换零空间方法。若新初始化失败，也不会回退使用旧模型。"

  ,"Real-space, surface, and topological calculations" -> "实空间、表面与拓扑计算"
  ,"Optional MSGCorep interface" -> "可选 MSGCorep 接口"
  ," — draw a static band plot with an automatic conventional path." -> " — 使用自动约定路径绘制静态能带图。"
  ," — select a built-in conventional path from the current lattice metric." -> " — 从当前晶格度量选取内置约定路径。"
  ," — draw atoms, primitive cells, and magnetic-moment arrows." -> " — 绘制原子、原胞和磁矩箭头。"
  ," — draw the first Brillouin zone and selected k path." -> " — 绘制第一 Brillouin 区及选定的 k 路径。"
  ," — plot Wilson-loop phase or Wannier-center branches along a parameter path." -> " — 绘制 Wilson loop 相位或沿参数路径变化的 Wannier 中心分支。"
  ," — plot semi-infinite surface spectral weight on a momentum-energy grid." -> " — 在动量-能量网格上绘制半无限表面谱权重。"
  ," — draw signed Berry curvature on an oriented momentum slice." -> " — 在定向动量切片上绘制带符号的 Berry 曲率。"
  ," — draw the {F23,F31,F12} Berry-curvature vector field." -> " — 绘制 {F23,F31,F12} Berry 曲率矢量场。"
  ," — map finite-system probability and phase onto matrix-aligned real-space basis records." -> " — 根据与矩阵顺序一致的实空间基记录绘制有限体系的概率与相位。"
  ," — read numerical real-space hopping blocks from solved shells." -> " — 从已求解的近邻 hopping 层读取数值实空间跃迁块。"
  ," — transform hopping blocks to an exact integer cell." -> " — 把跃迁块变换到精确整数晶胞。"
  ," — sum hopping blocks into a periodic Bloch Hamiltonian." -> " — 将跃迁块求和为周期 Bloch Hamiltonian。"
  ," — assemble a finite sparse real-space Hamiltonian." -> " — 组装有限稀疏实空间 Hamiltonian。"
  ," — assemble a slab, ribbon, or rod Hamiltonian." -> " — 组装薄膜、条带或棒状 Hamiltonian。"
  ," — compute a semi-infinite principal-layer surface Green function." -> " — 计算半无限主层表面 Green 函数。"
  ," — evaluate the semi-infinite surface spectral weight." -> " — 计算半无限表面谱权重。"
  ," — compute Wilson-loop phases with orbital-center sewing." -> " — 用轨道中心缝合计算 Wilson loop 相位。"
  ," — compute the Berry phase of an explicitly sampled closed path." -> " — 计算显式采样闭合路径的 Berry 相位。"
  ," — compute Berry curvature from an oriented Wilson plaquette." -> " — 从定向 Wilson 小方格计算 Berry 曲率。"
  ," — compute the integer charge enclosed around a gapless point." -> " — 计算无隙点周围包围的整数拓扑荷。"
  ," — run a bounded periodic grid-and-refinement search for direct-gap zeros." -> " — 对直接能隙零点执行有界周期网格与精化搜索。"
  ," — convert the ordered operations of an installed MSGCorep package." -> " — 转换已安装 MSGCorep 包产生的有序对称操作。"
  ," — print single- or double-valued band corepresentations at supplied k points." -> " — 在给定 k 点输出单值或双值能带共表示。"
  ,"Crystal structure, Brillouin zone, and k paths" -> "晶体结构、Brillouin 区与 k 路径"
  ,"Real-space Hamiltonians, surfaces, and topology" -> "实空间 Hamiltonian、表面与拓扑"
  ,"Band corepresentations with MSGCorep" -> "使用 MSGCorep 计算能带共表示"
  ,"When the program-generated input omits InitialBondShells, the testing-line default is 10 complete shells, matching init and initfromrep. An explicit value in the input Association always takes precedence." -> "当程序生成的输入省略 InitialBondShells 时，当前 testing 版本线默认准备 10 个完整近邻 hopping 层，与 init 和 initfromrep 一致。输入 Association 中的显式值始终优先。"
  ,"This complete call supplies every public initfromrep option explicitly. The exact onsite Hamiltonian is computed from the supplied ordered matrices:" -> "这个完整调用显式给出 initfromrep 的每个公开选项。精确 onsite Hamiltonian 由所提供的有序矩阵计算得到："
  ,"This complete call supplies every public init option explicitly. The displayed Grid confirms the basis prepared from those choices:" -> "这个完整调用显式给出 init 的每个公开选项。显示的 Grid 确认了由这些设置准备的基："
  ,"Start here: install, update, and first models" -> "从这里开始：安装、更新与首个模型"
  ,"Crystal, symmetry, and group theory" -> "晶体、对称性与群论"
  ,"Hamiltonians and band calculations" -> "Hamiltonian 与能带计算"
  ,"Plotting and reciprocal-space paths" -> "绘图与倒空间路径"
  ,"Physical properties and topology" -> "物性与拓扑"
  ,"Exact mathematics and cyclotomic fields" -> "精确数学与分圆域"
  ," — install, update, uninstall, load, and build a first model." -> " — 安装、更新、卸载并加载 MagneticTB，然后构造首个模型。"
  ," — follow complete workflows across representative models." -> " — 查看典型模型的完整工作流程。"
  ," — open independently validated model notebooks." -> " — 打开经过独立验证的模型 notebook。"
  ," — run MagneticTB 1.0 in a separate fresh kernel." -> " — 在单独的新内核中使用 MagneticTB 1.0 版。"
  ," — work through magnetic Wyckoff models across all seven crystal systems." -> " — 学习覆盖七大晶系的磁性 Wyckoff 模型。"
  ," — construct exact continuous-symmetry representations." -> " — 构造精确的连续对称性表示。"
  ," — use verified single- and double-valued MSGCorep workflows." -> " — 使用经过验证的 MSGCorep 单值与双值工作流程。"
  ," — connect crystal, magnetic-moment, k-path, Brillouin-zone, and band plots." -> " — 串联晶体、磁矩、k 路径、Brillouin 区和能带绘图。"
  ," — build real-space, surface, and topological workflows." -> " — 构建实空间、表面与拓扑工作流程。"
  ,StringJoin[
    "Expand the Hamiltonian about Gamma through total order three, then pass it to ",
    "pointChernNumber. The function calculates the occupied-band Berry flux through a cube ",
    "surrounding the point. The occupied band must remain separated from the other band on ",
    "the cube surface; the enclosed Chern number is given below."
  ] -> StringJoin[
    "先把 Hamiltonian 在 Gamma 点附近展开到总阶数三阶，再用 ",
    "pointChernNumber 计算包围该点的立方体表面上的占据带 Berry ",
    "通量。立方体表面上，占据带必须与另一条带分离。下面的结果给出立方体内节点的 Chern 数。"
  ]
  ,StringJoin[
    "Run the following command to display the parameter sliders. Move them to see how the ",
    "bands change:"
  ] -> "运行下面的命令，可以用滑块调整参数，观察能带如何变化："
  ,"Use the following command to adjust the parameters and view the bands:" -> "用下面的命令调整参数并查看能带："
  ,"Run the last line to display the sliders and plot the bands:" -> "运行最后一行，就可以通过滑块调整参数并画出能带："

  ,"plotRange -> {-1,2.2} restricts only the visible vertical range of the same tight-binding/reference overlay:" -> "plotRange -> {-1,2.2} 只限制同一紧束缚能带与参考能带叠加图的可见纵轴范围："

  ,"\"Hermitian\" -> False exports the shell previously solved with independent directed hoppings. The same setting is supplied to symham and hop:" -> "\"Hermitian\" -> False 导出先前按独立有向跃迁求解的近邻 hopping 层；symham 与 hop 必须使用相同设置："
  ,"\"KernelMethod\" -> \"Stacked\" exports the shell solved by the stacked exact null-space method. The selector must match in symham and hop:" -> "\"KernelMethod\" -> \"Stacked\" 导出由堆叠式精确零空间方法求解的近邻 hopping 层；symham 与 hop 中的选择必须一致："
  ,"\"ValidationLevel\" -> \"Full\" exports the shell solved with full exact validation. The selector must match in symham and hop:" -> "\"ValidationLevel\" -> \"Full\" 导出经过完整精确验证的近邻 hopping 层；symham 与 hop 中的选择必须一致："
  ,"\"hrExport\" -> file writes the generated Wannier90 text to that explicit path instead of returning it in the notebook:" -> "\"hrExport\" -> file 把生成的 Wannier90 文本写入指定路径，而不是在 notebook 中返回文本："
  ,"\"hrExport\" -> directory writes wannier90_hr.dat there instead of printing the text in the notebook:" -> "\"hrExport\" -> directory 把 wannier90_hr.dat 写入该目录，而不在 notebook 中打印文本："
  ,"\"hrExport\" writes the HR data to the specified file. With None, hop returns the text directly in the notebook:" -> "\"hrExport\" 把 HR 数据写入指定文件；设为 None 时，hop 直接在 notebook 中返回文本："
  ,"\"RealDigits\" -> 8 writes eight digits after the decimal point for both real and imaginary hopping components:" -> "\"RealDigits\" -> 8 使跃迁的实部和虚部都在小数点后写出八位："
  ,"\"wcc\" supplies one explicit fractional center per orbital when converting a Bloch Hamiltonian. This graphene Hamiltonian and both centers come from the model initialized above:" -> "转换 Bloch Hamiltonian 时，\"wcc\" 为每条轨道显式提供一个分数坐标中心。这里的石墨烯 Hamiltonian 和两个中心都来自上面初始化的模型："
  ,"When converting a supplied Hamiltonian, set \"wcc\" to the fractional center of each orbital. This example uses the graphene centers; for a periodic-gauge Hamiltonian, replace them with one zero vector per orbital:" -> "转换用户给出的 Hamiltonian 时，必须用 \"wcc\" 指定每条轨道的分数坐标中心。这里使用石墨烯的两个轨道中心；若 Hamiltonian 已写成周期规范，则改为每条轨道一个零向量："
  ,"\"TranslationTolerance\" -> 10^-6 accepts a 10^-8 numerical offset in the symmetry-generated tetragonal C4 chain phases and rounds them to the two nearest-neighbour cells:" -> "\"TranslationTolerance\" -> 10^-6 接受对称性生成的四方 C4 链相位中 10^-8 的数值偏差，并把它们舍入到两个最近邻晶胞："

  ,"\"BoundaryConditions\" -> {\"Open\",\"Open\",\"Periodic\"} wraps the chain periodically. MatrixPlot displays the actual finite Hamiltonian, including the corner couplings created by the periodic boundary:" -> "\"BoundaryConditions\" -> {\"Open\",\"Open\",\"Periodic\"} 把链首尾周期连接。MatrixPlot 显示实际有限 Hamiltonian，其中包括周期边界生成的首尾耦合："
  ,"\"HermiticityTolerance\" -> 10^-12 applies a stricter residual threshold to the same exact hopping data. MatrixPlot shows the accepted finite Hamiltonian:" -> "\"HermiticityTolerance\" -> 10^-12 对同一组精确 hopping 数据施加更严格的残差阈值。MatrixPlot 显示通过检查的有限 Hamiltonian："
  ,"\"Output\" -> \"Data\" returns the Hamiltonian together with its geometry and basis metadata. Extract the public Hamiltonian field and plot the matrix itself:" -> "\"Output\" -> \"Data\" 同时返回 Hamiltonian、几何信息和基底元数据；取出公开的 Hamiltonian 字段并直接绘制矩阵："

  ,"lattice -> IdentityMatrix[3] supplies the direct-lattice row vectors for the explicit-representation route. The evaluated lattice stored in the model is the identity matrix:" -> "lattice -> IdentityMatrix[3] 为显式表示入口提供按行排列的正格矢；模型中保存的求值后晶格为单位矩阵："
  ,"lattpar -> {a->1} evaluates the symbolic cubic lattice before the supplied matrices are installed. The session retains the exact replacement rule:" -> "lattpar -> {a->1} 在安装输入矩阵前对符号立方晶格求值；会话保留这条精确替换规则："
  ,"wyckoffposition -> {{{0,0,0},{0,0,0}}} supplies one nonmagnetic seed at the origin. With the identity group, its compiled site orbit contains one site:" -> "wyckoffposition -> {{{0,0,0},{0,0,0}}} 在原点提供一个非磁种子；对单位群，编译后的格点轨道只含一个格点："
  ,"symminformation -> {E} supplies the complete ordered finite group. The stored symmetry label is in the same order as the representation matrices:" -> "symminformation -> {E} 提供完整有序有限群；保存的对称标签与表示矩阵采用相同顺序："
  ,"repinformation -> {{IdentityMatrix[2]}} supplies one exact matrix for the one group operation. The installed full representation contains that same matrix:" -> "repinformation -> {{IdentityMatrix[2]}} 为唯一群操作提供一个精确矩阵；安装后的完整表示包含同一个矩阵："
  ,"orbitalLabels -> {{\"alpha\",\"beta\"}} names the two coordinates of the supplied representation. orbitalTable displays those labels in the fixed Hamiltonian order:" -> "orbitalLabels -> {{\"alpha\",\"beta\"}} 为输入表示的两个坐标命名；orbitalTable 按固定 Hamiltonian 顺序显示这些标签："
  ,"debugQ -> True uses the supported Boolean compatibility setting. The session records it without altering the supplied representation matrices:" -> "debugQ -> True 使用受支持的布尔兼容设置；会话记录该值，但不改变输入的表示矩阵："
  ,"InitialBondShells -> 1 prepares exactly the onsite shell for this explicit model. The session exposes that compiled shell count:" -> "InitialBondShells -> 1 为这个显式模型只准备 onsite 层；会话直接给出已编译层数："
  ,"RepresentationMode -> \"DirectProduct\" installs the supplied matrices as the full local action. The session records the selected construction mode:" -> "RepresentationMode -> \"DirectProduct\" 把输入矩阵作为完整局域作用安装；会话记录所选构造模式："
  ,"lattpar -> {a->2,c->4} evaluates the default symbolic lattice before the supplied matrices are installed. The stored numerical lattice shows the substitutions:" -> "lattpar -> {a->2,c->4} 在安装输入矩阵前对默认符号晶格求值；保存的数值晶格显示了这些替换："
  ,"orbitalLabels -> {{\"alpha\",\"beta\"}} names the two coordinates of the supplied representation. The session returns the labels in that same order:" -> "orbitalLabels -> {{\"alpha\",\"beta\"}} 为输入表示的两个坐标命名；会话按相同顺序返回这些标签："
  ,"RepresentationMode -> \"Induced\" accepts a reference-site stabilizer representation. For the identity group, only one 1 by 1 local matrix is needed:" -> "RepresentationMode -> \"Induced\" 接受参考点稳定子表示；对于单位群，只需一个 1×1 局域矩阵："

  ,"lattice -> {{a,0,0},{-a/2,Sqrt[3]a/2,0},{0,0,c}} gives the three primitive vectors of graphene by rows. The other four model inputs are supplied in the same call. The second bond shell then displays the three nearest-neighbour carbon bonds fixed by this lattice and Wyckoff position:" -> "lattice -> {{a,0,0},{-a/2,Sqrt[3]a/2,0},{0,0,c}} 按行给出石墨烯的三个原胞基矢。同一次调用还完整提供其余四项物理输入；第 2 键层显示由该晶格和 Wyckoff 位置确定的三条最近邻碳键："
  ,"lattpar -> {a->1,c->3} assigns numerical lengths to every symbol in lattice before bonds are found. In these units showbonds[2] reports the graphene nearest-neighbour distance and the three translated target positions:" -> "lattpar -> {a->1,c->3} 在搜索键之前为 lattice 中的每个符号指定数值长度。在这组单位下，showbonds[2] 显示石墨烯最近邻距离和三个平移后的目标位置："
  ,"wyckoffposition -> {{{1/3,2/3,0},{0,0,0}}} supplies one nonmagnetic carbon seed in fractional coordinates. The gray group generates the second carbon site; orbitalTable[] shows both sites in the Hamiltonian order:" -> "wyckoffposition -> {{{1/3,2/3,0},{0,0,0}}} 以分数坐标提供一个非磁碳原子种子。灰群生成第二个碳格点；orbitalTable[] 按 Hamiltonian 顺序显示这两个格点："
  ,"symminformation -> msgop[gray[191]] supplies the complete ordered magnetic-space-group records used for graphene. showSymmetryRepresentations[] displays the actual orbital matrices of the selected generators, including their unitary or antiunitary type:" -> "symminformation -> msgop[gray[191]] 提供石墨烯所用的完整有序磁空间群记录。showSymmetryRepresentations[] 显示所选生成元的实际轨道矩阵及其幺正或反幺正类型："
  ,"basisFunctions -> {{\"pz\"}} assigns one local pz orbital to the inequivalent carbon seed. Symmetry transports it to the second site, so orbitalTable[] contains two pz basis states and fixes their Hamiltonian row order:" -> "basisFunctions -> {{\"pz\"}} 为不等价碳原子种子指定一个局域 pz 轨道。对称性把它搬运到第二个格点，因此 orbitalTable[] 显示两个 pz 基态并固定其 Hamiltonian 行顺序："
  ,"debugQ -> True enables the compatibility diagnostics while the same complete graphene model is built. It does not change the mathematics: the onsite Hamiltonian remains the two-site diagonal matrix shown here:" -> "debugQ -> True 在构造同一个完整石墨烯模型时启用兼容诊断信息，但不改变数学结果；onsite Hamiltonian 仍是这里显示的双格点对角矩阵："
  ,"InitialBondShells -> 3 prepares onsite, nearest-neighbour, and next-nearest-neighbour shells during the complete graphene initialization. showbonds[3] displays the six next-nearest neighbours of each carbon site:" -> "InitialBondShells -> 3 在完整石墨烯初始化中准备 onsite、最近邻和次近邻三个键层。showbonds[3] 显示每个碳格点的六个次近邻："
  ,"GenerateSymmetryGroup -> True closes the supplied C4z generator to the complete fourfold group before constructing a px/py square-lattice model. The displayed matrices show the four generated operations acting on the two local orbitals:" -> "GenerateSymmetryGroup -> True 先把输入的 C4z 生成元闭合成完整四重群，再构造 px/py 方格模型。所显示的矩阵给出四个生成操作对两个局域轨道的作用："

  ,"lattice supplies the three tetragonal primitive vectors by rows. The same complete initfromrep call also supplies its lattice parameters, Wyckoff seed, ordered C4 group, exact matrices, and orbital labels. showbonds[2] displays the four in-plane nearest neighbours:" -> "lattice 按行给出三个四方原胞基矢。同一次完整 initfromrep 调用还给出晶格参数、Wyckoff 种子、有序 C4 群、精确矩阵和轨道标签；showbonds[2] 显示四个面内最近邻："
  ,"lattpar -> {a->1,c->2} assigns the numerical tetragonal lattice constants before the exact representation is installed. showbonds[2] then reports unit in-plane nearest-neighbour distance:" -> "lattpar -> {a->1,c->2} 在安装精确表示前指定四方晶格常数；showbonds[2] 随后显示面内最近邻距离为 1："
  ,"wyckoffposition supplies one nonmagnetic seed at the origin. orbitalTable[] shows that the two supplied representation coordinates px and py occupy this same physical site:" -> "wyckoffposition 在原点给出一个非磁种子；orbitalTable[] 显示输入表示的 px、py 两个坐标位于同一个物理格点："
  ,"symminformation supplies the complete ordered {E,C4z,C2z,C4z^-1} group. The representation table shows all four operations together with their supplied two-dimensional orbital matrices:" -> "symminformation 给出完整有序群 {E,C4z,C2z,C4z^-1}；表示表按同一顺序显示四个操作及其输入的二维轨道矩阵："
  ,"repinformation supplies the exact px/py matrices for E, C4z, C2z, and C4z^-1. showSymmetryRepresentations[All] displays all four matrices in the completed model:" -> "repinformation 按顺序给出 E、C4z、C2z 和 C4z^-1 的精确 px/py 矩阵；showSymmetryRepresentations[All] 显示完整模型中的四个矩阵："
  ,"orbitalLabels -> {{\"px\",\"py\"}} names the two coordinates of the supplied representation. orbitalTable[] displays those labels on the model orbitals:" -> "orbitalLabels -> {{\"px\",\"py\"}} 为输入表示的两个坐标命名；orbitalTable[] 在模型轨道上显示这些标签："
  ,"InitialBondShells -> 2 prepares the onsite and nearest-neighbour shells for this complete tetragonal model. showbonds[2] displays the four in-plane nearest neighbours:" -> "InitialBondShells -> 2 为完整四方模型准备 onsite 与最近邻两个 hopping 层；showbonds[2] 显示四个面内最近邻："
  ,"RepresentationMode -> \"DirectProduct\" keeps the supplied px and py coordinates on both C4-related sites. orbitalTable[] therefore contains four Hamiltonian basis states:" -> "RepresentationMode -> \"DirectProduct\" 在两个 C4 相关格点上都保留输入的 px、py 坐标，因此 orbitalTable[] 含有四个 Hamiltonian 基态："
  ,"RepresentationMode -> \"Induced\" needs only the px representation of the reference-site stabilizer {E,C2z}. C4 transports that one coordinate to the second site, so orbitalTable[] contains only the two symmetry-required basis states and records C4z as the transport operation:" -> "RepresentationMode -> \"Induced\" 只需输入参考点稳定子 {E,C2z} 上的 px 表示；C4 把这一个表示坐标搬运到第二个格点，因此 orbitalTable[] 只含对称性要求的两个基态，并把 C4z 记录为搬运操作："
  ,"MaxCyclotomicDegree bounds the exact field degree, not a numerical tolerance. The fourth-root field has degree two, so limit two succeeds and limit one fails:" -> "MaxCyclotomicDegree 限制的是精确数域次数，而不是数值容差。四次单位根数域的次数为 2，因此上限 2 成功，而上限 1 失败："

  ,"latticeold supplies the original implementation's direct-lattice row vectors. This complete call uses the graphene lattice and displays the resulting onsite Hamiltonian:" -> "latticeold 提供原版实现按行排列的正格矢；这个完整调用使用石墨烯晶格并显示所得 onsite Hamiltonian："
  ,"lattparold -> {aold->1,cold->3} evaluates the symbolic graphene lattice parameters before the original model is initialized:" -> "lattparold -> {aold->1,cold->3} 在初始化原版模型前对石墨烯的符号晶格参数求值："
  ,"wyckoffpositionold supplies the original implementation's fractional seed and magnetic moment. The graphene seed generates the two-site onsite Hamiltonian:" -> "wyckoffpositionold 提供原版实现的分数坐标种子与磁矩；石墨烯种子生成双格点 onsite Hamiltonian："
  ,"symminformationold supplies the complete ordered original-backend symmetry list returned by msgopold:" -> "symminformationold 提供 msgopold 返回的完整有序原版后端对称操作列表："
  ,"basisFunctionsold -> {{\"pz\"}} assigns one original-backend pz orbital to the graphene Wyckoff seed:" -> "basisFunctionsold -> {{\"pz\"}} 为石墨烯 Wyckoff 种子指定一个原版后端 pz 轨道："
  ,"debugQold -> True uses the original Boolean diagnostic setting while leaving the onsite Hamiltonian unchanged:" -> "debugQold -> True 使用原版布尔诊断设置，同时保持 onsite Hamiltonian 不变："

  ,"lattice -> DiagonalMatrix[{2,3,4}] supplies numerical direct-lattice vectors explicitly. The current model stores those three row vectors unchanged:" -> "lattice -> DiagonalMatrix[{2,3,4}] 显式提供数值正格矢；当前模型原样保存这三个行向量："
  ,"lattpar -> {a->2,b->3,c->4} evaluates the symbolic lattice before model construction. The stored numerical lattice shows the substitutions:" -> "lattpar -> {a->2,b->3,c->4} 在构造模型前对符号晶格求值；保存的数值晶格显示了这些替换："
  ,"wyckoffposition -> {{{1/4,0,0},{0,0,0}}} places the seed at fractional coordinate {1/4,0,0}. With only the identity operation, the resulting site orbit contains that one position:" -> "wyckoffposition -> {{{1/4,0,0},{0,0,0}}} 把种子放在分数坐标 {1/4,0,0}；只有单位操作时，所得格点轨道仅含该位置："
  ,"symminformation -> {E,C2z} supplies a complete ordered two-element group when GenerateSymmetryGroup -> False. The stored records retain that order and their unitary flags:" -> "当 GenerateSymmetryGroup -> False 时，symminformation -> {E,C2z} 提供完整有序二元群；保存的记录保留该顺序及其幺正标志："
  ,"basisFunctions -> {{\"s\",\"px\"}} assigns two local orbitals to the seed. The compiled representation dimension and stored basis specification both reflect that choice:" -> "basisFunctions -> {{\"s\",\"px\"}} 为种子指定两个局域轨道；编译后的表示维数和保存的基说明都会反映这一选择："
  ,"debugQ -> True uses the supported Boolean compatibility setting. The initialized session records the value without changing the representation algorithm:" -> "debugQ -> True 使用受支持的布尔兼容设置；初始化后的会话记录该值，但不改变表示算法："
  ,"InitialBondShells -> 2 prepares exactly two complete bond shells in this one-site model. The session exposes the compiled shell count directly:" -> "InitialBondShells -> 2 在这个单格点模型中准确准备两个完整近邻 hopping 层；会话直接给出已编译层数："
  ,"lattpar -> {a->2,c->4} evaluates the default symbolic lattice before model construction. The stored numerical lattice shows the substitutions:" -> "lattpar -> {a->2,c->4} 在构造模型前对默认符号晶格求值；保存的数值晶格显示了这些替换："
  ,"wyckoffposition -> {{{0,0,0},{0,0,0}}} places one nonmagnetic seed at the origin. With the default identity symmetry, the site orbit contains that one position:" -> "wyckoffposition -> {{{0,0,0},{0,0,0}}} 在原点放置一个非磁种子；在默认单位对称性下，格点轨道只包含这一个位置："
  ,"basisFunctions -> {\"s\",\"px\"} assigns two local orbitals to the default seed. The result shows the Hamiltonian dimension and stored basis:" -> "basisFunctions -> {\"s\",\"px\"} 为默认种子指定两个局域轨道；结果显示 Hamiltonian 维数和保存的基："
  ,"InitialBondShells -> 2 prepares exactly two complete bond shells. The session exposes the compiled shell count directly:" -> "InitialBondShells -> 2 准确准备两个完整近邻 hopping 层；会话直接给出已编译层数："
  ,"GenerateSymmetryGroup -> True closes the supplied generators. The default identity generator produces a one-element group:" -> "GenerateSymmetryGroup -> True 对输入生成元取闭包；默认单位生成元得到一个一元群："
  ,"RepresentationMode -> \"DirectProduct\" uses the direct-product construction. The result shows the selected mode and one-dimensional default representation:" -> "RepresentationMode -> \"DirectProduct\" 使用直积构造；结果显示所选模式和一维默认表示："
  ,"RepresentationMode -> \"Induced\" asks MagneticTB to derive the site representation automatically. The result shows the selected mode and dimension:" -> "RepresentationMode -> \"Induced\" 让 MagneticTB 自动导出格点表示；结果显示所选模式和维数："
  ,"SiteLocalData supplies the reference-site matrix explicitly in Induced mode. For the default one-site identity model, only one 1 by 1 matrix is needed:" -> "SiteLocalData 在 Induced 模式中显式提供参考点矩阵；对于默认单格点单位模型，只需一个 1×1 矩阵："

  ,"\"prec\" -> 10^-12 chops parsed real and imaginary matrix entries below that threshold. The 10^-13 off-diagonal entries become exact zeros:" -> "\"prec\" -> 10^-12 把解析后绝对值低于该阈值的矩阵实部和虚部截为零；10^-13 的非对角元变成精确零："
  ,"\"prec\" -> 10^-12 chops symmetry-generated hopping entries below that threshold. The 10^-13 non-onsite hopping becomes an exact zero:" -> "\"prec\" -> 10^-12 把低于该阈值的对称性生成 hopping 元截为零；10^-13 的非 onsite hopping 变成精确零："
  ,"\"ncell\" -> 0 keeps only the onsite translation after the complete HR file has been parsed and validated:" -> "\"ncell\" -> 0 在完整 HR 文件解析并验证后只保留 onsite 平移："

  ,"\"Hermitian\" -> False selects a non-Hermitian shell with independent forward and reverse hoppings instead of identifying them by Hermiticity:" -> "\"Hermitian\" -> False 选择前向与反向跃迁彼此独立的非厄米近邻 hopping 层，而不通过厄米性把二者等同："
  ,"\"KernelMethod\" -> \"Stacked\" selects the shell solved by the stacked exact null-space method. The same option is used when displaying its hopping origins:" -> "\"KernelMethod\" -> \"Stacked\" 选择由堆叠式精确零空间方法求解的近邻 hopping 层；显示跃迁来源时使用同一选项："
  ,"\"ValidationLevel\" -> \"Full\" selects the shell solved with the full residual, independence, and rank-nullity checks. The shell is solved explicitly before its origins are displayed:" -> "\"ValidationLevel\" -> \"Full\" 选择经过完整残差、线性独立性及秩-零度检查的近邻 hopping 层；显示来源前先显式求解该层："

  ,StringJoin[
    "\"wcc\" -> {{1/3,2/3,0},{2/3,1/3,0}} explicitly supplies the two graphene orbital ",
    "centers in the row and column order of grapheneHamiltonian. Use this form to convert a ",
    "user-supplied Hamiltonian that was not built by the current init session:"
  ] -> "\"wcc\" -> {{1/3,2/3,0},{2/3,1/3,0}} 按 grapheneHamiltonian 的行列顺序显式提供两个石墨烯轨道中心。转换并非由当前 init 会话构造的用户 Hamiltonian 时，应使用这种形式："
  ,"symmetrysetold -> All uses every original-backend symmetry operation supplied to initold when constructing the nearest-neighbour shell:" -> "symmetrysetold -> All 在构造最近邻 hopping 层时使用提供给 initold 的全部原版后端对称操作："
  ,"\"CartesianCoordinates\" -> True writes the original nearest-neighbour term in Cartesian momentum coordinates:" -> "\"CartesianCoordinates\" -> True 用笛卡尔动量坐标写出原版最近邻项："

  ,"\"Lattice\" -> DiagonalMatrix[{2 a,3 a,4 a}] supplies the symbolic direct-lattice rows used when the point-operation matrices are constructed. The valid call prepares two Wannier centers:" -> "\"Lattice\" -> DiagonalMatrix[{2 a,3 a,4 a}] 提供构造点操作矩阵时使用的符号正格矢行；这个有效调用准备两个 Wannier 中心："
  ,"\"LattPar\" -> {a->1} evaluates every symbolic lattice entry before the HR symmetry data are compiled. The resulting representation contains one 2 by 2 operation matrix:" -> "\"LattPar\" -> {a->1} 在编译 HR 对称数据前对每个符号晶格分量求值；所得表示包含一个 2×2 操作矩阵："
  ,"\"WyckoffPosition\" -> {{{1/4,0,0},{0,0,0}}} places both local orbitals on the fractional center {1/4,0,0}. The returned wcc list shows that ordering:" -> "\"WyckoffPosition\" -> {{{1/4,0,0},{0,0,0}}} 把两个局域轨道都放在分数坐标中心 {1/4,0,0}；返回的 wcc 列表显示该顺序："
  ,"\"Symminformation\" -> {E} supplies the complete ordered operation list. The returned symmetry records keep its label and unitary flag:" -> "\"Symminformation\" -> {E} 提供完整有序操作列表；返回的 symmetry 记录保留其标签和幺正标志："
  ,"\"BasisFunctions\" -> {\"s\",\"px\"} fixes a two-orbital VASP-ordered local basis. The operation matrices therefore have dimension 2 by 2:" -> "\"BasisFunctions\" -> {\"s\",\"px\"} 固定按 VASP 顺序排列的双轨道局域基；因此操作矩阵维数为 2×2："
  ,"\"Software\" -> \"VASP\" selects the currently supported external orbital-order convention. The result exposes the three compiled data fields:" -> "\"Software\" -> \"VASP\" 选择当前支持的外部轨道顺序约定；结果给出三个已编译数据字段："
  ,"\"Lattice\" -> DiagonalMatrix[{2 a,3 a,4 a}] supplies the symbolic direct-lattice rows used when the point-operation matrices are constructed. The complete call returns these two fractional Wannier centers:" -> "\"Lattice\" -> DiagonalMatrix[{2 a,3 a,4 a}] 提供构造点操作矩阵时使用的符号正格矢行；完整调用返回这两个分数坐标 Wannier 中心："
  ,"\"LattPar\" -> {a->1} evaluates every symbolic lattice entry before the HR symmetry data are compiled. Display the resulting operation matrix itself:" -> "\"LattPar\" -> {a->1} 在编译 HR 对称数据前对每个符号晶格分量求值；直接显示所得操作矩阵："
  ,"\"BasisFunctions\" -> {\"s\",\"px\"} fixes a two-orbital VASP-ordered local basis. The displayed symmetry matrix acts in that ordered basis:" -> "\"BasisFunctions\" -> {\"s\",\"px\"} 固定按 VASP 顺序排列的双轨道局域基；所显示的对称矩阵作用在这个有序基上："
  ,"\"Software\" -> \"VASP\" selects the currently supported external orbital-order convention. Display the resulting centers, matrices, and ordered operations together:" -> "\"Software\" -> \"VASP\" 选择当前支持的外部轨道顺序约定；同时显示所得中心、矩阵和有序操作："

  ,"returns the database identifier of the gray magnetic space group associated with ordinary space-group number n." -> "返回与普通空间群编号 n 对应的灰磁空间群数据库标识。"
  ,"gray is an Association whose integer keys run from 1 through 230." -> "gray 是一个 Association，其整数键为 1 到 230。"
  ,"Pass the returned integer to msgop or showMSGWyckoff. It is a database identifier, not an operation list." -> "把返回的整数传给 msgop 或 showMSGWyckoff；它是数据库标识，不是操作列表。"
  ,"Look up the gray magnetic-space-group identifier associated with ordinary space group 191:" -> "查询与普通空间群 191 对应的灰磁空间群标识："
  ,StringJoin[
    "Use the identifier returned by gray[191] directly with msgop. The selected records ",
    "show two unitary operations and pure time reversal:"
  ] -> "把 gray[191] 返回的标识直接传给 msgop；所选记录显示两个幺正操作和纯时间反演："
  ,"A number outside 1 through 230 is not a key of the association:" -> "1 到 230 之外的编号不是该 Association 的键："

  ,"returns the bundled database identifier of the magnetic space group with BNS number n.m." -> "返回 BNS 编号为 n.m 的磁空间群在内置数据库中的标识。"
  ,"bnsdict is an Association containing all 1651 BNS-number pairs in the bundled magnetic-space-group database." -> "bnsdict 是一个 Association，包含内置磁空间群数据库的全部 1651 个 BNS 编号对。"
  ,StringJoin[
    "The values of typeI, gray, typeIII, and typeIV partition the same 1651 identifiers. ",
    "Use bnsdict when the complete BNS pair is already known."
  ] -> "typeI、gray、typeIII 和 typeIV 的值将同一组 1651 个标识按类型划分；已知完整 BNS 编号对时，可直接使用 bnsdict。"
  ,"Look up the database identifier of the nontrivial type-III magnetic space group with BNS number 164.87:" -> "查询非平凡 III 型磁空间群 BNS 164.87 的数据库标识："
  ,"The BNS pair 191.234 and gray[191] select the same gray magnetic space group:" -> "BNS 编号对 191.234 与 gray[191] 选中同一个灰磁空间群："
  ,StringJoin[
    "Use the identifier returned by bnsdict[{164, 87}] directly with msgop. The selected ",
    "records show unitary and antiunitary operations of this type-III group:"
  ] -> "把 bnsdict[{164, 87}] 返回的标识直接传给 msgop；所选记录显示这个 III 型群的幺正和反幺正操作："

  ,"returns the database identifier of the type-I magnetic space group with BNS number n.m." -> "返回 BNS 编号为 n.m 的 I 型磁空间群数据库标识。"
  ,"typeI is an Association with 230 BNS-number pairs as keys." -> "typeI 是一个以 230 个 BNS 编号对为键的 Association。"
  ,"Pass the returned integer to msgop or showMSGWyckoff. The key must be the pair {n,m}, not n alone." -> "把返回的整数传给 msgop 或 showMSGWyckoff；键必须是二元组 {n,m}，不能只写 n。"
  ,"Look up the type-I magnetic-space-group identifier with BNS number 191.233:" -> "查询 BNS 编号为 191.233 的 I 型磁空间群标识："
  ,"A single integer is not a typeI key; use the complete BNS pair:" -> "单个整数不是 typeI 的键；请使用完整的 BNS 二元组："

  ,"returns the database identifier of the type-III magnetic space group with BNS number n.m." -> "返回 BNS 编号为 n.m 的 III 型磁空间群数据库标识。"
  ,"typeIII is an Association with 674 BNS-number pairs as keys." -> "typeIII 是一个以 674 个 BNS 编号对为键的 Association。"
  ,"Look up the type-III magnetic-space-group identifier with BNS number 164.87:" -> "查询 BNS 编号为 164.87 的 III 型磁空间群标识："
  ,"A BNS pair absent from the bundled association returns Missing:" -> "内置 Association 中不存在的 BNS 二元组返回 Missing："

  ,"returns the database identifier of the type-IV magnetic space group with BNS number n.m." -> "返回 BNS 编号为 n.m 的 IV 型磁空间群数据库标识。"
  ,"typeIV is an Association with 517 BNS-number pairs as keys." -> "typeIV 是一个以 517 个 BNS 编号对为键的 Association。"
  ,"Look up the type-IV magnetic-space-group identifier with BNS number 75.4:" -> "查询 BNS 编号为 75.4 的 IV 型磁空间群标识："

  ,"returns the database identifier of the gray magnetic layer group associated with layer-group number n." -> "返回与层群编号 n 对应的灰磁层群数据库标识。"
  ,"graylayer is an Association whose integer keys run from 1 through 80." -> "graylayer 是一个 Association，其整数键为 1 到 80。"
  ,"Pass the returned three-integer identifier directly to mlgop." -> "把返回的三整数标识直接传给 mlgop。"
  ,"Use mlgop[graylayer[n]] to read the complete ordered operation list; graylayer[n] alone returns only the three-integer database identifier." -> "使用 mlgop[graylayer[n]] 读取完整的有序操作列表；单独调用 graylayer[n] 只返回三整数数据库标识。"
  ,"Look up the gray magnetic-layer-group identifier associated with layer group 25:" -> "查询与层群 25 对应的灰磁层群标识："
  ,"Look up gray layer group 25 and immediately read its complete eight-operation list:" -> "查询灰磁层群 25，并立即读取其完整的八个操作："
  ,"A number outside 1 through 80 is not a key of the association:" -> "1 到 80 之外的编号不是该 Association 的键："

  ,"returns the database identifier of the gray magnetic rod group associated with rod-group number n." -> "返回与杆群编号 n 对应的灰磁杆群数据库标识。"
  ,"grayrod is an Association whose integer keys run from 1 through 75." -> "grayrod 是一个 Association，其整数键为 1 到 75。"
  ,"Pass the returned three-integer identifier directly to mrgop." -> "把返回的三整数标识直接传给 mrgop。"
  ,"Use mrgop[grayrod[n]] to read the complete ordered operation list; grayrod[n] alone returns only the three-integer database identifier." -> "使用 mrgop[grayrod[n]] 读取完整的有序操作列表；单独调用 grayrod[n] 只返回三整数数据库标识。"
  ,"Look up the gray magnetic-rod-group identifier associated with rod group 25:" -> "查询与杆群 25 对应的灰磁杆群标识："
  ,"Look up gray rod group 25 and immediately read its complete eight-operation list:" -> "查询灰磁杆群 25，并立即读取其完整的八个操作："
  ,"A number outside 1 through 75 is not a key of the association:" -> "1 到 75 之外的编号不是该 Association 的键："

  ,"Use gray[n] for a gray group. Type-I, type-III, and type-IV groups use complete BNS pairs: typeI[{n,m}], typeIII[{n,m}], and typeIV[{n,m}]." -> "灰群使用 gray[n]。I 型、III 型和 IV 型磁群使用完整 BNS 二元组：typeI[{n,m}]、typeIII[{n,m}] 和 typeIV[{n,m}]。"

  ,"look up gray magnetic-space-group database identifiers." -> "查询灰磁空间群数据库标识。"
  ,"look up magnetic-space-group database identifiers by complete BNS pair." -> "按完整 BNS 编号对查询磁空间群数据库标识。"
  ,"look up type-I magnetic-space-group identifiers by BNS pair." -> "按 BNS 二元组查询 I 型磁空间群标识。"
  ,"look up type-III magnetic-space-group identifiers by BNS pair." -> "按 BNS 二元组查询 III 型磁空间群标识。"
  ,"look up type-IV magnetic-space-group identifiers by BNS pair." -> "按 BNS 二元组查询 IV 型磁空间群标识。"
  ,"look up gray magnetic-layer-group database identifiers." -> "查询灰磁层群数据库标识。"
  ,"look up gray magnetic-rod-group database identifiers." -> "查询灰磁杆群数据库标识。"

  ,"symmetry-generated P4 Chern model" -> "对称性生成的 P4 Chern 模型"
  ,"symmetry-generated gray-P4 model" -> "对称性生成的灰 P4 模型"
  ,"C4 higher-order topological-insulator corner state" -> "C4 高阶拓扑绝缘体角态"
  ,"The symmetry-generated P4 Chern model has one chiral branch crossing the bulk gap." -> "对称性生成的 P4 Chern 模型有一条穿过体能隙的手性分支。"
  ,"The symmetry-generated gray-P4 model has a counterpropagating pair of surface branches." -> "对称性生成的灰 P4 模型有一对反向传播的表面态分支。"
  ,"The occupied band of the symmetry-generated P4 Chern model winds once across the transverse Brillouin zone." -> "对称性生成的 P4 Chern 模型的占据带沿横向 Brillouin 区绕行一次。"
  ,StringJoin[
    "The two occupied branches of the symmetry-generated gray-P4 model exchange partners between ",
    "the two time-reversal-invariant endpoints of half the Brillouin zone."
  ] -> StringJoin[
    "对称性生成的灰 P4 模型有两条占据分支；它们在半个 Brillouin 区的两个",
    "时间反演不变端点之间交换配对。"
  ]
  ,StringJoin[
    "For |s1/t1| = 0.2 < 1 the open symmetry-generated C4 breathing-square model has four ",
    "near-zero-energy corner modes. Projecting a corner orbital onto that four-state ",
    "subspace selects one representative mode without relying on an ",
    "arbitrary eigenvector inside the nearly degenerate subspace."
  ] -> StringJoin[
    "当 |s1/t1| = 0.2 < 1 时，具有开放边界、由对称性生成的 C4 呼吸方格模型",
    "有四个近零能角态。把一个角点轨道投影到这个四态子空间，可选出一个可复现的",
    "代表态，而不依赖近简并子空间中任意选取的本征矢。"
  ]
  ,StringJoin[
    "The following plotSurfaceSpectrum examples use numerical H(R) blocks produced by ",
    "hoppingData from symmetry-generated Hamiltonians. Their Graphics outputs are generated ",
    "by the current MagneticTB plotting API."
  ] -> StringJoin[
    "下面的 plotSurfaceSpectrum 范例使用 hoppingData 从对称性生成 Hamiltonian 得到的",
    "数值 H(R) 分块；图形输出由当前 MagneticTB 绘图接口直接生成。"
  ]
  ,StringJoin[
    "For the P4 Chern example, transformHoppings moves the original second in-plane direction to cell ",
    "axis 3; axis 1 is then the conserved edge momentum and axis 3 is semi-infinite. The gray-P4 ",
    "four-band example uses its stored cell with axis 3 as the stacking direction. The resulting plots ",
    "show one chiral gap-crossing branch and a counterpropagating surface-state pair, respectively."
  ] -> StringJoin[
    "在 P4 Chern 范例中，transformHoppings 把原来的第二个面内方向移到晶胞第 3 轴；",
    "第 1 轴因此是守恒的边界动量，第 3 轴是半无限方向。灰 P4 四带范例直接使用其",
    "保存的晶胞，并以第 3 轴为堆垛方向。两幅图分别显示一条穿越能隙的手性分支和",
    "一对反向传播的表面态分支。"
  ]

  ,"High-symmetry tetragonal antiferromagnet: spin symmetry and SOC" ->
    "高对称四方反铁磁体：自旋对称性与 SOC"
  ,StringJoin[
    "Consider the tetragonal antiferromagnet with magnetic space group 124.360, P_c4/mcc. ",
    "Its 32 operations generate two sites, at z=0 and z=1/2, with opposite moments. We ",
    "compare the Hamiltonians and bands with and without spin-orbit coupling (SOC). Keep all ",
    "px, py, and pz spin states so that onsite SOC can mix the different orbitals and spins."
  ] -> StringJoin[
    "下面考虑磁空间群为 124.360，P_c4/mcc 的四方反铁磁模型。群中 32 个操作生成 z=0 ",
    "和 z=1/2 两个磁矩相反的位置。我们比较有、无自旋轨道耦合（SOC）时的 Hamiltonian ",
    "和能带。由于 onsite SOC 会混合不同的 p 轨道和自旋，这里保留 px、py、pz ",
    "的全部自旋分量。"
  ]
  ,"High-symmetry MSG 124.360 with SOC" ->
    "含 SOC 的高对称 MSG 124.360"
  ,StringJoin[
    "With SOC, spin rotates together with the spatial operation. Pass the complete magnetic ",
    "space group to init, then calculate the onsite term and the nearest-neighbor hopping ",
    "between the two layers:"
  ] -> StringJoin[
    "有 SOC 时，自旋旋转与空间操作相联系。把完整磁空间群传给 init，再求 onsite ",
    "项和两层之间的最近邻跃迁："
  ]
  ,StringJoin[
    "The two layers have opposite moments. The 6 by 6 matrix below is the onsite block of ",
    "the z=0 layer, with rows and columns ordered as {pxup,pxdn,pyup,pydn,pzup,pzdn}:"
  ] -> StringJoin[
    "两层的磁矩方向相反。下面的 6 × 6 矩阵是 z=0 层的 onsite 子块，行列顺序为 ",
    "{pxup,pxdn,pyup,pydn,pzup,pzdn}："
  ]
  ,"Without SOC: add continuous C-infinity spin symmetry" ->
    "无 SOC：加入连续 C-infinity 自旋对称性"
  ,StringJoin[
    "Without SOC, spin can also rotate independently about the magnetic axis. Write the ",
    "finite operations as a spin-space group and add this continuous C_infinity symmetry. ",
    "The finite representation matrices below come from the preceding model; pass them to ",
    "initfromrep together with the continuous rotation."
  ] -> StringJoin[
    "无 SOC 时，还可以绕磁化轴独立旋转自旋。把有限操作写成自旋空间群，再加上这个连续的 ",
    "C_infinity 对称性。下面使用前一个模型得到的有限群表示矩阵，与连续旋转一起传给 ",
    "initfromrep。"
  ]
  ,StringJoin[
    "Compare the two onsite matrices. C_infinity sets the px/py-pz spin-mixing terms to ",
    "zero. With onsite and nearest-neighbor hopping retained, the numbers of real parameters ",
    "are 13 with SOC and 9 without SOC:"
  ] -> StringJoin[
    "对比两个 onsite 矩阵，可以看到加入 C_infinity 后，px/py 与 pz ",
    "之间混合自旋的项变为零。保留 onsite 和最近邻跃迁时，有 SOC 的模型含 13 个实参数，无 ",
    "SOC 时为 9 个："
  ]
  ,"Band splitting along G-Z" -> "沿 G-Z 的能带劈裂"
  ,StringJoin[
    "The nearest neighbors lie along c, so plot the bands on G-Z. For clarity, keep one ",
    "nonzero hopping parameter in each model. In the SOC model, also turn on the two onsite ",
    "spin-orbital mixing terms shown in the matrix above."
  ] -> StringJoin[
    "最近邻沿 c 方向分布，下面沿 G-Z 画能带。为了看清区别，两种模型都只取一个非零跃迁参数，有 ",
    "SOC 时再打开前面矩阵中的两个 onsite 自旋轨道混合项。"
  ]
  ,"The bands with SOC are:" -> "有 SOC 时的能带为："
  ,"With C_infinity symmetry and no SOC, the bands are:" -> "无 SOC，并加入 C_infinity 对称性后，能带为："
  ,StringJoin[
    "For these parameter values, the SOC model has six doubly degenerate levels at G. ",
    "Without SOC, the degeneracies at G are {4,2,4,2}. The fourfold degeneracies occur at ",
    "this high-symmetry point; they do not extend along the whole path."
  ] -> StringJoin[
    "在这组参数下，有 SOC 时 G 点有六组二重简并能级；无 SOC 时，G 点的简并度为 ",
    "{4,2,4,2}。这里的四重简并出现在高对称点，并不是整条路径都四重简并。"
  ]
  ,"Scope of the continuous interface" -> "连续接口的适用范围"
  ,StringJoin[
    "The continuous operation used by initfromrep is one pure spin C_infinity rotation with ",
    "no spatial action. Supply the complete finite group and its matrices in matching order. ",
    "This rotation restricts the local states; it does not generate more atoms, complete the ",
    "finite group, or add neighbor hoppings."
  ] -> StringJoin[
    "initfromrep 这里接受一个没有空间作用的纯内部 C_infinity ",
    "连续旋转。有限群必须完整，表示矩阵要与群元顺序一一对应。连续旋转只约束局域态，不会生成更多原子、补全有限",
    "群或增加近邻跃迁。"
  ]
  ,"compare a type-IV MSG with the continuous C-infinity symmetry required without SOC." ->
    "对比 IV 型 MSG 与无 SOC 时所需的连续 C-infinity 自旋对称性。"
  ," \[LongDash] compare a type-IV MSG with the continuous C-infinity symmetry required without SOC." ->
    " \[LongDash] 对比 IV 型 MSG 与无 SOC 时所需的连续 C-infinity 自旋对称性。"

  ,"gray magnetic space group" -> "灰磁空间群"
  ,"database identifier" -> "数据库标识"
  ,"type-I magnetic space group" -> "I 型磁空间群"
  ,"type-III magnetic space group" -> "III 型磁空间群"
  ,"type-IV magnetic space group" -> "IV 型磁空间群"
  ,"gray magnetic layer group" -> "灰磁层群"
  ,"gray magnetic rod group" -> "灰磁杆群"
|>
