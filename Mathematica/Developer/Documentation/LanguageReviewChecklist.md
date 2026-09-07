# MagneticTB 帮助文档逐页语言审校表

本表以原版《MagneticTB 使用手册》和 MagneticTB 论文为叙述基准。每页需同时检查英文和简体中文，并通过下列五项才可标记完成：

1. 首句先说函数解决的建模问题，不以内部数据结构开场。
2. 输入约定与物理意义都说清楚，不只改写 `Options[f]`。
3. 代码之间有自然的物理过渡，代码之后解释实际输出。
4. 中文不逐词直译，不连续堆叠“编译、会话、记录、缓存、认证”等实现词。
5. 保留完整可运行范例和真实 Output；语言修改不得削弱物理内容。

## 总览与教程

| 页面 | 首句 | 输入意义 | 结果解释 | 中文自然度 | 状态 |
|---|---:|---:|---:|---:|---|
| `guide/MagneticTB` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `tutorial/GettingStarted` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `tutorial/GeneralExamples` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `tutorial/ValidatedModels` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `tutorial/WyckoffCrystalSystems` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `tutorial/ContinuousSymmetry` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `tutorial/LegacyBackend` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `tutorial/CrystalAndKPaths` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `tutorial/RealSpaceAndTopology` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `tutorial/Corepresentations` | ✓ | ✓ | ✓ | ✓ | 已审 |

## 群、表示与初始化

| 页面 | 首句 | 输入意义 | 结果解释 | 中文自然度 | 状态 |
|---|---:|---:|---:|---:|---|
| `ref/msgop` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/init` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/initfromrep` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/CurrentModelSession` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/CompileMagneticTBInput` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/pointMatrix` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/GenerateGroup` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/getGenerator` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/brokenSymmetryInitRules` | ✓ | ✓ | ✓ | ✓ | 已审 |

## Hamiltonian、键与基底

| 页面 | 首句 | 输入意义 | 结果解释 | 中文自然度 | 状态 |
|---|---:|---:|---:|---:|---|
| `ref/symham` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/unsymham` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/symhamII` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/orbitalTable` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/hoppingData` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/transformHoppings` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/buildBlochHamiltonian` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/buildRealSpaceHamiltonian` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/buildSlabHamiltonian` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/surfaceGreenFunction` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/surfaceSpectralFunction` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/showHamiltonianBasis` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/showbonds` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/showHoppingParameters` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/showSymmetryRepresentations` | ✓ | ✓ | ✓ | ✓ | 已审 |

## 能带、数据库与接口

| 页面 | 首句 | 输入意义 | 结果解释 | 中文自然度 | 状态 |
|---|---:|---:|---:|---:|---|
| `ref/bandManipulate` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/bandplot` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/showband` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/standardKPath` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/showCrystalStructure` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/showBrillouinZone` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/compareBand` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/banddata` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/hop` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/readHR` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/symmetrizationHRInit` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/showMSGWyckoff` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/mlgop` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/mrgop` | ✓ | ✓ | ✓ | ✓ | 已审 |

## 拓扑性质与 Corep

| 页面 | 首句 | 输入意义 | 结果解释 | 中文自然度 | 状态 |
|---|---:|---:|---:|---:|---|
| `ref/wilsonLoop` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/berryPhase` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/berryCurvature` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/pointChernNumber` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/findGaplessPoints` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/getMSGElemFromMSGCorep` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/getTBBandCorep` | ✓ | ✓ | ✓ | ✓ | 已审 |

## 分圆域精确线性代数

| 页面 | 首句 | 输入意义 | 结果解释 | 中文自然度 | 状态 |
|---|---:|---:|---:|---:|---|
| `ref/CompileCyclotomicMatrices` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/CyclotomicRREF` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/CyclotomicNullSpace` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/CyclotomicCommonKernel` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/CyclotomicCommonNullSpace` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/RestoreCyclotomicExpression` | ✓ | ✓ | ✓ | ✓ | 已审 |

## 旧版独立后端

| 页面 | 首句 | 输入意义 | 结果解释 | 中文自然度 | 状态 |
|---|---:|---:|---:|---:|---|
| `ref/initold` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/symhamold` | ✓ | ✓ | ✓ | ✓ | 已审 |
| `ref/bandManipulateold` | ✓ | ✓ | ✓ | ✓ | 已审 |

## 完成条件

- 64 个英文页和 64 个中文页已纳入审校表；本次 23 个新增/变更函数与 3 个新教程成对复核。
- 所有修改从 `GenerateDocumentationSources.wls`、`GenerateChineseDocumentationSources.wls` 及翻译表生成，不手工修改最终 `.nb`。
- 中英文页面重新生成并通过 canonicalize、expected text、Wolfram 解析、实际示例回放和渲染检查。
- 独立检查智能体逐页审读语言与物理逻辑，而不只运行自动门禁。
