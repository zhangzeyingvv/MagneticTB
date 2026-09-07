# MagneticTB documentation rules

本目录中的规则专门约束 MagneticTB 的 Mathematica 官方帮助、教程、示例和文档生成脚本。

## 文档的定位

- 文档面向使用 MagneticTB 建立物理模型的用户，而不是面向程序内部开发者。
- 按物理建模顺序讲解：晶格、Wyckoff 位置、对称操作、轨道或表示、键壳层、Hamiltonian、能带与结果检查。
- 不得用内部数据结构、缓存字段或实现细节代替用户真正需要输入的模型信息。
- 英文与简体中文帮助必须同步；代码、函数名、选项名和输出保持一致，只翻译解释文字。

## 语言与叙述风格

本文档延续原版《MagneticTB 使用手册》和 MagneticTB 论文的写法。语言不追求“像程序员”，而要让熟悉紧束缚物理的读者自然地看懂怎么建模。

- 每引入一个函数或选项，先说它解决什么物理问题，再说要输入什么，最后解释计算结果怎么看。
- 优先使用紧束缚、磁空间群、Wyckoff 位置、局域轨道、键层、跃迁矩阵、Bloch Hamiltonian 等物理词汇。“编译”、“schema”、“session”、“record”、“cache entry”、“static directed-bond orbit”等内部实现词汇只能出现在高级接口页或确有调试必要的段落中。
- 不得把测试报告的口吻直接放进用户帮助，例如“已认证”、“验证为 True”、“紧凑生成集”或“规范模型 Association”。这些信息可以放在页面背后的自动测试中。
- 中文必须是自然的中文叙述，不得逐词直译英文语序。专业名词保留必要的英文，例如 Hamiltonian、onsite、hopping、DirectProduct 和 Induced；普通叙述不中英文夹杂。
- 句子尽量短而完整，一段只说一个物理观点。避免连续堆叠四个以上抽象名词，也不要用“该”、“相应的”、“进行”等空泛词把简单动作写得僵硬。
- 代码后面必须解释具体输出：哪个格点或轨道对应哪一行，哪个参数是 onsite 或 hopping，相位因子从哪条键来，对称性为什么使若干矩阵元相同或为零。
- 英文页与中文页保持相同的物理逻辑，但允许根据两种语言的自然表达调整语序，不要强求逐句字面对应。

### 原手册的写法是语言基准

遇到一句话“内容正确但读起来生硬”时，以项目根目录中的《MagneticTB 使用手册》和 `MagneticTB manual.pdf` 为准。新版功能可以增加，但叙述顺序和语气应尽量保持原手册的特点：

- 先写用户眼前的物理任务，例如“要构造石墨烯的紧束缚模型，只需给出……”，不要先写程序内部经过哪些模块。
- 常用句式为“例如……”“其中……表示……”“有了这些输入，就可以……”“检查无误后，再……”。避免把几个名词压成一句抽象定义。
- 代码之后逐项解释实际输出。可以写“第一、二行分别对应两个碳原子”“`e1` 是 onsite 参数”“三个指数项来自三条最近邻键”，不要只写“结果通过验证”或“返回一个对象”。
- 用户只需要知道程序是否已经准备好键、是否已经求过某一键层。除 `CurrentModelSession` 等高级页面外，不主动介绍“静态编译”“会话 schema”“缓存记录”“规范数据”等内部名称。
- “测试通过”“已认证”“回归模型”“机器结果”“严格验收”等属于发布检查，不属于物理说明。帮助页应改写成相应的物理结论；详细数值残差留在自动测试中。
- 中文正文尽量用主动句和具体主语。优先写“`symham` 求出这一键层允许的 hopping”，少写“相应约束被施加并完成求解”。
- 英文正文沿用原英文手册的简洁写法：先说 what the user wants to calculate，再用 `For example`, `Here`, `With ... we can ...` 解释输入和结果；避免 release-note 或 internal-API 口吻。

每次大规模修改文档，都要维护 `LanguageAudit.md`。其中必须逐页写明是否修改、保留理由以及人工复查结果；自动重放成功不能代替语言复查。

## 官方函数页的内容结构

每个正式函数页应根据函数能力组织以下内容，不得只放一个最小调用：

1. **基本范例（Basic Examples）**
   - 给出用户第一次接触该函数时真正能够理解和运行的完整范例。
   - 从模型输入开始，直接计算并显示该函数最主要的实际结果。
   - 基本范例必须自洽，但不能为了缩短代码而退化为没有物理内容的平凡输入。
2. **选项（Options）**
   - 列出所有公开选项的默认值、允许值、物理意义和计算影响。
   - 对会改变结果的主要选项分别给出完整、可运行的范例，直接展示改变后的结果。
   - 不得只打印 `Options[f]` 就算完成选项说明。
   - 不得介绍未公开的内部选项或暗示函数会执行其职责范围外的操作。
3. **应用（Applications）**
   - 使用真实或具有明确物理意义的模型展示函数如何用于实际研究。
   - 应包含多个轨道、多个 Wyckoff 位置、磁性、反幺正、自旋空间群、矩形 hopping、能带等适用场景。
   - 若模型来自论文、书籍或公开数据库，应核对原始文献并在页面中给出明确引用。
4. **性质与关系（Properties & Relations，可适用时加入）**
   - 展示新旧实现、DirectProduct/Induced、不同选项或不同构造方式之间的实质关系。
   - 关系必须通过实际 Hamiltonian、表示矩阵、能带或物理约束体现，不能只比较数组长度。
5. **可能的问题（Possible Issues，可适用时加入）**
   - 明确说明壳层越界、未初始化、表示不幺正、输入群不闭合、反幺正处理等真实失败边界。
   - 示例必须显示真实失败消息或 `$Failed`，不能制造成功假象。

章节不适用于某个简单函数时可以省略，但不得用空章节、占位文字或无意义输出凑齐结构。

## `MagneticTB 入门` 核心教程

`MagneticTB 入门` 是新用户的主入口，必须是整套文档中讲解最详细的页面，不得写成函数名列表、极简速查表或对旧 notebook 的链接集合。

- 假定读者理解基本 tight-binding 物理和 Wolfram Language 语法，但不知道 MagneticTB 的数据约定和工作流。
- 先说清楚 MagneticTB 解决什么问题：从晶格、Wyckoff 位置、对称性和局域轨道出发，生成对称性允许的 real-space hopping 与 Bloch Hamiltonian，再用于能带和其他物性计算。
- 用一个完整且不平凡的基本模型贯穿主线；后续每引入一个新概念，都必须给出可独立运行的代码和真实 Output。
- 不得在入门教程中突然出现未定义的符号、Hamiltonian、动量路径、参数规则、表示矩阵或内部 session 字段。

入门教程至少必须按以下顺序详细讲解：

1. **安装与加载**
   - 说明 paclet 安装、更新、版本检查、卸载与帮助系统打开方式。
   - 解释 `Needs["MagneticTB`"]` 只负责加载，不是建模范例。
   - 说明 2.0 与同包发布的 Old 后端如何选择，并明确新旧版必须使用不同内核。
2. **第一个完整模型**
   - 逐项解释 `lattice`、`lattpar`、`wyckoffposition`、`symminformation`、`basisFunctions` 的物理含义、数据层次、单位和顺序。
   - 给出完整 `init[...]`，不使用 `init[]`、`IdentityMatrix[3]` 或隐含的前置状态。
   - 立即显示 `orbitalTable[]`、`showbonds[...]` 或其他能让用户确认模型的实际结果，并逐列说明其含义。
3. **键层与 Hamiltonian**
   - 解释键层编号是什么、`InitialBondShells` 的默认值及范围、为什么 `symham[i]` 不允许超出 `init` 已编译范围。
   - 分开调用各个 `symham[i]`，再使用 `Sum[MagneticTB`symham[i], {i, n}]` 构造总 Hamiltonian。
   - 直接显示符号 Hamiltonian 矩阵，并解释对角项、跃迁项、参数和 `Exp[I k.d]` 相位的来源。
   - 解释 `init` 之后已完成的静态编译、`symham` 首次求解与第二次缓存读取的区别，但不向入门用户暴露无需手工操作的内部字段。
4. **能带**
   - 从已构造的 Hamiltonian 出发，完整定义参数初值、高对称路径和标签。
   - 实际调用 `bandManipulate`，显示交互式输出，并说明参数全为零时和存在 onsite/hopping 时应如何理解能带。
5. **对称性与物理检查**
   - 说明对称操作、表示矩阵、幺正/反幺正标志与 Hamiltonian 协变关系之间的联系。
   - 至少展示一个具体表示矩阵和它对 Hamiltonian 的作用；反幺正关系必须明确包含复共轭。
   - 说明为什么 Gamma 点检查不能代替泛 k 点和非 Gamma 小群的验证。
6. **从基础模型到实际模型**
   - 紧接着给出 Graphene、CsCl 矩形 hopping 和一个含自旋/反幺正的模型，说明多原子、多轨道和局域维数不同时输入层次如何变化。
   - 介绍 DirectProduct 与 Induced 的适用情形和物理差别，但不在入门主线中堆砌内部算法细节；详细理论链接到专门教程。
7. **常见错误与调试**
   - 给出未运行 `init`、键层越界、轨道输入层次错误、非幺正表示等真实失败例子。
   - 说明如何根据消息停下检查输入；不得建议用户依赖隐式 fallback 或继续使用失败后的 session。

入门教程的每个主要步骤都应同时回答三个问题：“用户要输入什么”、“这个输入的物理意义是什么”、“正确的实际输出应该是什么”。

## 每个例子必须完整

- 每个函数页和教程中的模型例子都必须能够在全新的 Wolfram 内核中从头独立运行。
- `Needs["MagneticTB`"]` 只是加载语句，不能单独构成一个范例，也不能作为范例的唯一 Input。
- 若范例需要写 `Needs["MagneticTB`"]`，应把它放在完整计算的开头，并让同一范例最终产生实际结果；不得只显示返回的 `Null`。
- 禁止用 `init[]` 作为正式模型例子，也禁止用“运行前一个页面的初始化”“假设已有 h”之类的隐含前提。
- 使用新版 `init` 时，例子至少要显式给出：
  - `lattice`
  - `lattpar`
  - `wyckoffposition`
  - `symminformation`
  - `basisFunctions`
- 使用 `initfromrep` 时，例子至少要显式给出：
  - `lattice`
  - `lattpar`
  - `wyckoffposition`
  - 完整且有序的 `symminformation`
  - 与群元顺序对应的 `repinformation`
  - `orbitalLabels`
- 使用旧版时必须完整调用 `initold`，显式给出对应的 `...old` 选项；旧版与新版不得在同一个内核中混合加载。
- 若例子调用 `symham`、`showbonds`、`orbitalTable`、`bandManipulate` 或导出函数，必须在同一例子中先给出完整初始化。
- 能带例子必须在同一例子中构造 Hamiltonian、定义完整路径并调用 `bandManipulate`；不得引用未定义的 `h`、`path` 或参数规则。
- 一个 Input 单元中的多条 Wolfram Language 表达式必须全部保留，不能在生成 notebook 时只解析第一条表达式。

## 用户页面必须直接展示实际结果

- 用户帮助页的 Output 应直接回答“这个函数算出了什么”，而不是回答“返回对象有多大或是什么类型”。
- 不得把以下诊断输出当作主要范例结果：
  - `Dimensions[...]`
  - 参数个数或参数名列表
  - `Head[...]`
  - `Length[...]`
  - `Keys[...]`
  - 仅含 `True`/`False` 的成功标志
- 上述量可以用于内部自动测试和页面背后的认证，但通常不应出现在用户可见的基本范例或应用中。
- 应根据函数性质直接显示有内容的结果，例如：
  - `symham`：明确的符号 Hamiltonian 矩阵或经过适当排版的矩阵表达式；
  - `init`：初始化后直接展示模型的轨道表、键结构、对称表示或随后得到的 Hamiltonian，而不是只显示初始化成功；
  - `orbitalTable`：直接渲染完整 `Dataset`；
  - `showbonds`：直接渲染键表或键的几何图；
  - `bandManipulate`：直接显示可操作的 `Manipulate` 能带界面，必要时补充固定参数下的静态能带图；
  - 表示相关函数：直接显示表示矩阵以及它们对 Hamiltonian 的作用；
  - 导出函数：显示生成的真实文件内容或核心数据，而不是只显示文件名字符串。
- 若完整符号结果非常长，应选择一个仍有物理意义、但结果足够清晰的模型，或使用分块矩阵、表格、图形和必要解释；不得退回到只显示 `Dimensions`。
- 参数数目、矩阵维数、对易残差等认证信息可以放在“验证说明”中作为补充，不能替代实际结果。

## 晶格与对称性

- 晶格矩阵必须符合程序数据库采用的 Bravais 晶格约定。
- 不得为了省事用 `IdentityMatrix[3]` 代替非立方晶系的真实晶格。
- 七个晶系的例子必须使用各自正确的符号晶格矩阵，并通过 `lattpar` 给出 `a`、`b`、`c`、`alpha`、`beta`、`gamma` 中实际需要的参数。
- 简单立方晶格也应写成带晶格常数的矩阵，例如 `{{a,0,0},{0,a,0},{0,0,a}}`，并设置 `lattpar -> {a -> ...}`。
- 六方晶格采用 MagneticTB 的约定：`{{a,0,0},{-a/2,Sqrt[3] a/2,0},{0,0,c}}`。
- 使用 `msgop` 数据库对称操作时，应明确给出所选磁空间群及相应 Wyckoff 位置，不得把平凡群示例伪装成真实材料模型。

## 示例覆盖范围

完整教程至少应覆盖以下相互独立的模型类型：

- 简单单轨道模型，用于解释基本调用流程。
- 多个等价原子的 Wyckoff 轨道。
- 多个 Wyckoff 轨道和每个位置多个局域轨道。
- 两端局域维数不同的矩形 hopping，例如 CsCl 的 `s`–`p` 模型。
- 至少一个包含多个 `d`、`p` 轨道的材料模型。
- 含自旋的磁群模型。
- 共线与非共线离散自旋空间群模型。
- 含反幺正操作的模型。
- `DirectProduct` 与 `Induced` 两种模式；适用同一物理轨道时要比较参数数目和能带。
- 三斜、单斜、正交、四方、三方、六方、立方七个晶系。
- 完全隔离的旧版 `initold`、`symhamold`、`bandManipulateold` 示例。

## 强制收录的现有范例

以下清单是正式说明书的最低覆盖基线，不是可选素材。“收录”指将完整、可运行的 Input、实际 Output 和必要的物理解释写入帮助系统，不能只给出外部 `.nb` 文件的路径或链接。

### `Examples/GeneralExamples.nb`

- [`Examples/GeneralExamples.nb`](../../Examples/GeneralExamples.nb) 中的所有可运行模型和工作流都必须进入正式帮助，不得只挑最简单的一部分。
- 当前该 notebook 的盘点基线是 52 个 Input 单元和 35 个 Output 单元；后续源 notebook 增加范例时，强制收录清单也自动扩展，不得固定停留在这个数量。
- 当前必须覆盖的主题为：
  - Graphene；
  - MoS2 三带 tight-binding 模型，包括笛卡尔坐标形式和 `band.dat` 导出；
  - magnetic C3 Weyl point，包括 Chern 数计算；
  - magnetic cubic nodal line；
  - C4 topological-insulator 模型；
  - magnetic space group 198.11 的 spinless 模型和 spinful C4 Weyl point；
  - Wyckoff 位置显示；
  - magnetic layer group 的对称操作显示；
  - magnetic rod group 的对称操作显示。
- 迁移这些范例时必须保留原来的物理任务和实际结果；因 2.0 API 需要改写语法时，应明确记录变更，不得把复杂模型换成平凡演示。

### 已运行的模型回归范例

下列已经用于新旧算法、表示代数、Hamiltonian 协变性、能带或缓存验证的模型，必须纳入官方教程或相应函数页：

- 经验证的模型应直接作为“基本范例”、“选项”、“应用”或“性质与关系”中的正式模型，不要另建一批与函数说明脱节的“测试范例”。
- 每个函数优先复用最能说明它物理作用的已验证模型：例如 CsCl 用于矩形 hopping，P4 用于 DirectProduct/Induced，`E + PT` 用于反幺正诱导，共线/非共线 SSG 用于自旋空间群输入。
- 测试中的参数数目、残差和 `Verified -> True` 只用于页面背后认证；用户可见范例应展示该模型的 Hamiltonian、表示矩阵、键、能带或导出数据。

- `SimpleCubicS`：简单立方单 `s` 轨道；
- `GraphenePz`：Graphene `pz` 模型，含反幺正操作；
- `CsClRectangular1x3`：CsCl `s`–`p` 的 `1 x 3` 矩形 hopping；
- `SpinfulSquarePxPy`：有自旋的方格 `px/py` 模型；
- `P4Direct` 与 `P4Induced`：同一 P4 物理轨道的 DirectProduct/Induced 对照；
- `FeSeTwoWyckoffRectangular3x2`：两个 Wyckoff 轨道、多轨道和 `3 x 2` 矩形 hopping；
- `CollinearDiscreteSSG`：8 元共线离散自旋空间群；
- `NoncollinearOctahedralSSG`：24 元非共线八面体自旋空间群；
- MoS2 11 带多轨道模型和 FeSe 10 带多 Wyckoff 模型；
- 纯内部 C4 自旋群模型和半平移与时间反演组合的自旋空间群模型；
- 只有 `E + PT` 的两点 Induced 模型，以及 `C3 x Z2^T / C3` 反幺正诱导的复局域表示范例；
- 含一个 `C_infinity` 连续生成元的 `initfromrep` 范例，包括 14.1.2.3.L.1 的 000 位置和 DirectProduct/Induced 分支。

### 七晶系 Wyckoff 数据库范例

下列 10 个数据库范例必须全部收录，并从正式 `wyckoffMSG.mx` 在运行时读取，不得把位置手写成与数据库脱节的另一套数据：

- 三斜：BNS 1.3, Wyckoff `a`；BNS 2.6, Wyckoff `i`；
- 单斜：BNS 3.3, Wyckoff `e`；BNS 3.4, Wyckoff `a`；
- 正交：BNS 16.3, Wyckoff `q`；BNS 16.3, Wyckoff `r`；
- 四方：BNS 75.3, Wyckoff `c`；
- 三方：BNS 143.3, Wyckoff `a`；
- 六方：BNS 168.111, Wyckoff `b`；
- 立方：BNS 195.3, Wyckoff `a`。

上述每个七晶系范例都必须同时提供 DirectProduct 与 Induced 的完整 `init`、`symham` 和 `bandManipulate`，并在页面背后自动验证参数子空间一致、表示矩阵精确幺正、泛 k 点 Hamiltonian 协变、非 Gamma 小群和 reciprocal sewing。可见页面应显示具体 Hamiltonian、表示矩阵和能带，而不是把一串 `True` 当成主要结果。

### 收录完整性检查

- 正式文档必须维护一份“源范例 -> 帮助 URI -> 已执行 Output -> 物理验证”对照表；清单中任何一项未映射时，文档不得视为完成。
- 同一物理模型在多个测试中重复时可以合并文档叙述，但只有在所有已测功能、选项、DirectProduct/Induced 分支、新旧版分支和实际结果都保留且可追溯时才允许合并。
- 文档发布前必须在全新内核中重跑上述全部范例。任何跳过、超时、消息、`$Failed` 或缺失 Output 都必须阻止构建。

## 输出与物理检查

- 官方帮助必须包含实际计算得到的 `Output` 单元；不得手工猜测结果。
- 每个显示的 Output 都必须由相同 Input 在全新内核中重新计算并核对。
- Hamiltonian 维数和独立参数数目只作为内部验证或实际结果后的补充信息，不得替代 Hamiltonian 本身、能带、矩阵或其他主要输出。
- 合适时还要检查：
  - Hamiltonian 的精确厄米性；
  - 表示矩阵的精确幺正性；
  - 幺正操作下的协变或对易关系；
  - 反幺正操作下包含复共轭的半线性协变关系；
  - DirectProduct 与 Induced 的参数子空间、能带和简并是否一致；
  - 非 Gamma 高对称点的小群约束和 reciprocal sewing；
  - 新旧实现的参数数目、Hamiltonian 子空间和能带。
- 不得只验证 Gamma 点就声称完整验证了空间群对称性。
- 验证失败、消息、超时或 `$Failed` 必须明确暴露；不得静默跳过、`Quiet` 掩盖或改用另一算法兜底。

## 文档生成与 Notebook 工作流

- 可编辑的唯一来源是 `Developer/Documentation/GenerateDocumentationSources.wls`、`GenerateChineseDocumentationSources.wls` 及其直接读取的示例源文件。
- `Documentation/English` 和 `Documentation/ChineseSimplified` 中的 `.nb` 是生成结果，不应手工修改。
- 修改 `.nb` 时必须遵守 Mathematica Notebook 技能的 canonicalize、validate 和一次 Wolfram 检查流程。
- 中文和 Unicode 内容必须经过 canonical ASCII notebook 转换，避免乱码。
- 文档页数、URI、Guide 链接和 Paclet 构建脚本中的预期页数必须同步。
- 构建过程中出现任何非预期 Mathematica 消息都应视为失败，不能输出“带错成功”。

### Mathematica 启动方式

- 本机已经确认存在 `wolframscript` SharedMemory/WSTP 启动器故障。Mathematica 文档生成、示例执行和构建不得先尝试 `wolframscript`，必须直接使用 `/Applications/Mathematica.app/Contents/MacOS/WolframKernel -noinit -script ...`。
- 调用时固定设置 `OMP_NUM_THREADS=1` 和 `MKL_NUM_THREADS=1`。例如：

  ```sh
  env OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    /Applications/Mathematica.app/Contents/MacOS/WolframKernel \
    -noinit -script /absolute/path/to/script.wls
  ```

- Mathematica 命令无输出、长时间高 CPU 或疑似卡住时，立即使用 `diagnose-mathematica-hangs` skill；不得重复启动同一个长命令。先区分启动器、内核、项目代码和子进程，再决定是否停止。
- 若进程表现为高 CPU 的 `wolframscript`、没有 `WolframKernel` 子进程，且采样栈位于 `mldevice_main`、`MLSharedMemoryWorld` 或 `shm_open`，说明 Wolfram Language 脚本尚未开始执行。此时只终止已经确认归属本任务的启动器进程，并改用上述直接内核命令。
- 单页文档修改应优先使用已有的定向生成入口；不得为了生成一个页面启动无关的完整文档计算。

## Notebook 排版与换行

- 官方帮助中的 Wolfram Language 代码必须按语义换行，不得把完整 `init`、长列表、长 `Association` 或能带路径挤在一条横行中。
- `init` 与 `initfromrep` 的函数名单独起始，每个公开选项单独一行，嵌套的轨道、Wyckoff 位置和表示数据按层级缩进，右括号与起始表达式对齐。例如：

  ```wl
  MagneticTB`init[
    MagneticTB`lattice -> {
      {a, 0, 0},
      {-a/2, Sqrt[3] a/2, 0},
      {0, 0, c}
    },
    MagneticTB`lattpar -> {a -> 1, c -> 3},
    MagneticTB`wyckoffposition -> {
      {{1/3, 2/3, 0}, {0, 0, 0}}
    },
    MagneticTB`symminformation -> MagneticTB`msgop[MagneticTB`gray[191]],
    MagneticTB`basisFunctions -> {{"pzup", "pzdn"}}
  ]
  ```

- 一个物理范例中，完整初始化、单壳层 Hamiltonian、总 Hamiltonian、参数赋值、能带绘图和验证应分成清晰的独立 Input/Output 单元，便于用户逐步运行和定位错误。
- 不得把相互依赖的数十条语句藏在一个巨大 Input 单元中；也不得把一个完整思路拆成大量只定义一个无法独立理解的临时符号的碎片单元。
- 长的格点矩阵、动量路径、参数规则和键列表应每个项目一行；短且结构清楚的列表不强制拆行。
- 符号矩阵输出优先使用合适的矩阵排版；宽矩阵应分块、选取有物理意义的子块或转为可读表格，不得依赖过长的水平滚动。
- 文字说明按一个完整物理观点分段，避免整页只有一个超长段落，也避免每句都分成独立 Text 单元。

## 发布前的真实渲染检查

- 自动重放 Input/Output 不能代替真实排版检查。发布前必须用 Mathematica FrontEnd 把英文和简体中文的 `init`、`GettingStarted`、`GeneralExamples`、`WyckoffCrystalSystems` 页面打印成 PDF。
- 必须用 Poppler 渲染全部 PDF 页面并检查完整页面缩略图；对宽矩阵、Dataset、表格和图形还要用较高分辨率单独复查。
- 任何被裁切的矩阵、横向滚动的 Dataset、超出页边界的表格、重叠文字、乱码或不可读图形都必须先修复；不能因为执行测试通过就忽略排版错误。
- 完成人工/视觉检查后，将所检查 notebook 的 SHA256、PDF 页数和渲染记录写入 `RenderQAReport.wl`。`ValidateDocumentationRenderQA.wls` 必须属于正式发布门禁；检查后的 notebook 一旦变化，旧记录立即失效并阻止构建。
- `RenderQAReport.wl` 只能在实际重新渲染并检查之后更新，不能由构建脚本自动把 `Passed` 写成 `True`。
- 英文与中文页必须使用相同的代码换行、单元分组、输出尺寸和图形尺寸，不得只整理其中一种语言。
- 文档生成后必须渲染检查主要页面，确认没有代码被截断、括号挤成一行、输出超出页宽、图形过大/过小或中英文布局不一致。

## 文献与引用

- 必要时可以查阅文献，但必须优先使用原始论文、官方数据库或正式出版物。
- 文献模型必须核对晶格约定、Wyckoff 位置、轨道顺序、对称群、参数定义和能带路径，不能只照抄二手代码。
- 页面中应在相关模型附近给出作者、题目、期刊、卷页、年份和 DOI 或稳定链接；不能只写“见文献”。
- 若示例只是受某文献启发而没有逐项复现，应明确写明是示意模型，不能声称复现论文结果。
- 引用不能代替完整输入；读者不打开论文也应能运行文档中的例子。

## 修改原则

- 先运行现有完整例子确认基线，再修改文档。
- 文档中发现程序错误时，先报告程序问题；未经明确授权，不得借写说明书修改正式算法。
- 复用已经存在的模型、验证函数和生成工具，不重复实现同类辅助函数。
- 删除失效示例和临时探针；保留一个清晰、可重复运行的文档示例验证入口。
- 每次交付时说明修改了哪些源文件、生成了哪些帮助页、运行了哪些例子以及实际结果。
