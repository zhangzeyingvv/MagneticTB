# MagneticTB

[English](README.md) | [简体中文](README.zh-CN.md)

MagneticTB 是用于构造满足对称性约束的紧束缚 Hamiltonian 的 Wolfram Language 程序包，支持磁性与非磁性晶体的对称性、局域轨道基、符号 hopping、能带计算及相关分析和绘图。

当前 paclet 要求 Wolfram Language 12.1 或更高版本，内含英文和简体中文帮助文档。

QQ 交流群：**625192239**。

## 安装、更新与卸载

MagneticTB 从本地 `.paclet` 文件安装；本仓库不提供自动在线更新渠道。

在 Mathematica 中，将下面的路径换成本地 paclet 文件的绝对路径，然后安装：

```wl
pacletArchive =
  "/absolute/path/to/MagneticTB/build/MagneticTB-2.0.10.paclet";

PacletInstall[pacletArchive]
```

安装后先保存当前工作，完整退出并重新打开 Mathematica（不只是重启内核），再查看已安装版本并加载程序包：

```wl
PacletFind["MagneticTB"]
Needs["MagneticTB`"]
```

更新时，安装同名程序包的新版本文件，然后完整退出并重新打开 Mathematica：

```wl
PacletInstall["/absolute/path/to/MagneticTB-<new-version>.paclet"]
```

如果出现 `PacletInstall::samevers`，说明相同版本已经安装。需要重新安装这个版本时，先卸载，再安装：

```wl
PacletUninstall["MagneticTB"]
PacletInstall["/absolute/path/to/MagneticTB-<version>.paclet"]
```

如果只需卸载 MagneticTB：

```wl
PacletUninstall["MagneticTB"]
PacletFind["MagneticTB"]
```

更新、重装或卸载后，也应保存当前工作，完整退出并重新打开 Mathematica，避免继续使用之前加载的函数定义或前端中仍打开的旧帮助页。

## 打开 Mathematica 帮助文档

安装 paclet 后，先按上述步骤完整重启 Mathematica，再在文档中心搜索 `MagneticTB`，或直接打开帮助主页：

```wl
SystemOpen["paclet:MagneticTB/guide/MagneticTB"]
```

也可以直接打开某个函数的帮助页，例如：

```wl
SystemOpen["paclet:MagneticTB/ref/init"]
```

在 Notebook 中，将光标放在 MagneticTB 函数名上，按 Mathematica 的函数帮助快捷键，也可以打开对应参考页。Mathematica 根据前端语言选择已安装的英文或简体中文帮助。

## 基本示例：石墨烯 Hamiltonian 与能带

下面从 `msgop` 返回的完整有序对称操作列表出发，构造不考虑自旋的石墨烯模型。Hamiltonian 包含 onsite、最近邻和第二近邻 hopping。

```wl
Needs["MagneticTB`"];

grapheneOperations = msgop[gray[191]];

init[
  lattice -> {
    {a, 0, 0},
    {-a/2, Sqrt[3] a/2, 0},
    {0, 0, c}
  },
  lattpar -> {a -> 1, c -> 3},
  wyckoffposition -> {
    {{1/3, 2/3, 0}, {0, 0, 0}}
  },
  symminformation -> grapheneOperations,
  basisFunctions -> {{"pz"}}
];

grapheneHamiltonian = Total[symham /@ Range[3]];
MatrixForm[grapheneHamiltonian]
```

下面的路径采用倒空间分数坐标，沿 Γ–M–K–Γ 绘图，其中代码中的 `G` 表示 Γ 点。指定参数后即可画出能带：

```wl
graphenePath = {
  {{{0, 0, 0}, {0, 1/2, 0}}, {"G", "M"}},
  {{{0, 1/2, 0}, {1/3, 1/3, 0}}, {"M", "K"}},
  {{{1/3, 1/3, 0}, {0, 0, 0}}, {"K", "G"}}
};

bandplot[
  graphenePath,
  60,
  grapheneHamiltonian,
  {e1 -> 0.05, r1 -> 0.02, t1 -> 0.5}
]
```

更多使用方法见 paclet 中的 MagneticTB 帮助主页和“MagneticTB 入门”教程。

## 从源码构建 paclet

在仓库根目录运行：

```bash
sh Developer/Paclet/BuildRelease.sh
```

在 macOS 上，脚本默认使用 `/Applications/Mathematica.app/Contents/MacOS/WolframKernel`。如果内核安装在其他位置，请指定其绝对路径：

```bash
WOLFRAM_KERNEL="/absolute/path/to/WolframKernel" \
  sh Developer/Paclet/BuildRelease.sh
```

构建程序将运行所需文件和中英文帮助文档打包，生成可安装的文件：

```text
build/MagneticTB-<version>.paclet
```

其中 `<version>` 为程序包版本号。生成的 `build/` 和 `tmp/` 目录是本地构建产物，不属于人工维护的源码。

## 引用

如果在发表的工作中使用了 MagneticTB，请引用：

Zeying Zhang, Zhi-Ming Yu, Gui-Bin Liu, and Yugui Yao,
“MagneticTB: A package for tight-binding model of magnetic and non-magnetic
materials,” *Computer Physics Communications* **270**, 108153 (2022),
[doi:10.1016/j.cpc.2021.108153](https://doi.org/10.1016/j.cpc.2021.108153)。

```bibtex
@article{Zhang2022MagneticTB,
  author  = {Zhang, Zeying and Yu, Zhi-Ming and Liu, Gui-Bin and Yao, Yugui},
  title   = {MagneticTB: A package for tight-binding model of magnetic and non-magnetic materials},
  journal = {Computer Physics Communications},
  volume  = {270},
  pages   = {108153},
  year    = {2022},
  doi     = {10.1016/j.cpc.2021.108153}
}
```
