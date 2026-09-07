# MagneticTB Python / Rust / Web

[English](README.md) | 简体中文

MagneticTB 根据晶格、Wyckoff 位置、磁空间群和局域轨道，构造满足对称性约束的紧束缚 Hamiltonian。

- **Rust** 是唯一的数学内核，负责群、表示、键约束、精确零空间、Hamiltonian 和物性计算。
- **Python** 提供用户接口，不重复实现 Rust 数学算法。
- **Web** 随 Python 包一起安装，在本机浏览器中提供图形界面。

当前包版本为 `0.1.0`，要求 CPython 3.9 或更高版本。详细使用方法见 [Python 网页帮助](docs/user-guide/web/index.html)。

## 安装 wheel（推荐）

普通用户只需 Python 和匹配的 wheel，不需要安装 Rust、Cargo 或系统编译工具。

下面的 `PLATFORM` 是占位符，请替换为实际 wheel 文件名中的平台标记。

安装 wheel 时会一次装好 Python 接口和 Web 运行依赖，无需选择额外安装选项。

macOS / Linux / WSL：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install ./magnetictb-0.1.0-cp39-abi3-PLATFORM.whl
```

Windows PowerShell：

```powershell
py -3 -m venv .venv
.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\python.exe -m pip install .\magnetictb-0.1.0-cp39-abi3-PLATFORM.whl
```

安装完成后检查。macOS / Linux / WSL：

```bash
.venv/bin/python -c "import magnetictb; print(magnetictb.__version__)"
```

Windows PowerShell：

```powershell
.venv\Scripts\python.exe -c "import magnetictb; print(magnetictb.__version__)"
```

### 使用 Conda 安装 wheel

已有 Conda 可以直接使用；没有时先安装与你的系统和 CPU 匹配的 [Miniforge](https://github.com/conda-forge/miniforge#install)。Windows 打开 Miniforge Prompt 或 Anaconda Prompt；macOS / Linux 使用已启用 Conda 的终端。

以下命令替代上面的 `.venv` 步骤。请在 wheel 所在目录执行安装命令，并把 `PLATFORM` 替换为实际 wheel 文件名中的平台标记：

```bash
# 已有同名环境时跳过创建，只执行 activate。
conda create -n magnetictb --override-channels -c conda-forge python=3.12 pip
conda activate magnetictb

# 一次安装 Python 接口和 Web 依赖，不需要 Rust。
python -m pip install "./magnetictb-0.1.0-cp39-abi3-PLATFORM.whl"
python -c "import magnetictb; print(magnetictb.__version__)"

# 启动本地网页。
magnetictb-web
```

打开 <http://127.0.0.1:8000/>。以后重新打开终端，先激活 `magnetictb`，再运行 `python` 或 `magnetictb-web`。不要在这个环境里再创建 `.venv`。这里由 Conda 管理 Python，pip 安装本地 MagneticTB wheel。

## 启动 Web 界面

wheel 会一并安装 Web 运行依赖。在创建 `.venv` 的目录中，使用同一环境启动：

Conda 用户先运行 `conda activate magnetictb`，再运行 `magnetictb-web`。下面带路径的命令适用于 `.venv` 用户。

macOS / Linux / WSL：

```bash
.venv/bin/magnetictb-web
```

Windows PowerShell：

```powershell
.venv\Scripts\magnetictb-web.exe
```

浏览器打开：

- Web 界面：<http://127.0.0.1:8000/>
- API 文档：<http://127.0.0.1:8000/api/docs>
- 健康检查：<http://127.0.0.1:8000/api/health>

服务默认只监听本机 `127.0.0.1:8000`，使用 `Ctrl+C` 停止。

## 基本示例

### 查询磁空间群操作

```python
from magnetictb import gray, msgop

operations = msgop(gray[191])
print(len(operations))
print(operations[0])
```

`msgop` 返回完整、有序的普通 Python 列表。每个操作的结构为：

```text
[label, rotation, translation, "F" 或 "T"]
```

### 构造 Hamiltonian

下面是一个单点、单 `s` 轨道的最小三维模型：

```python
from magnetictb import SymmetryOperation, init, symham

identity = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
)

init(
    lattice=identity,
    wyckoffposition=(((0, 0, 0), (0, 0, 0)),),
    symminformation=(SymmetryOperation("E", identity),),
    basis_functions=("s",),
    initial_bond_shells=2,
)

onsite = symham(1)
nearest = symham(2)

print(onsite)
print(nearest)
```

`init(...)` 成功时返回 `None` 并安装当前模型。`symham(shell)` 返回指定单一物理键层的完整二维符号矩阵；shell 使用从 1 开始的编号。

Graphene、三带 MoS₂、精确输入、多 shell 组合和显式表示示例见 [Python 用户指南](docs/user-guide/web/index.html)。Web 首页也内置了可直接载入和运行的 Graphene、MoS₂ 及其他晶系示例。

## 从源码构建（可选）

普通用户可以跳过这一节：安装 wheel 不需要 Rust 工具链。

三种系统的详细步骤见[从源码构建](docs/user-guide/web/guide/BuildFromSource.html)。

### 1. 准备环境

从源码编译需要：

- CPython 3.9 或更高版本；
- Rust 1.85 或更高版本；
- Cargo（使用 `rustup` 安装 Rust 时会同时安装）；
- `maturin==1.14.1`；
- Windows 上的 MSVC C++ Build Tools、macOS 上的 Xcode Command Line Tools，或 Linux 上的 GCC/Clang 与 linker。

先检查 Python、Rust 和 Cargo：

```text
python3 --version
rustc --version
cargo --version
```

Windows 可用 `py -3 --version` 代替 `python3 --version`。各系统 toolchain 的详细安装步骤见[源码构建指南](docs/user-guide/web/guide/BuildFromSource.html)。

### 2. macOS / Linux / WSL

进入本仓库的 `python/` 目录：

```bash
cd python
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install "maturin==1.14.1"
.venv/bin/maturin develop --release --locked
```

检查安装：

```bash
.venv/bin/python -c "import magnetictb; print(magnetictb.__version__)"
```

应输出：

```text
0.1.0
```

### 3. Windows PowerShell

进入本仓库的 `python` 目录：

```powershell
cd python
py -3 -m venv .venv
.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\python.exe -m pip install "maturin==1.14.1"
.venv\Scripts\maturin.exe develop --release --locked
```

检查安装：

```powershell
.venv\Scripts\python.exe -c "import magnetictb; print(magnetictb.__version__)"
```

`maturin develop` 会编译 Rust core，并把 Python 包安装到当前 `.venv`。它不会生成远程发布，也不会上传文件。

### 在 Conda 中编译

如果选择从源码安装，仍需先安装前文的系统编译工具，激活这个 Conda 环境，再从仓库根目录运行：

```bash
# 将源码构建安装到当前 Conda 环境。
cd python
python -m pip install "maturin==1.14.1"
python -m maturin develop --release --locked
```

### 构建 wheel

在已经创建好的虚拟环境中，从 `python/` 目录运行：

#### macOS / Linux / WSL

```bash
.venv/bin/maturin build --release --locked --out dist
```

#### Windows PowerShell

```powershell
.venv\Scripts\maturin.exe build --release --locked --out dist
```

生成的 wheel 位于 `python/dist/`。文件名包含操作系统、CPU 架构和 ABI 标记，例如：

```text
magnetictb-0.1.0-cp39-abi3-PLATFORM.whl
```

wheel 已包含编译后的 Rust crate 代码和 Web 静态资源。安装 wheel 的普通用户不需要 Rust、Cargo 或 maturin；但 wheel 必须与用户的操作系统和 CPU 架构匹配。

## 引用

如果在研究中使用 MagneticTB，请引用：

> Zeying Zhang, Zhi-Ming Yu, Gui-Bin Liu, Yugui Yao, “MagneticTB: A package for tight-binding model of magnetic and non-magnetic materials,” *Computer Physics Communications* **270**, 108153 (2022). <https://doi.org/10.1016/j.cpc.2021.108153>

```bibtex
@article{ZHANG2022108153,
  title   = {MagneticTB: A package for tight-binding model of
             magnetic and non-magnetic materials},
  journal = {Computer Physics Communications},
  volume  = {270},
  pages   = {108153},
  year    = {2022},
  doi     = {10.1016/j.cpc.2021.108153},
  author  = {Zeying Zhang and Zhi-Ming Yu and
             Gui-Bin Liu and Yugui Yao}
}
```
