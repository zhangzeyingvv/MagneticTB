# Rust 精确分圆域公共核

本实现直接从 Mathematica 原型移植，不经过 C++：

| Rust | Mathematica 来源 |
|---|---|
| `src/exact.rs` 有理数部分 | `Wolfram/RationalArithmetic.wl` |
| `src/exact.rs` 多项式部分 | `Wolfram/PolynomialArithmetic.wl` |
| `src/exact.rs` 分圆域部分 | `Wolfram/CyclotomicField.wl` |
| `src/matrix.rs` 矩阵部分 | `Wolfram/ExactMatrix.wl` |
| `src/matrix.rs` RREF/NullSpace | `Wolfram/ExactNullSpace.wl` |
| `src/matrix.rs` 公共核 | `Wolfram/CommonKernel.wl` |
| `src/native_backend.rs` raw 后端 | `RationalFastPath.wl`、`CyclotomicLinearAlgebraFastPath.wl`、`AdaptiveBlocking.wl`、`RealSubfieldFastPath.wl`、`CommonKernel.wl` |
| `src/fixture.rs` | `Wolfram/FixtureIO.wl` |

## 依赖

- Rust 1.85 或更新版本；
- `num-bigint`、`num-rational`、`num-traits`；
- `serde_json`。

版本已锁定在 `Cargo.lock`。精确核心不使用浮点数。

## 构建

```sh
cd rust
cargo build --release
```

Rust 不解析 `Sqrt`、`Sin`、`Cos` 等 Wolfram 表达式；这些由 Mathematica 前端编译为规范、语言无关的分圆域 JSON。

`matrix.rs` 只保留公开矩阵契约和结果收尾；`native_backend.rs` 中的 `CoefficientWorkPlan`、泛型 `raw_common_kernel`、有理/系数专用 RREF、稀疏性自适应分块和最大实子域投影直接对齐上表的 Mathematica 模块。系数乘法先提取共同分母和大整数系数，再做整数卷积与首一多项式原位约简；矩阵乘加和 RREF 乘减采用只规范化一次的融合运算。内循环不构造公开 `Element`，最终仍严格验证每个原始约束的残差。
