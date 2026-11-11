# FE_solver_2D_elasticity

二维线弹性方程（Lamé 方程组）的有限元求解器，一次三角元 / 二次三角元，MATLAB 实现。

算例：Ω = [0,1]²，−∇·σ(u) = f，∂Ω 上 u1 = u2 = 0，λ = 1、μ = 2；
解析解 u1 = sin(πx)sin(πy)、u2 = x(x−1)y(y−1)。

## 使用方法

在 MATLAB 中：

```matlab
addpath('D:\matlab\FE_solver_2D_elasticity');
run_error_table            % 输出两张收敛性表格（一次元 201 / 二次元 202）
```

单独求解某个网格：

```matlab
[solution, P, T, Pb, Tb] = FE_solver_2D_elasticity(0, 1, 0, 1, 201, 201, 16, 16);
[L2, H1, Linf] = compute_error_2D(solution, P, T, Pb, Tb, 201);
```

## 结果

| 单元类型 | ‖u−u_h‖_∞ | ‖u−u_h‖_0 | \|u−u_h\|_1 |
|---|---|---|---|
| 一次元（201） | O(h²) | O(h²) | O(h) |
| 二次元（202） | O(h³) | O(h³) | O(h²) |

完整数值见 `收敛性表格运行结果.txt`；每一处代码改动的问题诊断与说明见 `代码修改说明.md`
（源码内以 `% [修改1]` … `% [修改17]` 就地标注）。
