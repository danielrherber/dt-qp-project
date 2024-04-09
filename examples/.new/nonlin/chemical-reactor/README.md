## chemical-reactor

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u,x_1,x_2 \\
\text{minimize :} \quad & -x_2(t_f) \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
- u x_1 \\ u x_1 - \rho u^k x_2
\end{bmatrix} \\
& -0.1 \leq x_1(t) \leq 1.1 \\
& -0.1 \leq x_2(t) \leq 1.1 \\
& a_L \leq u(t) \leq a_U \\
& \mathbf{x}(0) = (1,0.01)
\end{align*}
```

### Reference
S. J. Citron, *Elements of Optimal Control*, Holt, Rinehart and Winston, New York, 1969.