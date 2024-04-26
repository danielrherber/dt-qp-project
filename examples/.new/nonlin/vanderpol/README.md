## vanderpol

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u, x_1, x_2 \\
\text{minimize :} \quad & \int_0^{t_f} \left[ x_1^2+x_2^2+u^2 \right] \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix} x_2 \\ - x_1 + x_2 - x_1^2x_2 + u \end{bmatrix} \\
& \mathbf{x}(0) = (1,0) \\
& -0.3 \leq u(t) \leq 1
\end{align*}
```

### Reference
E. B. Canto, et al., "Restricted second order information for the solution of optimal control problems using control vector parameterization", *Journal of Process Control*, vol. 12, no. 2002, pp. 243-255, 2002, doi: 10.1016/S0959-1524(01)00008-7