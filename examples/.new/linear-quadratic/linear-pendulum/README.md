## linear pendulum

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u,\mathbf{x} \\
\text{minimize :} \quad & -x_{1}(t_{f})\\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 \\ -k/m & 0 \end{bmatrix} \mathbf{x} + \begin{bmatrix} 0 \\ 1/m \end{bmatrix} u \\
& \mathbf{x}(0) = (x_0,v_{0})\\
& |u(t)| \leq u_{\text{max}}
\end{align*}
```

### Reference
Based on Example 7.2 of A. Bressan, "Viscosity Solutions of Hamilton-Jacobi Equations and Optimal Control Problems", Lecture Notes, pp. 30-31, 2011

### Solution
A closed-form solution is available for this problem.