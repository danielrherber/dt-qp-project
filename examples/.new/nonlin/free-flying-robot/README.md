## free-flying-robot

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,u_2,u_3,u_4,x_1,x_2,x_3,x_4,x_5,x_6 \\
\text{minimize :} \quad & \int_{0}^{t_f} [u_1 + u_2 + u_3 + u_4] \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
x_4 \\ x_5 \\ x_6 \\
[u_1 - u_2 + u_3 - u_4] \cos(x_3) \\
[u_1 - u_2 + u_3 - u_4] \sin(x_3) \\
\alpha [u_1 - u_2] - \beta [u_3 - u_4]
\end{bmatrix} \\
& \mathbf{x}(0) = (-10,-10,\pi/2,0,0,0) \\
& \mathbf{x}(t_f) = (0,0,0,0,0,0) \\
& 0 \leq u_1 \\
& 0 \leq u_2 \\
& u_1 + u_2 \leq 1 \\
& u_3 + u_4 \leq 1
\end{align*}
```

### Reference
326-330 of J. T. Betts, "*Practical Methods for Optimal Control and Estimation Using Nonlinear Programming*." Society for Industrial and Applied Mathematics, Jan. 2010, doi: 10.1137/1.9780898718577