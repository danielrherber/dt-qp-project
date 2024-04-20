## train

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1, u_2, x_1, x_2  \\
\text{minimize :} \quad & \int_{0}^{t_f} [u_1x_2 + \rho(u_1^2 +u_2^2)] \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
x_2\\
h(x_1) - (a + bx_2 + cx_2^2) + u_1 -u_2
\end{bmatrix} \\
& \mathbf{x}(0) = \mathbf{x}_0 \\
& \mathbf{x}(t_f) = \mathbf{x}_f\\
& u_{1,\min} \leq u_1(t) \leq u_{1,\max} \\
& u_{2,\min} \leq u_2(t) \leq u_{2,\max} \\
\text{where:} \quad & h(x_1) = \sum_{j=1}^{2}\left[ \frac{s_{j+1} - s_j}{\pi} \tan^{-1}\frac{x_1 - z_j}{\epsilon }\right]
\end{align*}
```

### Reference
R. J. Vanderbei, "*Case Studies in Trajectory Optimization: Trains, Planes, and Other Pastimes*," Optimization and Engineering, vol. 2, no. 2, pp. 215-243, 2001, doi: 10.1023/a:1013145328012