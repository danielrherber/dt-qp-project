## second-order-singular

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,x_1,x_2,x_3 \\
\text{minimize :} \quad & x_3(t_f) \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
x_2 \\
u_1 \\
x_2^2/2 + x_3^2/2
\end{bmatrix} \\
& \mathbf{x}(0) = (0,1,0) \\
& -1 \leq u_1 \leq 1
\end{align*}
```

### Reference
G. M. Aly, "*The computation of optimal singular control*", International Journal of Control, vol. 28, no. 5, pp. 681–688, Nov. 1978, doi: 10.1080/00207177808922489. [Online]. Available: http://dx.doi.org/10.1080/00207177808922489