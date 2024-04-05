## alp-rider

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,u_2,x_1,x_2,x_3,x_4 \\
\text{minimize :} \quad & \int_{0}^{20} \left[10^{2}(x_1^2 + x_2^2 + x_3^2 + x_4^2) + 10^{-2} (u_1^2 +u_2^2)\right] \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
-10 x_1 + u_1 + u_2 \\
-2 x_2 + u_1 + 2 u_2 \\
-3 x_3 + 5 x_4 + u_1 - u_2 \\
5 x_3 - 3 x_4 + u_1 + 3 u_2
\end{bmatrix} \\
& \mathbf{x}(0) = (2,1,2,1) \\
& \mathbf{x}(20) = (2,3,1,-2) \\
& x_1^2 + x_2^2 + x_3^2 + x_4^2 \geq 3 p(t,3,12) + 3 p(t,6,10) + 3 p(t,10,6) + 8 p(t,15,4) + 0.01 \\
\text{where :} \quad & p(t,a,b) = e^{-b(t-a)^2}
\end{align*}
```

### Reference
pp. 163-165 in J. T. Betts, "*Practical Methods for Optimal Control and Estimation Using Nonlinear Programming.*" Society for Industrial and Applied Mathematics, Jan. 2010, doi: 10.1137/1.9780898718577