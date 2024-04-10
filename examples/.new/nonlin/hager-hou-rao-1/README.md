## hager-hou-rao-1

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,u_2,x_1,x_2 \\
\text{minimize :} \quad & \int_{0}^{1} \left[ 2 x_1^2 x_2^2 + 1.25/x_2^2 + u_2/x_2 + u_1^2 + u_2^2 \right] \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
x_1 + u_1/x_2 + x_1 x_2 u_2 \\
-x_2 [0.5 + x_2 u_2]
\end{bmatrix} \\
& \mathbf{x}(0) = (1,1)
\end{align*}
```

### Reference
W. W. Hager, H. Hou, and A. V. Rao, "*Convergence Rate for a Gauss Collocation Method Applied to Unconstrained Optimal Control*", J Optim Theory Appl, vol. 169, no. 3, pp. 801–824, Mar. 2016, doi: 10.1007/s10957-016-0929-7. [Online]. Available: http://dx.doi.org/10.1007/s10957-016-0929-7