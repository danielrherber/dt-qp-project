## container-crane

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,u_2,x_1,x_2,x_3,x_4,x_5,x_6 \\
\text{minimize :} \quad & \frac{1}{2} \ \int_{0}^{t_f} [x_3^2 + x_6^2 + \rho(u_1^2 +u_2^2)] \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
x_4 \\
x_5\\
x_6 \\
u_1 + c_4x_3\\
u_2\\
-\frac{1}{x_2}[u_1 + c_5x_3 + 2x_5 x_6]
\end{bmatrix} \\
& \mathbf{x}(0) = \mathbf{x}_0 \\
& \mathbf{x}(t_f) = \mathbf{x}_f \\
& x_{4,\min} \leq x_4(t) \leq x_{4,\max} \\
& x_{5,\min} \leq x_5(t) \leq x_{5,\max} \\
& u_{1,\min} \leq u_1(t) \leq u_{1,\max} \\
& u_{2,\min} \leq u_2(t) \leq u_{2,\max}
\end{align*}
```

### Reference
D. Augustin and H. Maurer, "*Sensitivity Analysis and Real-Time Control of a Container Crane under State Constraints*," in Online Optimization of Large Scale Systems, Springer Berlin Heidelberg, 2001, pp. 69–82