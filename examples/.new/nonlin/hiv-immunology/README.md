## hiv-immunology

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,u_2,x_1,x_2 \\
\text{minimize :} \quad & \int_{0}^{t_f} [-x_1 + (A_1u_i^2 +A_2u_2^2)] \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
 s_1 - \frac{s_2x_2}{b_1+x_2} - \mu x_1 - k x_1 x_2 + u_1 x_1 \\
\frac{g(1-u_2)x_2}{b_2 + x_2} - c x_1 x_2
\end{bmatrix} \\
& \mathbf{x}(0) = \mathbf{x}_0 \\
& x_{1,\min} \leq x_1(t) \leq x_{1,\max} \\
& x_{2,\min} \leq x_2(t) \leq x_{2,\max} \\
& u_{1,\min} \leq u_1(t) \leq u_{1,\max} \\
& u_{2,\min} \leq u_2(t) \leq u_{2,\max}
\end{align*}
```

### Reference
H. R. Joshi, "*Optimal control of an HIV immunology model*," Optimal Control Applications and Methods, vol. 23, no. 4, pp. 199–213, 2002, doi: 10.1002/oca.710