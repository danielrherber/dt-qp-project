## tumor

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u, x_1, x_2 \\
\text{minimize :} \quad & x_1(t_f) \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
- \zeta x_1 \log(x_1/x_2) \\
x_2 \left( b - \mu - d x_1^{2/3} + G u \right)
\end{bmatrix} \\
& \mathbf{x}(0) = \mathbf{x}_0 \\
& 0 \leq u(t) \leq u_{\max} \\
& \int_0^{t_f} u(t) \mathrm{d}t \leq A
\end{align*}
```

### Reference
U. Ledzewicz and H. Schättler, "Analysis of optimal controls for a mathematical model of tumour anti‐angiogenesis", *Optim. Control Appl. Meth.*, vol. 29, pp. 41-57, 2008, doi:10.1002/oca.814