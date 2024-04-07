## bryson-denham

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,x_1,x_2 \\
\text{minimize :} \quad & \int_0^{1} u_1^2 \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix} x_2 \\ u_1 \end{bmatrix} \\
& \mathbf{x}(0) = (0,0) \\
& \mathbf{x}(1) = (0,-1) \\
& x_1(t) \leq \ell
\end{align*}
```

### Reference
pp. 120-123 of A. E. Bryson Jr. and Y.-C. Ho, *Applied Optimal Control*. Taylor & Francis, 1975, isbn: 9780891162285.

### Solution
A closed-form solution is available for this problem at the reference above.