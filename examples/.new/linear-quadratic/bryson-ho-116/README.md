## Bryson-Ho-116

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u, \mathbf{x} \\
\text{minimize :} \quad & \int_0^{t_{f}} |u|\mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix} x_2 \\ u \end{bmatrix} \\
& \mathbf{x}(0) = (x_{0}, v_{0})\\
& \mathbf{x}(t_{f}) = (0,0)\\
& |u(t)| \leq 1
\end{align*}
```

### Reference
pp. 116-117 of A. E. Bryson Jr. and Y.-C. Ho, *Applied Optimal Control*. Taylor & Francis, 1975, isbn: 9780891162285.

### Solution
A closed-form solution is available for this problem at the reference above.