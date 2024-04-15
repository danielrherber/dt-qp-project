## bryson-ho-156

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u, \mathbf{x} \\
\text{minimize :} \quad & \frac{c}{2}\left[x_{1}(t_{f}) \right]^2 + \frac{1}{2}\int_{t_{0}}^{t_{f}}u^2 \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix} x_2 \\ -w^2 x_{1} + u \end{bmatrix} \\
& x(t_{0}) = (x_{0},v_{0})\\
\end{align*}
```

### Reference
p. 156 of A. E. Bryson Jr. and Y.-C. Ho, *Applied Optimal Control*. Taylor & Francis, 1975, isbn: 9780891162285.

### Solution
A closed-form solution is available for this problem at the reference above.