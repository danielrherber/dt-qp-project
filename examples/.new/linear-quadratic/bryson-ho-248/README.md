## bryson-ho-248

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u, \mathbf{x} \\
\text{minimize :} \quad & \frac{1}{2}\int_{0}^{t_{f}}x_{1}^2 \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 \\ 0 & 0 \end{bmatrix} \mathbf{x} + \begin{bmatrix} 1 \\ -1 \end{bmatrix} u  \\
& \mathbf{x}(0) = (\alpha,\beta)\\
& \mathbf{x}(t_{f}) = (0, 0)\\
& |u(t)| \leq \gamma
\end{align*}
```

### Reference
pp. 248-250 of A. E. Bryson Jr. and Y.-C. Ho, *Applied Optimal Control*. Taylor & Francis, 1975, isbn: 9780891162285.

### Solution
A relatively straightforward implicit solution is available for this problem at the reference above.