## output tracking

### Formulation
```math
\begin{align*}
\text{changing :} \quad & \mathbf{u}, \mathbf{x} \\
\text{minimize :} \quad & \int_{t_{0}}^{t_{f}} \left [ (\mathbf{y}-\mathbf{C}\mathbf{x})^T \mathbf{Q}(\mathbf{y}-\mathbf{C}\mathbf{x}) + \mathbf{u}^{T}\mathbf{R}\mathbf{u} \right ]~\mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u} \\
& \mathbf{x}(t_{0}) = \mathbf{x}_{0}\\
\end{align*}
```

### Reference
p. 175 of A. E. Bryson Jr. and Y.-C. Ho, *Applied Optimal Control*. Taylor & Francis, 1975, isbn: 9780891162285.

### Solution
A equivalent boundary value problem is available for this problem at the reference above.