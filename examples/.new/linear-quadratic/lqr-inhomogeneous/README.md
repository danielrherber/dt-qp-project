## LQR inhomogeneous

### Formulation
```math
\begin{align*}
\text{changing :} \quad & \mathbf{u},\mathbf{x} \\
\text{minimize :} \quad & \mathbf{x}^{T}(t_{f})\mathbf{M}\mathbf{x}(t_{f}) + \int_{t_{0}}^{t_{f}} \left [ \mathbf{x}^{T}\mathbf{Q}(t)\mathbf{x} + \mathbf{u}^{T} \mathbf{R}(t)\mathbf{u} \right ] \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \mathbf{A}(t) \mathbf{x} + \mathbf{B}(t)\mathbf{u} + \mathbf{d}(t) \\
& \mathbf{x}(t_{0}) = \mathbf{x}_0\\
\end{align*}
```

### Reference
Adapted from pp. 175-176 of A. E. Bryson Jr. and Y.-C. Ho, *Applied Optimal Control*. Taylor & Francis, 1975, isbn: 9780891162285.

### Solution
An equivalent boundary value problem is available for this problem at the reference above.