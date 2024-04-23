## minimum energy transfer

### Formulation
```math
\begin{align*}
\text{changing :} \quad & \mathbf{u}, \mathbf{x} \\
\text{minimize :} \quad & \int_{t_{0}}^{t_{f}} \mathbf{u}^T\mathbf{I}\mathbf{u}~\mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u} \\
& \mathbf{x}(t_{0}) = \mathbf{x}_{0}\\
& \mathbf{x}(t_{f}) = \mathbf{x}_{f}
\end{align*}
```

### Solution
Solving the finite-time controllability Gramian integral provides the solution to this problem.