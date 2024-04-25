## Tavallaei Tousi 1

### Formulation
```math
\begin{align*}
\text{changing :} \quad & \mathbf{u}, \mathbf{x} \\
\text{minimize :} \quad & \int_{0}^{t_{f}} \left [ t^{k}\mathbf{x}^{T}\mathbf{Q}\mathbf{x} + \mathbf{u}^{T}\mathbf{R}\mathbf{u} \right ]~\mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u} \\
& \mathbf{x}(t_{0}) = \mathbf{x}_{0}\\
\end{align*}
```

### Reference
Based on the numerical example in M. A. Tavallaei and B. Tousi, "Closed Form Solution to an Optimal Control Problem by Orthogonal Polynomial Expansion," American J. of Engineering and Applied Sciences, vol. 1, no. 2, pp. 104-109, 2008. doi: 10.3844/ajeassp.2008.104.109


### Solution
A equivalent boundary value problem is available for this problem.