## cart

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u, \mathbf{x} \\
\text{minimize :} \quad & -x_{1}(t_{f}) + \frac{1}{2}\int_{0}^{t_{f}}u^2 \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix} x_{2} \\ -x_{2} + u \end{bmatrix}\\
& \mathbf{x}(0) = (0,0)
\end{align*}
```

### Reference
pp. 124-126 of D. H. Ballard, *An Introduction to Natural Computation*. MIT Press, 1999, isbn: 9780262522588.

### Solution
A closed-form solution is available for this problem at the reference above.