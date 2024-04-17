## gas absorber

### Formulation
```math
\begin{align*}
\text{changing :} \quad & \mathbf{u}, \mathbf{x} \\
\text{minimize :} \quad & \int_{0}^{t_{f}} \left [ \mathbf{x}^{T}\mathbf{I}\mathbf{x}+ \mathbf{u}^{T}\mathbf{I}\mathbf{u} \right ] \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \mathrm{tridiag}_{n}(\begin{bmatrix} 0.538998 & -1.173113 & 0.634115 \end{bmatrix})\mathbf{x} + \dots\\
& \begin{bmatrix} 0.538998 & 0 & \dots & 0 & 0 \\ 0 & 0 & \dots & 0 &0.634115   \end{bmatrix}^{T} \mathbf{u} \\
& \mathbf{x}(0) = -0.0307 - \frac{i-1}{n-1}(0.1273 - 0.0307) \quad i = 1, \dots, n \\
& \mathbf{0} \leq \mathbf{u}
\end{align*}
```

### Reference
Section 7.4.4 of R. Luus, *Iterative Dynamic Programming*. CRC Press, 2000, isbn: 1584881488



### Solution
There is currently no exact solution.