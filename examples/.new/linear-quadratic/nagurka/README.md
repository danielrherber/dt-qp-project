## nagurka

### Formulation
```math
\begin{align*}
\text{changing :} \quad & \mathbf{u}, \mathbf{x} \\
\text{minimize :} \quad & 10 \left [ x_{1}(1) \right]^{2} + \int_{0}^{1} \left [ \mathbf{x}^{T}\mathbf{I}_{n}\mathbf{x} + \mathbf{u}^{T} \mathbf{I}_{n}\mathbf{u} \right ]~\mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 & 0 & \cdots & 0 \\
																0 & 0 & 1 & \ddots & \vdots \\
																\vdots & \ddots & \ddots & \ddots & 0\\
																0 & \cdots & 0 & 0 & 1\\
																1 & -2 & 3 & \cdots & (-1)^{n+1}n	\\
												\end{bmatrix} \mathbf{x} + \mathbf{I}_{n}\mathbf{u}\\
& \mathbf{x}(0) = \begin{bmatrix} 1 & 2 & \cdot & n \end{bmatrix} \\
\end{align*}
```

### Reference
Section 6.4 of R. Luus, *Iterative Dynamic Programming*. CRC Press, 2000, isbn: 1584881488


### Solution
There is currently no exact solution.