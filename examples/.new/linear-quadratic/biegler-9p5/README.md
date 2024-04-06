## biegler-9p5

### Formulation
```math
\begin{align*}
\text{changing :} \quad & x_1,x_2,p_1 \\
\text{minimize :} \quad & [x_1(1)]^2 \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
x_2 \\
\tau^2 x_1 - (\pi^2+\tau^2) \sin(\pi t)
\end{bmatrix} \\
& \mathbf{x}(0) = (0,p_1)
\end{align*}
```

### Reference
Example 9.5 from L. T. Biegler, "*Nonlinear Programming*." Society for Industrial and Applied Mathematics, Jan. 2010. doi: 10.1137/1.9780898719383. Available: http://dx.doi.org/10.1137/1.9780898719383

### Solution
A closed-form solution is available for this problem at the reference above.