## bryson-ho-63

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u,x,y \\
\text{minimize :} \quad & \int_0^{t_f} y \left[ V \cos(\theta) + u \right] \mathrm{d}t \\
\text{subject to :} \quad & \begin{bmatrix} \dot{x} \\ \dot{y} \end{bmatrix} = \begin{bmatrix} V \cos(\theta) + u \\ V \sin(\theta) \end{bmatrix} \\
& x(0) = x(t_f) = 0 \\
& y(0) = y(t_f) = 0
\end{align*}
```

### Reference
p. 63 of A. E. Bryson Jr. and Y.-C. Ho, *Applied Optimal Control*. Taylor & Francis, 1975, isbn: 9780891162285

### Solution
A partial closed-form solution is available for this problem at the reference above.