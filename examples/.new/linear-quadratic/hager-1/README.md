## hager-1

### Reference
W. W. Hager, Runge-Kutta Methods in Optimal Control and the Transformed Adjoint System," *Numerische Mathematik*, vol. 87, no. 2, pp. 247-282, Dec. 2000. doi: 10.1007/s002110000178

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u, x \\
\text{minimize :} \quad & \frac{1}{2} \int_{0}^{1} \left [ \frac{5}{4}x^{2} + xu + u^2 \right ] \mathrm{d}t \\
\text{subject to :} \quad & \dot{x} = \frac{x}{2} + u\\
& x(0) = 1 \\
\end{align*}
```

### Solution
A closed-form solution is available for this problem at the reference above.