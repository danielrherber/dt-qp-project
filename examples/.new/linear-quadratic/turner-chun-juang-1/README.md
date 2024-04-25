## Turner Chun Juang 1

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u, x \\
\text{minimize :} \quad & S(x(t_{f}) - \eta)^2 + \int_{t_{0}}^{t_{f}} \left [ \frac{1}{\epsilon_{m}^{2}}(x(t)-\eta)^2 + \frac{1}{u_{m}^{2}}u^{2}(t) \right ]~\mathrm{d}t \\
\text{subject to :} \quad & \dot{x} = -\frac{1}{a}x(t) + u(t) \\
& x(0) = x_{0}\\
\end{align*}
```

### Reference
Example 6.1 of J. D. Turner, H. M. Chun, J. N. Juang, "Closed-Form Solutions for a Class of Optimal Quadratic Tracking Problems", Journal of Optimization Theory and Applications, vol. 47, no. 4, pp. 465-481, Dec. 1985, doi: 10.1007/BF00942192

### Solution
A closed-form solution is available for this problem at the reference above.