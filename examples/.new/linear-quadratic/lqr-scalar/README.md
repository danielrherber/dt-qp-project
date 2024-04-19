## LQR scalar

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u, x \\
\text{minimize :} \quad & \frac{m}{2}[x(t_{f})]^{2} + \frac{1}{2} \int_{t_{0}}^{t_{f}} \left [ qx^{2} + r u^{2}  \right ] \mathrm{d}t \\
\text{subject to :} \quad & \dot{x} = ax + bu \\
& x(t_{0}) = x_0\\
\end{align*}
```

### Solution
A closed-form solution is available for this problem.