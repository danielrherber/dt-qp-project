## hdae: heat diffusion process with inequality

### Reference
pp. 192-195 of J. T. Betts, *Practical Methods for Optimal Control and Estimation Using Nonlinear Programming*. SIAM, 2010, isbn: 9780898716887.

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_{0}, u_{\pi}, \mathbf{x} \\
\text{minimize :} \quad & \int_{0}^{5} \left [ (\frac{1}{2}\delta + q_{1}) u_{0}^{2} + \sum_{k = 1}^{n-1} \delta x_{k}^{2} + (\frac{1}{2}\delta + q_{2}) u_{\pi}^{2} \right ] \mathrm{d}t \\
\text{subject to :} \quad & \dot{x}_{1} = \frac{1}{\delta^2} (x_{2}-2x_{1} + u_{0})\\
& \dot{x}_{k} = \frac{1}{\delta^2}(x_{k+1} - 2x_{k} + x_{k-1}) \quad k = 2, \dots, n-2 \\
& \dot{x}_{n-1} = \frac{1}{\delta^2}(u_{\pi} - 2x_{n-1} + x_{n-2})\\
& g_{0}(t) \leq u_{0}(t)\\
& g_{k}(t) \leq x_{k}(t) \quad k = 1, \dots, n-1 \\
& g_{n}(t) \leq u_{\pi}(t)\\
&u_{0}(0) = 0\\
&x_{k}(0) = 0 \quad k = 1,\dots, n-1\\
&u_{\pi}(0) = 0\\
\text{where:} \quad & g_{k}(t) = c \left [  \sin(k\delta)\sin(\frac{\pi t}{5}) - a \right ] - b \\
\delta = \pi/n
\end{align*}
```

### Solution
A boundary value problem is discussed in the reference above.