## nonlinear-1d

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u,x \\
\text{minimize :} \quad & \frac{1}{2} \int_{0}^{5} \left[ x + u^2 \right] \mathrm{d}t \\
\text{subject to :} \quad & \dot{x} = 2 x + 2 u \sqrt{x} \\
& x(0) = 2 \\
& x(5) = 1
\end{align*}
```

### Reference
D. Garg et al., "*"Direct Trajectory Optimization and Costate Estimation of General Optimal Control Problems Using a Radau Pseudospectral Method*", in AIAA Guidance, Navigation, and Control Conference, 2009, doi: 10.2514/6.2009-5989 [Online]. Available: http://dx.doi.org/10.2514/6.2009-5989