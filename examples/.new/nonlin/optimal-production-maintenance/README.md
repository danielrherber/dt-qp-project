## optimal-production-maintenance

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u,x,p \\
\text{minimize :} \quad & -ax(t_f)e^{-\rho t_f} - bp(t_f)e^{-\rho t_f} - \int_{0}^{t_f} e^{-\rho t} \left[ ws(t) - hx(t) - ru^2(t) - qu(t) - d - cm(t) \right] \mathrm{d}t \\
\text{subject to :} \quad & \dot{x}(t) = p(t) u(t) - s(t) \\
& \dot{p}(t) = - \left[ \alpha(t) + m(t) \right] p(t) + m(t)  \\
& (x,p)(0) = (x_0,p_0) \\
& 0 \leq u(t) \\
& 0 \leq m(t) \leq M \\
& 0 \leq x(t)
\end{align*}
```

### Reference
D. I. Cho, P. L. Abad, and M. Parlar, "*Optimal production and maintenance decisions when a system experience age-dependent deterioration*", Optim. Control Appl. Meth., vol. 14, no. 3, pp. 153–167, Jul. 1993, doi: 10.1002/oca.4660140302. [Online]. Available: http://dx.doi.org/10.1002/oca.4660140302