## transfer-min-fuel

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_r, u_t, r, \theta, v_r, v_t \\
\text{minimize :} \quad & \int_0^{t_f} \| \mathbf{u}(t) \|_p \mathrm{d}t \\
\text{subject to :} \quad & \dot{r} = v_r \\
& \dot{\theta} = \frac{v_t}{r} \\
& \dot{v}_r = \frac{v_t^2}{r} - \frac{1}{r^2} + u_r \\
& \dot{v}_t = -\frac{v_r v_t}{r} + u_t \\
& (r,\theta,v_r,v_t)(0) = (1, 0, 0, 1) \\
& (r,v_r,v_t)(t_f) = (4,0,0.5) \\
& \| \mathbf{u}(t) \|_q \leq u_{\max} \\
\text{where :} \quad & \frac{1}{p} + \frac{1}{q} = 1
\end{align*}
```

### Reference
I. M. Ross, Q. Gong, and P. Sekhavat, "*Low-Thrust, High-Accuracy Trajectory Optimization*," Journal of Guidance, Control, and Dynamics, vol. 30, no. 4, pp. 921–933, Jul. 2007, doi: 10.2514/1.23181