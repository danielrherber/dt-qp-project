## transfer-min-time

### Formulation
```math
\begin{align*}
\text{changing :} \quad & \phi, r, \theta, v_r, v_\theta, t_f  \\
\text{minimize :} \quad & t_f \\
\text{subject to :} \quad & \begin{bmatrix} \dot{r} \\ \dot{\theta} \\ \dot{v}_r \\ \dot{v}_\theta \end{bmatrix}
=
\begin{bmatrix} v_r \\ \frac{v_\theta}{r} \\ \frac{v_\theta^2}{r} - \frac{\mu}{r^2} + \frac{T}{m_0 + \dot{m}t} \sin(\phi) \\ - \frac{v_r v_\theta}{r} + \frac{T}{m_0 + \dot{m}t} \cos(\phi) \end{bmatrix} \\
& (r,\theta,v_r,v_\theta)(0) = (1,0,0,1) \\
& (r,v_r,v_\theta)(t_f) = \left(r_f, 0, \sqrt{\frac{\mu}{r_f}} \right)
\end{align*}
```

### Reference
p. 66-69 of A. E. Bryson Jr. and Y.-C. Ho, *Applied Optimal Control*. Taylor & Francis, 1975, isbn: 9780891162285