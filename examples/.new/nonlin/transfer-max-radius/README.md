## transfer-max-radius

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1, u_2, r, \theta, v_r, v_\theta  \\
\text{minimize :} \quad & -r(t_f) \\
\text{subject to :} \quad & \begin{bmatrix} \dot{r} \\ \dot{\theta} \\ \dot{v}_r \\ \dot{v}_\theta \end{bmatrix}
=
\begin{bmatrix} v_r \\ \frac{v_\theta}{r} \\ \frac{v_\theta^2}{r} - \frac{\mu}{r^2} + \frac{T}{m_0 + \dot{m}t} u_1 \\ - \frac{v_r v_\theta}{r} + \frac{T}{m_0 + \dot{m}t} u_2 \end{bmatrix} \\
& (r,\theta,v_r,v_\theta)(0) = (1,0,0,1) \\
& (v_r,v_\theta)(t_f) = \left( 0, \sqrt{\frac{\mu}{r(t_f)}} \right) \\
& u_1^2 + u_2^2 = 1
\end{align*}
```

### Reference
p. 66-69 of A. E. Bryson Jr. and Y.-C. Ho, *Applied Optimal Control*. Taylor & Francis, 1975, isbn: 9780891162285