## earth-launch

### Formulation
```math
\begin{align*}
\text{changing :} \quad & \theta,x,y,v_x,v_y,m,t_f \\
\text{minimize :} \quad & t_f \\
\text{subject to :} \quad & \begin{bmatrix}
\dot{x} \\ \dot{y}  \\ \dot{v}_x  \\ \dot{v}_y \\ \dot{m}
\end{bmatrix} = \begin{bmatrix}
v_x \\
v_y \\
\frac{1}{m}\left( T \cos(\theta) - 0.5 C_D S D v_x^2 \right) \\
\frac{1}{m}\left( T \sin(\theta) - 0.5 C_D S D v_y^2 \right) - g \\
-\frac{T}{g I_{\mathrm{sp}}}
\end{bmatrix} \\
& (x,y,v_x,v_y,m)(0) = (0,0,0,0,117000) \\
& (y,v_x,v_y)(t_f) = (185000,7796.6961,0) \\
& -\pi/2 \leq \theta(t) \leq \pi/2 \\
\text{where :} \quad & D = \rho_{\mathrm{ref}}e^{-y/h_\mathrm{ref}}
\end{align*}
```

### Reference
pp. 222-226 in J. M. Longuski, J. J. Guzmán, and J. E. Prussing, *Optimal Control with Aerospace Applications*. Springer New York, 2014 [Online]. Available: http://dx.doi.org/10.1007/978-1-4614-8945-0