## space-shuttle-reentry

### Formulation
```math
\begin{align*}
\text{changing :} \quad & \alpha,\beta,r,\theta,\phi,v,\gamma,\psi,t_f \\
\text{minimize :} \quad & -\theta(t_f) \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix} \dot{r} \\ \dot{\theta} \\ \dot{\phi} \\ \dot{v} \\ \dot{\gamma} \\ \dot{\psi} \end{bmatrix} = \begin{bmatrix}
v \sin(\gamma) \\
v  \cos(\gamma) \sin(\psi) \cos^{-1}(\phi) r^{-1} \\
v \cos(\gamma) \cos(\psi) r^{-1} \\
-D - g \sin(\gamma) \\
L \cos(\beta) v^{-1} - \cos(\gamma)\left[ g v^{-1} - v r^{-1} \right] \\
L \sin(\beta) \cos^{-1}(\gamma) v^{-1} + v \cos(\gamma) \sin(\psi) \tan(\phi) r^{-1}
\end{bmatrix} \\
& \mathbf{x}(0) = (r_0, \theta_0, \phi_0, v_0, \gamma_0, \psi_0) \\
& (r,v,\gamma)(t_f) = (r_f,v_f,\gamma_f) \\
& \mathbf{x}_{\min} \leq \mathbf{x}(t) \leq \mathbf{x}_{\max} \\
& \mathbf{u}_{\min} \leq \mathbf{u}(t) \leq \mathbf{u}_{\max} \\
& q(t) \leq q_U \\
\text{where :} \quad & D = Q S C_D m^{-1}, \ C_D = C_{D0} + C_{D1} \alpha + C_{D2} \alpha^2 \\
& L = Q S C_L m^{-1}, \ C_L = C_{L0} + C_{L1} \alpha \\
& Q = \frac{1}{2} \rho v^2 \\
& \rho = \rho_0 e^{-[r - R_e]/H} \\
& g = \mu r^{-2} \\
& \hat{\alpha} = 180 \alpha \pi^{-1} \\
& q_r = 17700 \sqrt{\rho} \left[ 0.0001 v \right]^{3.07} \\
& q_a = c_{0} + c_{1}\hat{\alpha} + c_{2}\hat{\alpha}^12 + c_{3}\hat{\alpha}^3 \\
& q = q_a q_r
\end{align*}
```

### Reference
J. T. Betts, "*Practical Methods for Optimal Control and Estimation Using Nonlinear Programming*." Society for Industrial and Applied Mathematics, Jan. 2010, doi: 10.1137/1.9780898718577