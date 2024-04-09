## dynamic-soaring

### Formulation
```math
\begin{align*}
\text{changing :} \quad & C_L,\phi,x,y,z,V,\gamma,\Psi,\beta,t_f \\
\text{minimize :} \quad & \beta \\
\text{subject to :} \quad & \dot{\mathbf{\xi}} = \begin{bmatrix} \dot{x} \\ \dot{y} \\ \dot{z} \\ \dot{V} \\ \dot{\gamma} \\ \dot{\Psi} \end{bmatrix}
= \begin{bmatrix}
V \cos(\gamma) \sin(\Psi) + (\beta z + W_0) \\
V \cos(\gamma) \cos(\Psi) \\
V \sin(\gamma) \\
-\frac{\rho S}{2m} (C_{D0} + K C_L^2) V^2 - g \sin(\gamma) - \beta V \sin(\gamma) \sin(\Psi) \cos(\gamma) \\
\frac{\rho S}{2m} C_L V \cos(\phi) - g \cos(\gamma)/V + \beta \sin(\gamma) \sin(\Psi) \sin(\gamma) \\
\left(\frac{\rho S}{2m} C_L V \sin(\phi) - \beta \sin(\gamma) \cos(\Psi)\right)/\cos(\gamma)
\end{bmatrix}
\\
& \ell_{\min} \leq \frac{1}{2}\frac{\rho S}{mg} C_L V^2 \leq \ell_{\max} \\
& \mathbf{\xi}_l \leq \mathbf{\xi} \leq \mathbf{\xi}_u \\
& \mathbf{u}_l \leq \mathbf{u} \leq \mathbf{u}_u \\
& (x,y,z)(0) = (x_0,y_0,z_0) \\
& (V,\gamma,\Psi)(t_f) - (V,\gamma,\Psi)(0) = (0,0,-2\pi)
\end{align*}
```

### Reference
Y. J. Zhao, "*Optimal Pattern of Glider Dynamic Soaring*," Optimal Control Applications and Methods, vol. 25, no. 2, pp. 67-89, 2004, doi: 10.1002/oca.710