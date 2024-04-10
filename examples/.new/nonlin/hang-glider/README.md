## hang-glider

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u,x,y,v_x,v_y,t_f \\
\text{minimize :} \quad & -x(1) \\
\text{subject to :} \quad & \begin{bmatrix}
x^\prime \\ y^\prime \\ v_x^\prime \\ v_y^\prime
\end{bmatrix} = t_f \begin{bmatrix}
 v_x \\
 v_y \\
\frac{1}{m}\left( -L \sin(\eta) - D \cos(\eta) \right) \\
\frac{1}{m}\left( L \cos(\eta) - D \sin(\eta) - mg \right)
\end{bmatrix} \\
& (x,y,v_x,v_y)(0) = (0,1000,13.2275675,-1.28750052) \\
& (y,v_x,v_y)(1) = (1000,13.2275675,-1.28750052) \\
& 0 \leq u(t) \leq 1.4 \\
\text{where:} \quad & t_f = p_1 \\
& C_D = C_0 + k u^2 \\
& D = \frac{1}{2} C_D \rho S v_r^2 \\
& L = \frac{1}{2} C_L \rho S v_r^2 \\
& X = \left(\frac{x}{R} - 2.5 \right)^2 \\
& V_y = v_y - u_M\left(1 - X \right) e^{-X} \\
& v_r = \sqrt{v_x^2 + V_y^2} \\
& \sin\eta = \frac{V_y}{v_r} \\
& \cos\eta = \frac{V_x}{v_r}
\end{align*}
```

### Reference
R. Bulirsch, E. Nerz, H. J. Pesch, and O. von Stryk, "*Combining Direct and Indirect Methods in Optimal Control: Range Maximization of a Hang Glider,*" in Optimal Control, Birkhäuser Basel, 1993, pp. 273–288 [Online]. Available: http://dx.doi.org/10.1007/978-3-0348-7539-4_20