## two-link-robot

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1, u_2, x_1, x_2, x_3, x_4, t_f \\
\text{minimize :} \quad & t_f \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
\frac{ \sin(x_3)[9/4\cos(x_3) x_1^2 + 2x_2^2 ] + 4/3 [u_1-u_2] - 3/2 \cos(x_3)u_2 }{31/36 + 9/4 \sin^2(x_3)} \\
\frac{ -[\sin(x_3)[7/2 x_1^2 + 9/4\cos(x_3) x_2^2] - 7/3u_2 + 3/2 \cos(x_3) [u_1 - u_2] ]    }{31/36 + 9/4 \sin^2(x_3)} \\
x_2 - x_1 \\
x_1
\end{bmatrix} \\
& \mathbf{x}(0) = (0,0,0.5,0) \\
& \mathbf{x}(t_f) = (0,0,0.5,0.522) \\
& \mathbf{0} \leq \mathbf{x} \\
& -1 \leq u_1 \leq 1 \\
& -1 \leq u_2 \leq 1
\end{align*}
```

### Reference
Section 12.4.2 of R. Luus, *Iterative Dynamic Programming*. CRC Press, 2000, isbn: 1584881488