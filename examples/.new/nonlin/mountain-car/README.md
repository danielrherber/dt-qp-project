## mountain-car

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u,x,v,t_f \\
\text{minimize :} \quad & t_f \\
\text{subject to :} \quad & \begin{bmatrix}
\dot{x} \\ \dot{v}
\end{bmatrix} = \begin{bmatrix}
v \\
0.001 u - 0.0025 \cos(3x)
\end{bmatrix} \\
& (x,v)(0) = (-0.5,0) \\
& x(t_f) = 0.5 \\
& |u(t)| \leq 1 \\
& v(t_f) \geq 0
\end{align*}
```

### Reference
A. A. Melnikov, A. Makmal, H. J. Briegel, "*Projective Simulation Applied to the Grid-World and the Mountain-Car Problem*", ArXiv:1405.5459 [Cs], May 2014. arXiv.org, http://arxiv.org/abs/1405.5459