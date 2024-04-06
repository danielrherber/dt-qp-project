## moon-landing

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,x_1,x_2,x_3,t_f \\
\text{minimize :} \quad & - x_3(t_f) \\
\text{subject to :} \quad & \begin{bmatrix} \dot{x}_1 \\ \dot{x}_2 \\ \dot{x}_3 \end{bmatrix} = \begin{bmatrix}
x_2 \\
-1 + u_1/x_3 \\
-u_1/2.3
\end{bmatrix} \\
& 0 \leq u_1 \leq 1.1 \\
& \mathbf{x}(0) = (1,-0.783,1) \\
& (x_1,x_2)(t_f) = (0,0)
\end{align*}
```

### Reference
Example 2 in P. Williams, "*Comparison of Differentiation and Integration Based Direct Transcription Methods.*" AAS/AIAA Space Flight Mechanics Meeting, AAS 05-128, Jan. 2005