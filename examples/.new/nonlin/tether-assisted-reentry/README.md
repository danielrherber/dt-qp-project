## tether-assisted-reentry

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,x_1,x_2,x_3,x_4 \\
\text{minimize :} \quad & \frac{1}{2} \int_0^5 \left[ x_3 \left( (x_2+1)^2 + 3 \cos^2(x_1) - 1 \right) - u_1 \right]^2 \mathrm{d}t \\
\text{subject to :} \quad & \begin{bmatrix} \dot{x}_1 \\ \dot{x}_2 \\ \dot{x}_3 \\ \dot{x}_4 \end{bmatrix} = \begin{bmatrix}
x_2 \\
-2(x_2+1) \frac{x_4}{x_3} - 3\sin(x_1)\cos(x_1) \\
x_4 \\
x_3 \left( (x_2+1)^2 + 3 \cos^2(x_1) - 1 \right) - u_1
\end{bmatrix} \\
& x_4 \geq 0 \\
& 0.01 \leq u_1 \leq 12 \\
& \mathbf{x}(0) = (30 \text{ deg}, 0, 0.05, 0.05) \\
& \mathbf{x}(5) = (0, \sqrt{6}/2, 1, 0)
\end{align*}
```

### Reference
Example 1 in P. Williams, "*Comparison of Differentiation and Integration Based Direct Transcription Methods.*" AAS/AIAA Space Flight Mechanics Meeting, AAS 05-128, Jan. 2005