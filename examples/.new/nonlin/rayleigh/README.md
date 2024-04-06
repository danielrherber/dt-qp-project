## rayleigh

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,x_1,x_2 \\
\text{minimize :} \quad & \int_0^{4.5} \left[ u_1^2 + x_1^2 \right] \mathrm{d}t \\
\text{subject to :} \quad & \begin{bmatrix} \dot{x}_1 \\ \dot{x}_2 \\ \end{bmatrix} = \begin{bmatrix}
x_2 \\
-x_1 + x_2 (1.4 - 0.14 x_2^2) + 4 u_1
\end{bmatrix} \\
& -1 \leq u_1 \leq 1 \\
& \mathbf{x}(0) = (-5,-5) \\
& \mathbf{x}(4.5) = (0,0)
\end{align*}
```

### Reference
H. Maurer and D. Augustin, "*Sensitivity Analysis and Real-Time Control of Parametric Optimal Control Problems Using Boundary Value Methods*," Online Optimization of Large Scale Systems. Springer Berlin Heidelberg, pp. 17–55, 2001. doi: 10.1007/978-3-662-04331-8_2. Available: http://dx.doi.org/10.1007/978-3-662-04331-8_2

Also, Example 4.6 in J. T. Betts, "*Practical Methods for Optimal Control and Estimation Using Nonlinear Programming*." Society for Industrial and Applied Mathematics, Jan. 01, 2010. doi: 10.1137/1.9780898718577. Available: http://dx.doi.org/10.1137/1.9780898718577 