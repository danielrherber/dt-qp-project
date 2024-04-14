## simple-sasa

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,x_1,x_2,p_1 \\
\text{minimize :} \quad & -x_1(t_f) \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
x_2 \\
-\frac{p_1}{J} x_1 + \frac{1}{J} u_1
\end{bmatrix} \\
& \mathbf{x}(0) = (0,0) \\
& x_2(t_f) = 0 \\
& 0 \leq p_1 \\
& -u_{\max} \leq u_1(t) \leq u_{\max}
\end{align*}
```

### Reference
D. R. Herber and J. T. Allison, "*Unified Scaling of Dynamic Optimization Design Formulations*", in Volume 2A: 43rd Design Automation Conference, 2017, doi: 10.1115/detc2017-67676 [Online]. Available: http://dx.doi.org/10.1115/DETC2017-67676