## simple-co-design-transfer

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u,x_1,x_2,p_1 \\
\text{minimize :} \quad & \int_0^{t_f} u^2 \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}x_2 \\ -p_1 x_1 + u \end{bmatrix} \\
& \mathbf{x}(0) = (x_0,v_0) \\
& \mathbf{x}(t_f) = (0,0)
\end{align*}
```

### Reference
D. R. Herber, J. T. Allison. 'Nested and simultaneous solution strategies for general combined plant and control design problems.' *ASME Journal of Mechanical Design*, 141(1), p. 011402, Jan 2019. doi: 10.1115/1.4040705

### Solution
A closed-form solution for fixed parameter values is available for this problem at the reference above.