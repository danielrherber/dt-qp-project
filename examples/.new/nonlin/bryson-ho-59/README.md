## bryson-ho-59

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,x_1,x_2,x_3,x_4 \\
\text{minimize :} \quad & -x_1(t_f) \\
\text{subject to :} \quad & \begin{bmatrix} \dot{x}_1 \\ \dot{x}_2 \\ \dot{x}_3 \\ \dot{x}_4 \end{bmatrix} = \begin{bmatrix}
a \cos(u_1) \\ a \sin(u_1) \\ x_1 \\ x_2
\end{bmatrix} \\
& 0 \leq x_4(t) \\
& \mathbf{x}(0) = \mathbf{0} \\
& x_2(t_f) = 0 \\
& x_4(t_f) = h
\end{align*}
```

### Reference
pp. 59-63 of A. E. Bryson Jr. and Y.-C. Ho, *Applied Optimal Control*. Taylor & Francis, 1975, isbn: 9780891162285

### Solution
A closed-form solution is available for this problem at the reference above.