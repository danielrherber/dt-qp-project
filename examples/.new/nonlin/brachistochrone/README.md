## brachistochrone

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u,x,y,v,t_f \\
\text{minimize :} \quad & t_f \\
\text{subject to :} \quad & \begin{bmatrix} x^\prime \\ y^\prime \\ v^\prime \end{bmatrix} = t_f \begin{bmatrix}
 v \sin(u) \\ v \cos(u) \\ g \cos(u)
\end{bmatrix} \\
& (x,y,v)(0) = (0,0,0) \\
& (x,y)(1) = (x_f,y_f)
\end{align*}
```

### Reference
pp. 133-134 of A. E. Bryson Jr. and Y.-C. Ho, *Applied Optimal Control*. Taylor & Francis, 1975, isbn: 9780891162285

### Solution
A solution is available for this problem.