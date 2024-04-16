## dt-qp-1

### Formulation
```math
\begin{align*}
\text{changing :} \quad & \mathbf{u}, \mathbf{x} \\
\text{minimize :} \quad & \int_{0}^{1} \left[ u_{1}^{2}/10 + u_{2}^{2}/10 + u_{1}x_{1} + u_{1}x_{2} + 5(x_{2}-g(t))^2 \right] \mathrm{d}t + \underset{0 \leq t \leq 1}{\mathrm{max}} \quad x_{3}(t) \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix} -1 & 2 & 0 \\ 3 & -4 & 0 \\ 1 & 2 & -1 \end{bmatrix} \mathbf{x} + \begin{bmatrix} 1 & 0 \\ -1 & 0 \\ 0 & 1/20 \end{bmatrix} \mathbf{u} \\
& x_{1}(0) = 2\\
& x_{3}(0) = 1/2\\
& x_{2}(0) - x_{2}(1) = 0\\
& \int_{0}^{1} x_{1}(t) \mathrm{d}t = 0\\
& - x_{1}(t) + u_{2}(t)/12 \leq 0\\
& x_{2}(t) \leq g(t) \\
& |u_{2}(t)| \leq 10
\end{align*}
```

### Reference
pp. 130-131 of D. R. Herber. *Advances in Combined Architecture, Plant, and Control Design.* PhD Dissertation, University of Illinois at Urbana-Champaign, Urbana, IL, USA, Dec. 2017.

### Solution
There is currently no exact solution.