## jaddu shimemura

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u, \mathbf{x} \\
\text{minimize :} \quad & \int_{0}^{1} \left [ x_{1}^2 + x_{2}^2 + 0.005 u^{2} \right ]\mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 \\ 0 & -1 \end{bmatrix} \mathbf{x} + \begin{bmatrix} 0 \\ 1  \end{bmatrix} u \\
& \mathbf{x}(0) = (0,-1)\\
& \begin{cases}
			   & \text{Examaple 0}\\
            x_{2} - 8(t-0.5)^2 +0.5 \leq 0, & \text{Examaple 1}\\
			x_{1} - 8(t-0.5)^2 +0.5 \leq 0,	& \text{Examaple 2}\\
			x_{1}(0.5) = 0.5	& \text{Examaple 3}
		 \end{cases}
\end{align*}
```


### Reference
Multiple examples from H. Jaddu and E. Shimemura, "Solution of Constrained Linear Quadratic Optimal Control Problem Using State Parameterization," *Trans. of the Society of Instrument and Control Engineers*, vol. 34, no. 9, pp. 1164-1169, 1998. doi: 10.9746/sicetr1965.34.1164


### Solution
Closed-form solutions or finite-dimensional optimization problems are available for examples in the reference above.