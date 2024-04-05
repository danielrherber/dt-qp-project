## betts-biehn-campbell-1

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u,x_1,x_2 \\
\text{minimize :} \quad & \int_{34/15}^{4}\left[ x_1^2 + 10^{-3} u^2 \right] \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}x_2 \\ u \end{bmatrix} \\
& x_1(34/15) = 302399/50625 \\
& x_2(34/15) = 70304/3375 \\
& -x_1(t) + 15 - (t-4)^4 \leq 0
\end{align*}
```

### Reference
J. T. Betts, N. Biehn, and S. L. Campbell, Convergence of Nonconvergent IRK Discretizations of Optimal Control Problems with State Inequality Constraints," *SIAM Journal on Scientific Computing*, vol. 23, no. 6, pp. 1981-2007, 2002. doi: 10.1137/S1064827500383044

### Solution
A closed-form solution is available for this problem at the reference above.