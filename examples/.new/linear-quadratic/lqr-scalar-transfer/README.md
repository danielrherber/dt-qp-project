## LQR scalar transfer

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u, x \\
\text{minimize :} \quad & \int_{0}^{t_{f}} \left [ ru^{2} + qx^{2} \right ] \mathrm{d}t\\
\text{subject to :} \quad & \dot{x} = ax + bu\\
& x(0) = c\\
& x(t_{f}) = d\\
\end{align*}
```


### Reference
- A special case of this problem is in H. P. Geering, "Optimal Control with Engineering Applications", Springer, 2007, pp. 46-48, doi: 10.1007/978-3-540-69438-0
- With certain problem parameters, this problem has similar behavior to the "Hyper-Sensitive" problem in A. V. Rao and K. D. Mease, Eigenvector Approximate Dichotomic Basis Methods for Solving Hyper-Sensitive Optimal Control Problems, Optimal Control Applications and Methods, Vol. 21, No. 1., January-February, 2000, pp. 1-17.


### Solution
A closed-form solution is available for this problem.