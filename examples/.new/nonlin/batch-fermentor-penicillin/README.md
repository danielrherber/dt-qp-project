## batch-fermentor-penicillin

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,x_1,x_2,x_3,x_4,t_f \\
\text{minimize :} \quad & -x_2(t_f) x_4(t_f) \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
h_1\left(x_1,x_3\right)x_1 - u_1\left(\frac{x_1}{500 x_4} \right) \\
h_2\left(x_3\right) x_1 - 0.01x_2 - u_1\left(\frac{x_2}{500 x_4}\right) \\
-h_1\left(x_1,x_3\right)\frac{x_1}{0.47} - h_2\left(x_3\right)\frac{x_1}{1.2} - x_1\left(\frac{0.029x_3}{0.0001+x_3}\right) + \frac{u_1}{x_4}\left(1 - \frac{x_3}{500}\right) \\
\frac{u_1}{500}
\end{bmatrix} \\
& \mathbf{x}(0) = (1.5,0,0,7) \\
& 0 \leq x_1(t) \leq 40 \\
& 0 \leq x_2(t) \leq 50 \\
& 0 \leq x_3(t) \leq 25 \\
& 0 \leq x_4(t) \leq 10 \\
& 0 \leq u_1(t) \leq 50 \\
\text{where :} \quad & h_1\left(x_1,x_3\right) = 0.11\left(\frac{x_3}{0.006x_1 + x_3}\right) \\
& h_2\left(x_3\right) = 0.0055\left(\frac{x_3}{0.0001 + x_3(1+10x_3)}\right)
\end{align*}
```

### Reference
Case Study I in J. R. Banga, E. Balsa-Canto, C. G. Moles, and A. A. Alonso, "*Dynamic optimization of bioprocesses: Efficient and robust numerical strategies*", Journal of Biotechnology, vol. 117, no. 4, pp. 407–419, Jun. 2005, doi: 10.1016/j.jbiotec.2005.02.013. [Online]. Available: http://dx.doi.org/10.1016/j.jbiotec.2005.02.013