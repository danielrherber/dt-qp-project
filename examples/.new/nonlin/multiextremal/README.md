## multiextremal

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,x_1,x_2 \\
\text{minimize :} \quad & \int_{0}^{t_f} L(\cdot) \mathrm{d}t + M(\cdot) \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix}
f_1(\cdot) \\
f_2(\cdot)
\end{bmatrix} \\
& (x_1,x_2)(0) = (x_{10},x_{20}) \\
& u_{\min} \leq u_1(t) \leq u_{\max} \\
\text{Problem 4 :} \quad & L = 0, \quad M = x_1^2(t_f) + x_2^2(t_f) \\
& f_1 = x_2, \quad f_2 = u_1 - \sin(x_1) \\
& x_{10} = 5, \quad x_{20} = 0 \\
&  u_{\min} = -1, \quad u_{\max} = 1 \\
& t_f = 5 \\
\text{Problem 8 :} \quad & L = x_1^2 + x_2^2 + 0.1 u_1^2, \quad M = 0 \\
& f_1 = -[2+u_1][x_1+0.25] + [x_2+0.5]e^{25 x_1/[x_1+2]} \\
& f_2 = 0.5 - x_2 - [x_2+0.5]e^{25 x_1/[x_1+2]} \\
& x_{10} = 0.09, \quad x_{20} = 0.09 \\
&  u_{\min} = 0, \quad u_{\max} = 5 \\
& t_f = 0.78
\end{align*}
```

### Reference
Problems 4-12 in A. Yu. Gornov, T. S. Zarodnyuk, T. I. Madzhara, A. V. Daneeva, and I. A. Veyalko, "*A Collection of Test Multiextremal Optimal Control Problems*", in Optimization, Simulation, and Control, Springer New York, 2012, pp. 257–274 [Online]. Available: http://dx.doi.org/10.1007/978-1-4614-5131-0_16