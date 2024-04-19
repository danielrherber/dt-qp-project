## greenhouse glimate

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u,\mathbf{x} \\
\text{minimize :} \quad & -p_{5}x_{1}(t_{f}) + p_{4}\int_{t_{0}}^{t_f} u_{1} \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}} = \begin{bmatrix} p_{1}d_{1}(t)x_{2}(t) \\ p_{2} \left [ d_{2}(t) - x_{2}(t)\right ] + p_{3}u_{1}(t)  \end{bmatrix} \\
& x(0) = x_0 \\
& 0 \leq u_{1}(t) \leq u_{\text{max}}
\end{align*}
```

### Reference
pp. 15-24 of G. van Straten, G. van Willigenburg, E. van Henten, and R. van Ooteghem, *Optimal Control of Greenhouse Cultivation*. CRC Press, 2010 [Online]. Available: http://dx.doi.org/10.1201/b10321


