## time-energy-interceptor

### Formulation
```math
\begin{align*}
\text{changing :} \quad & a_M, {\alpha}_T, {R}_{T1}, {R}_{T1}, {\alpha}_M, {R}_{M1}, {R}_{M1}, t_f \\
\text{minimize :} \quad & \int_{0}^{t_f} \left( a_M^2 + W_t \right) \mathrm{d}t \\
\text{subject to :} \quad & \begin{bmatrix}
\dot{\alpha}_T \\ \dot{R}_{T1} \\ \dot{R}_{T1}  \\ \dot{\alpha}_M \\ \dot{R}_{M1} \\ \dot{R}_{M1}
\end{bmatrix} = \begin{bmatrix}
a_T/V_T \\
V_T\cos(\alpha_T) \\
V_T\sin(\alpha_T) \\
a_M/V_M \\
V_M\cos(\alpha_M) \\
V_M\sin(\alpha_M)
\end{bmatrix} \\
& ({\alpha}_T, {R}_{T1}, {R}_{T1}, {\alpha}_M, {R}_{M1}, {R}_{M1})(0) = ({\alpha}_{T0}, {R}_{T10}, {R}_{T10}, {\alpha}_{M0}, {R}_{M10}, {R}_{M10}) \\
& (R_{T1}(t_f)-R_{M1}(t_f))^2 + (R_{T2}(t_f)-R_{M2}(t_f))^2 = r^2\\
& |u(t)| \leq \bar{a}_M
\end{align*}
```

### Reference
A. Banerjee, M. Nabi, and T. Raghunathan, "*Time-energy optimal guidance strategy for realistic interceptor using pseudospectral method*", Transactions of the Institute of Measurement and Control, vol. 42, no. 13. SAGE Publications, pp. 2361–2371, Apr. 21, 2020 [Online]. Available: http://dx.doi.org/10.1177/0142331220910919