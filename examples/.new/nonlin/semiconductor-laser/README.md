## semiconductor-laser

### Formulation
```math
\begin{align*}
\text{changing :} \quad & I,S,N,t_f \\
\text{minimize :} \quad & t_f \\
\text{subject to :} \quad & \dot{S} = - \frac{S}{\tau_p} + \Gamma G(N,S) S + \beta B N\left[ N + P_0 \right] \\
& \dot{N} = \frac{I}{q} - R(N) - \Gamma G(N,S) S \\
& (S,N)(0) = (S_0,N_0) \\
& (S,N)(t_f) = (S_f,N_f) \\
& I_{\min} \leq I(t) \leq I_{\max} \\
\text{where:} \quad & G(N,S) = G_p \left[ N - N_{\mathrm{tr}} \right] \left[ 1 - \epsilon  S \right] \\
& R(N) = AN + BN \left[ N + P_0 \right] + CN \left[ N + P_0 \right]^2
\end{align*}
```

### Reference
J.-H. R. Kim, G. L. Lippi, and H. Maurer, "*Minimizing the transition time in lasers by optimal control methods*", Physica D: Nonlinear Phenomena, vol. 191, no. 3–4, pp. 238–260, May 2004, doi: 10.1016/j.physd.2003.12.002. [Online]. Available: http://dx.doi.org/10.1016/j.physd.2003.12.002