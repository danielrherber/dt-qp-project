## tuberculosis

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1, u_2, S, T, L_1, L_2, I_1, I_2 \\
\text{minimize :} \quad & \int_0^{t_f} \left[ L_2 + I_2 + \frac{1}{2} B_1 u_1^2 + \frac{1}{2} B_2 u_2^2 \right] \mathrm{d}t \\
\text{subject to :} \quad & \dot{S} = \Lambda - \beta_1 S \frac{I_1}{N} - \beta^* S \frac{I_2}{N} - \mu S \\
& \dot{T} = u_1 r_1 L_1 - \mu T + \left( 1 - \left(1 - u_2 \right) \left( p + q \right) \right) r_2 I_1 - \beta_2 T \frac{I_1}{N} - \beta^* T \frac{I_2}{N} \\
& \dot{L}_1 = \beta_1 S \frac{I_1}{N} - \left(\mu + k+1 \right) L_1 - u_1 r_1 L_1 + \left(1 - u_2 \right) p r_2 I_1 + \beta_2 T \frac{I_1}{N} - \beta^* L_1 \frac{I_2}{N} \\
& \dot{L}_2 = \left(1 - u_1 \right)q r_2 I_1 - \left(\mu + k_2 \right)L_2 + \beta^*\left(S + L_1 + T \right) \frac{I_2}{N} \\
& \dot{I}_1 = k_1 L_1 - \left( \mu + d_1 \right)I_1 - r_2 I_1 \\
& \dot{I}_2 = k_1 L_2 - \left( \mu + d_2 \right)I_2 \\
& N = S + L_1 + I_1 + L_2 + I_2 + T \\
& u_{\min} \leq u_1(t) \leq u_{\max} \\
& u_{\min} \leq u_2(t) \leq u_{\max} \\
& \mathbf{x}(0) = \mathbf{x}_0
\end{align*}
```

### Reference
Sec. 6.16 in J. T. Betts, *Practical Methods for Optimal Control and Estimation Using Nonlinear Programming*. SIAM, 2010, isbn: 9780898716887