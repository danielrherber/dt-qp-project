## simple-suspension

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u_1,x_1,x_2,x_3,x_4,p_1,p_2 \\
\text{minimize :} \quad & \int_{0}^{t_f} \left[ w_1 x_1^2 + w_2 \dot{x}_4^2 + w_3 u^2 \right] \mathrm{d}t \\
\text{subject to :} \quad & \dot{\mathbf{x}}(t) =
\begin{bmatrix} 0 & 1 & 0 & 0 \\\frac{-k_t}{m_{us/4}} & \frac{-[p_1 + c_t]}{m_{us/4}}& \frac{p_2}{m_{us/4}}& \frac{p_1}{m_{us/4}} \\ 0 & -1 & 0 & 1 \\ 0 & \frac{p_1}{m_{s/4}} & \frac{-p_2}{m_{s/4}}& \frac{-p_1}{m_{s/4}}
\end{bmatrix} \mathbf{x}(t) +
\begin{bmatrix}
0 \\ \frac{-1}{m_{us/4}} \\0 \\ \frac{1}{m_{s/4}}
\end{bmatrix}u(t) +
\begin{bmatrix}
-1 \\ \frac{c_t}{m_{us/4}} \\ 0 \\ 0
\end{bmatrix}\dot{z}_0(t) \\
& \mathbf{x}(0) = \mathbf{0} \\
& | x_3(t) | \leq r_{\max} \\
& \mathbf{p}_{\min} \leq \mathbf{p} \leq \mathbf{p}_{\max}
\end{align*}
```

### Reference
D. R. Herber and A. K. Sundarrajan, "*On the uses of linear-quadratic methods in solving nonlinear dynamic optimization problems with direct transcription*", in ASME International Mechanical Engineering Congress & Exposition, 2020