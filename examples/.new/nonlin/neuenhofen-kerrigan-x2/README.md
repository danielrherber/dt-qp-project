## neuenhofen-kerrigan-x2

### Formulation
```math
\begin{align*}
\text{changing :} \quad & u,x_1,x_2 \\
\text{minimize :} \quad & x_2(1) \\
\text{subject to :} \quad & \dot{\mathbf{x}}(t) = \begin{bmatrix} \frac{u(t)}{2 x_1(t)} \\ 4 [x_1(t)]^4 + [u(t)]^2 \end{bmatrix} \\
& \mathbf{x}(0) = (1,0) \\
&\sqrt{0.4} \leq  x_1(t) \\
& -1 \leq u(t)
\end{align*}
```

### Reference
X2 from M. P. Neuenhofen and E. C. Kerrigan, "An integral penalty-barrier direct transcription method for optimal control," 2020, *arXiv*:2009.06217