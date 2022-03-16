BrysonHo166
===========

Problem
-------

Consider the following :term:`LQDO` problem from pp. 166--167 in [Bry75]_:

.. math::

	\underset{u,\mathbf{x}}{\text{minimize:}} \quad & F = \int_0^{t_f} u^2 \mathrm{d}t \\
	\text{subject to:} \quad & \dot{\mathbf{x}} = \begin{bmatrix}x_2 \\ -x_1 + u \end{bmatrix} \\
	& \mathbf{x}(0) = (x_0,v_0) \\
	& \mathbf{x}(t_f) = (0,0)

Code
----

The :term:`DTQP` MATLAB code for this problem is:

.. code-block:: matlab
	:linenos:
	:name: BrysonHo166


	% problem parameters
	tf = 20; % final time
	auxdata.x0 = -0.5; auxdata.v0 = 1;

	% system dynamics
	A = [0 1;-1 0]; B = [0;1];

	% Lagrange term
	L(1).left = 1; L(1).right = 1; L(1).matrix = 1/2; % 1/2*u^2

	% simple bounds
	LB(1).right = 4; LB(1).matrix = [auxdata.x0;auxdata.v0]; % initial states
	UB(1).right = 4; UB(1).matrix = [auxdata.x0;auxdata.v0]; % initial states
	LB(2).right = 5; LB(2).matrix = [0;0]; % final states
	UB(2).right = 5; UB(2).matrix = [0;0]; % final states

	% combine structures
	setup.A = A; setup.B = B; setup.L = L; setup.LB = LB; setup.UB = UB;
	setup.t0 = 0; setup.tf = tf; setup.auxdata = auxdata;

	% solve
	[T,U,Y,P,F,in,opts] = DTQP_solve(setup,[]]);

Please see `BrysonHo166.m <https://github.com/danielrherber/dt-qp-project/blob/master/examples/linear-quadratic/bryson-ho-166/BrysonHo166.m>`_ for this example in the repository with different options and comparisons to the closed-form solution below.

Solution
--------

It can be shown that the control trajectory that minimizes the objective while satisfying the constraints is:

.. math::

	u^*(t) = -\frac{2}{t_f^2 - \sin^2(t_f)}
	\begin{bmatrix}
	x_0 \\ v_0
	\end{bmatrix}^T \begin{bmatrix}
	\sin(t_f -t) \sin(t_f) - t_f \sin(t) \\
	-\cos(t_f -t) \sin(t_f) + t_f \cos(t)
	\end{bmatrix}

with an optimal objective function value of:

.. math::

	F^* = \frac{t_f \left({v_{0}}^2+{x_{0}}^2\right)+2 {t_f}^2 v_{0} x_{0}-\cos\left(t_f\right) \sin\left(t_f\right) \left({v_{0}}^2-{x_{0}}^2\right)}{{{t_f}^2 - \sin\left(t_f\right)}^2} -2 v_{0} x_{0}

The problem parameters used are :math:`t \in [0, 20]`, :math:`x_0 = -1/2`, and :math:`v_0 = 1`.
With these parameter values, :math:`F^* = 0.059842`.
The optimal trajectories for both the control and states is shown in X.

Results
-------