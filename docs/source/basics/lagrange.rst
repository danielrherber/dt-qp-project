Lagrange Objective Terms
========================

Linear-Quadratic Lagrange Terms
-------------------------------

Lagrange terms are formed using the objective structure format.
There can be an arbitrary :math:`n_L` individual structures that are then summed to construct the complete Lagrange term that is to be minimized:

.. math::

	L = \sum_{k=1}^{n_L} L_k = \sum_{k=1}^{n_L} \left[ \int_{t_0}^{t_f} \mathbf{z}_{\mathrm{left}}^T [\mathrm{matrix}] \mathbf{z}_{\mathrm{right}} \mathrm{d}t \right]_k

The properly-constructed structures then need to be assigned to ``setup.L``.

.. admonition:: Example --- Three Term Example

	.. math::

		L = \int_{t_0}^{t_f} \left( -2x_1 + u_1^2 + \sin(t) \right) \mathrm{d}t, \quad n_{x} = 1, n_{u} = 1

	.. code-block:: Matlab

		L(1).left = 0; L(1).right = 2; L(1).matrix = −2;
		L(2).left = 1; L(2).right = 1; L(2).matrix = 1;
		L(3).left = 0; L(3).right = 0; L(3).matrix = {@(t)sin(t)};
		setup.L = L;

.. admonition:: Example --- General Constant Quadratic State and Control Terms

	.. math::

		L = \int_{t_0}^{t_f} \left( \mathbf{x}^T \mathbf{Q} \mathbf{x} + \mathbf{u}^T \mathbf{R} \mathbf{u} \right) \mathrm{d}t

	.. code-block:: Matlab

		L(1).left = 2; L(1).right = 2; L(1).matrix = Q;
		L(2).left = 1; L(2).right = 1; L(2).matrix = R;
		setup.L = L;

.. admonition:: Example --- Time-varying Quadratic Control Term

	.. math::

		L_k = \int_{t_0}^{t_f} \sin(t) u_1^2 \mathrm{d}t, \quad n_u = 1

	.. code-block:: Matlab

		L(k).left = 1; L(k).right = 1; L(k).matrix = {@(t)sin(t)};

.. admonition:: Example --- Mixed Quadratic Term

	.. math::

		L_k = \int_{t_0}^{t_f} x_1(t_f) x_2(t) \mathrm{d}t, \quad n_{x} = 2

	.. code-block:: Matlab

		L(k).left = 4; L(k).right = 2; L(k).matrix = [0,0;1,0];

.. admonition:: Example --- Time-varying Linear State Term

	.. math::

		L_k = \int_{t_0}^{t_f} \left( e^{-t} x_1(t) + x_2(t) \right) \mathrm{d}t, \quad n_{x} = 3

	.. code-block:: Matlab

		L(k).left = 0; L(k).right = 2; L(k).matrix = {@(t)exp(−t), 1, 0};

Nonlinear Lagrange Terms
-------------------------------