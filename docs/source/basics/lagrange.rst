Lagrange Objective Terms
============


Linear-Quadratic Lagrange Terms
-------------------------------

Lagrange terms are formed using the objective structure.
There can be an arbitrary :math:`n_L` individual structures that are then summed to construct the complete Lagrange term:

.. math::

	L_k = \sum_{k=1}^{n_L} L_k

The properly-constructed structures then need to be assigned to ``setup.L``.

 .. admonition:: Example --- General Quadratic State Term

	.. math::

		L_k = \int_{t_0}^{t_f} \mathbf{\xi}^T \mathbf{Q} \mathbf{\xi} \mathrm{d}t

	.. code-block:: Matlab

		L(k).left = 2; L(k).right = 2; L(k).matrix = Q;


 .. admonition:: Example --- Time-varying Quadratic Control Term

	.. math::

		L_k = \int_{t_0}^{t_f} \sin(t) u_1^2 \mathrm{d}t, \quad n_u = 1

	.. code-block:: Matlab

		L(k).left = 1; L(k).right = 1; L(k).matrix = {@(t)sin(t)};


 .. admonition:: Example --- 

	.. math::

		L_k = \int_{t_0}^{t_f} \xi_1(t_f) \xi_2(t) \mathrm{d}t, \quad n_{\xi} = 2

	.. code-block:: Matlab

		L(k).left = 4; L(k).right = 2; L(k).matrix = [0,0;1,0];


 .. admonition:: Example --- Time-varying Linear State Term

	.. math::

		L_k = \int_{t_0}^{t_f} \left( e^{-t} \xi_1(t) + \xi_2(t) \right) \mathrm{d}t, \quad n_{\xi} = 3

	.. code-block:: Matlab

		L(k).left = 0; L(k).right = 2; L(k).matrix = {@(t)exp(−t), 1, 0};


 .. admonition:: Example --- Three Terms

	.. math::

		L = \int_{t_0}^{t_f} \left( -2\xi_1 + u_1^2 + \sin(t) \right) \mathrm{d}t, \quad n_{\xi} = 1, n_{u} = 1

	.. code-block:: Matlab

		L(1).left = 0; L(1).right = 2; L(1).matrix = −2;
		L(2).left = 1; L(2).right = 1; L(2).matrix = 1;
		L(3).left = 0; L(3).right = 0; L(3).matrix = {@(t)sin(t)};

Nonlinear Lagrange Terms
-------------------------------