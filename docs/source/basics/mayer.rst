Mayer Objective Terms
============

Linear-Quadratic Mayer Terms
-------------------------------

Mayer terms are formed using the objective structure format.
There can be an arbitrary :math:`n_M` individual structures that are then summed to construct the complete Mayer term that is to be minimized:

.. math::

	M = \sum_{k=1}^{n_M} M_k = \sum_{k=1}^{n_M} \left[ \mathbf{z}_{\mathrm{left}}^T [\mathrm{matrix}] \mathbf{z}_{\mathrm{right}} \right]_k

The properly-constructed structures then need to be assigned to ``setup.M``.

.. warning:: ``left``/``right`` should represent time-independent variables (i.e., they should not be valued at 1 and 2).

.. admonition:: Example --- Complete Example

	.. math::

		M_k = \left(x_1(t_f)\right)^2 - 2x_1(t_0) + 1, \quad n_x = 1

	.. code-block:: Matlab

		M(1).left = 5; M(1).right = 5; M(1).matrix = 1;
		M(2).left = 0; M(2).right = 4; M(2).matrix = -2;
		M(3).left = 0; M(3).right = 0; M(3).matrix = 1;
		setup.M = M;

.. admonition:: Example --- Linear Parameter

	.. math::

		M_k = p_1, \quad n_p = 1

	.. code-block:: Matlab

		M(k).left = 3; M(k).right = 3; M(k).matrix = 1;

.. admonition:: Example --- Quadratic Initial States

	.. math::

		M_k = \left(x_1(t_0)\right)^2 + \left(x_2(t_0)\right)^2, \quad n_x = 2

	.. code-block:: Matlab

		M(k).left = 4; M(k).right = 4; M(k).matrix = eye(2);

.. admonition:: Example --- Linear Final States

	.. math::

		M_k = x_1(t_f) - 2x_2(t_f), \quad n_x = 2

	.. code-block:: Matlab

		M(k).left = 0; M(k).right = 5; M(k).matrix = [1;−2];



Nonlinear Mayer Terms
-------------------------------