Additional Constraints
======================

Additional Linear Constraints
-----------------------------

The constraint structure format is used to represent additional linear constraints.
Use ``Y`` and ``Z`` to account for equality and inequality terms, respectively.
The field ``linear`` can have multiple values to represent the summation needed for certain constraints.
The fields ``right`` and ``matrix`` are analogous to their use in the objective function terms.
The value for ``b`` is the potentially time-varying function.

.. list-table::
   :align: center
   :header-rows: 0

   * - Linear equality constraints ``Y``
     - :math:`\sum_{i} \left[\mathrm{matrix}\right]_i \left[\mathbf{z}_{\mathrm{right}} \right]_i = \mathrm{b}_k`
   * - Linear inequality constraints ``Z``
     - :math:`\sum_{i} \left[\mathrm{matrix}\right]_i \left[\mathbf{z}_{\mathrm{right}} \right]_i \leq \mathrm{b}_k`

The properly-constructed structures then need to be assigned to ``setup.Y`` and ``setup.Z``, respectively

.. tip:: Simple bounds can be constructed using additional linear constraints but should be implemented as described in :ref:`Simple Bound Constraints`.

.. admonition:: Example --- Control/Parameter Equality Constraint

	.. math::

		u_1 - p_1 = \sin(t), \quad n_u = 1, n_p = 1

	.. code-block:: Matlab

		Y(k).linear(1).right = 1; Y(k).linear(1).matrix = 1;
		Y(k).linear(2).right = 3; Y(k).linear(2).matrix = −1;
		Y(k).b = @(t)sin(t);

.. admonition:: Example --- Control/Parameter Inequality Constraint

	.. math::

		u_1 - p_1 \leq \sin(t), \quad n_u = 1, n_p = 1

	.. code-block:: Matlab

		Z(k).linear(1).right = 1; Z(k).linear(1).matrix = 1;
		Z(k).linear(2).right = 3; Z(k).linear(2).matrix = −1;
		Z(k).b = @(t)sin(t);

.. admonition:: Example --- Linear State Inequality Constraint

	.. math::

		x_1 - 2x_2 \leq 0, \quad n_x = 2

	.. code-block:: Matlab

		Z(k).linear(1).right = 2; Z(k).linear(1).matrix = [1;−2];
		Z(k).b = 0;

Additional Nonlinear Constraints
--------------------------------