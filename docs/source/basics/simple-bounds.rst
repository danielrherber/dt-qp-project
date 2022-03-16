Simple Bound Constraints
========================

The bound structure format is used to represent additional linear constraints that can be written as simple upper and lower bounds.
Use ``UB`` and ``LB`` to account for upper and lower terms, respectively.
The fields ``right`` and ``matrix`` are used similarly as the other setup structures, with the exception that the values of matrix can be :math:`\pm\infty` to indicate no bounds when appropriate:

.. list-table::
   :align: center
   :header-rows: 0

   * - Upper bounds ``UB``
     - :math:`\left[ \mathbf{z}_{\mathrm{right}} \right]_k \leq \left[ \mathrm{matrix} \right]_k`
   * - Lower bounds ``LB``
     - :math:`\left[ \mathbf{z}_{\mathrm{right}} \right]_k \geq \left[ \mathrm{matrix} \right]_k`

The properly-constructed structures then need to be assigned to ``setup.UB`` and ``setup.LB``, respectively.

.. tip:: You can define multiple simple bounds on the same variables and the code with automatically include the strictest constraint at each point in time.

.. admonition:: Example --- Single State Upper Bound

	.. math::

		x_2 \leq \pi, \quad n_x = 2

	.. code-block:: Matlab

		UB(k).right = 2; UB(k).matrix = [inf;pi];

.. admonition:: Example --- Single State Lower Bound

	.. math::

		x_2(t_0) \geq 0 \quad n_x = 3

	.. code-block:: Matlab

		LB(k).right = 4; LB(k).matrix = [−inf;0;−inf];

.. admonition:: Example --- Zero Initial States

	.. math::

		\mathbf{x}(t_0) = \mathbf{0}

	.. code-block:: Matlab

		UB(k).right = 4; UB(k).matrix = zeros(1,nx);
		LB(k).right = 4; LB(k).matrix = zeros(1,nx);

.. admonition:: Example --- Time-varying Control Lower Bounds

	.. math::

		u_1 \geq \sin(t), u_2 \geq 0, \quad n_u = 2

	.. code-block:: Matlab

		LB(k).right = 1; LB(k).matrix = {@(t)sin(t);-inf};