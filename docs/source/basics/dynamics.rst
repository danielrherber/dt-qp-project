Dynamic Constraints
===================

.. warning:: Currently, you can only have either linear dynamics or nonlinear dynamics included in ``setup``, not both.

Linear Dynamic Constraints
-------------------------------

The explicit first-order differential equation is in the following form:

.. math::

	\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u} + \mathbf{G}\mathbf{p} + \mathbf{d}

with states :math:`\mathbf{x}\in\mathbb{R}^{n_x}`, controls :math:`\mathbf{u}\in\mathbb{R}^{n_u}`, and parameters :math:`\mathbf{p}\in\mathbb{R}^{n_p}`, state matrix :math:`\mathbf{A}\in\mathbb{R}^{n_x \times n_x}`, input matrix :math:`\mathbf{B}\in\mathbb{R}^{n_x \times n_u}`, parameter matrix :math:`\mathbf{G}\in\mathbb{R}^{n_x \times n_p}`, and disturbances :math:`\mathbf{d}\in\mathbb{R}^{n_x \times 1}`.
Each of the matrices :math:`(\mathbf{A}, \mathbf{B}, \mathbf{G}, \mathbf{d})` has an associated field in ``setup``.
All of the matrices can be time-varying and can be omitted if they are not present.

.. admonition:: Example --- Complete Example

	.. math::

		\dot{\mathbf{x}} = \begin{bmatrix} 0 & 1 \\ -3 & -2 \end{bmatrix} \mathbf{x} + \begin{bmatrix} 1 \\ 2 \end{bmatrix} \mathbf{u} + \begin{bmatrix} 0 \\ \sin(t) \end{bmatrix}

	.. code-block:: Matlab

		setup.A = [0,1;−3,−2];
		setup.B = [1;2];
		setup.d = {0;@(t)sin(t)};

.. admonition:: Example --- Time-varying State Matrix

	.. math::

		\mathbf{A} = \begin{bmatrix} \sin(t) & 1 \\ 0 & e^{-t} \end{bmatrix}

	.. code-block:: Matlab

		setup.A = {@(t)sin(t), 1; 0, @(t)exp(−t)};

.. admonition:: Example --- Time-varying Input Matrix

	.. math::

		\mathbf{B} = \begin{bmatrix} \cos(t) \\ \sin(t) \end{bmatrix}

	.. code-block:: Matlab

		setup.B = {@(t)cos(t); @(t)sin(t)};

.. admonition:: Example --- Constant Parameter Matrix

	.. math::

		\mathbf{G} = \begin{bmatrix} 1 \end{bmatrix}

	.. code-block:: Matlab

		setup.G = 1;

.. admonition:: Example --- Time-varying Disturbance Matrix

	.. math::

		\mathbf{d} = \begin{bmatrix} 1 \\ 0 \\ \sin(t) \end{bmatrix}

	.. code-block:: Matlab

		setup.d = {1;0; @(t)sin(t)};


Nonlinear Dynamic Constraints
-------------------------------