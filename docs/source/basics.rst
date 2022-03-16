*************
Problem Setup
*************

There are three primary optimization variable groupings:

.. math::

   \mathbf{z} = \begin{bmatrix} \mathbf{u}(t) \\ \mathbf{x}(t) \\ \mathbf{p}(t) \end{bmatrix}

with states :math:`\mathbf{x}(t)\in\mathbb{R}^{n_x}`, controls :math:`\mathbf{u}(t)\in\mathbb{R}^{n_u}`, and parameters :math:`\mathbf{p}(t)\in\mathbb{R}^{n_p}`.
This list can be extended with the pseudo-variables for the initial state values :math:`\mathbf{x}_0\in\mathbb{R}^{n_x}` and final state values :math:`\mathbf{x}_f\in\mathbb{R}^{n_x}` often used directly in the problem formulation and the unit vector:

.. math::

   \bar{\mathbf{z}} = \begin{bmatrix} \mathbf{1} \\ \mathbf{u} \\ \mathbf{x} \\ \mathbf{p} \\ \mathbf{x}_0 \\ \mathbf{x}_f \end{bmatrix}

These groupings are summarized in the following table with the number column corresponding to the value needed to reference a particular variable grouping:

.. list-table::
   :align: center
   :header-rows: 1

   * - Original Variable
     - :math:`\bar{\mathbf{z}}`
     - Number
   * - Unit vector :math:`\mathbf{1}`
     - :math:`\mathbf{z}_0`
     - 0
   * - Controls :math:`\mathbf{u}(t)`
     - :math:`\mathbf{z}_1(t)`
     - 1
   * - States :math:`\mathbf{x}(t)`
     - :math:`\mathbf{z}_2(t)`
     - 2
   * - Parameters :math:`\mathbf{p}`
     - :math:`\mathbf{z}_3`
     - 3
   * - Initial States :math:`\mathbf{x}(t_0)`
     - :math:`\mathbf{z}_4`
     - 4
   * - Final States :math:`\mathbf{x}(t_f)`
     - :math:`\mathbf{z}_5`
     - 5

For linear-quadratic dynamic optimization (:term:`LQDO`) problems, all problems can be written as follows:

.. math::

   \underset{\bar{\mathbf{z}}(t)}{\text{minimize:}} \quad & o = \int_{t_0}^{t_f} L(t,\bar{\mathbf{z}}(t)) \mathrm{d}t + M(\mathbf{p}, \mathbf{x}_0, \mathbf{x}_f) \\
   \text{subject to:} \quad & \dot{\mathbf{x}}(t) = \mathbf{A}(t)\mathbf{x}(t) + \mathbf{B}(t)\mathbf{u}(t) + \mathbf{G}(t)\mathbf{p}(t) + \mathbf{d}(t) \\
   & \mathbf{Y}^T(t) \bar{\mathbf{z}}(t) = \mathbf{y}_{0}(t) \\
   & \mathbf{Z}^T(t) \bar{\mathbf{z}}(t) \leq \mathbf{z}_{0}(t) \\
   & \bar{\mathbf{z}}(t) \leq \mathbf{z}_{\max}(t) \\
   & \bar{\mathbf{z}}(t) \geq \mathbf{z}_{\min}(t) \\
   \text{where:} \quad & \mathbf{x}_0 = \mathbf{x}(t_0), \mathbf{x}_f = \mathbf{x}(t_f)

with


.. list-table::
   :align: center
   :header-rows: 0

   * - Name
     - Expression/Equation
   * - :ref:`Linear-Quadratic Lagrange Terms`
     - :math:`\int_{t_0}^{t_f} L(t,\bar{\mathbf{z}}(t)) \mathrm{d}t`
   * - :ref:`Linear-Quadratic Mayer Terms`
     - :math:`M(\mathbf{p}, \mathbf{x}_0, \mathbf{x}_f)`
   * - :ref:`Linear Dynamic Constraints`
     - :math:`\dot{\mathbf{x}}(t) = \mathbf{A}(t)\mathbf{x}(t) + \mathbf{B}(t)\mathbf{u}(t) + \mathbf{G}(t)\mathbf{p}(t) + \mathbf{d}(t)`
   * - :ref:`Additional Linear Constraints` (Equality)
     - :math:`\mathbf{Y}^T(t) \bar{\mathbf{z}}(t) = \mathbf{y}_{0}(t)`
   * - :ref:`Additional Linear Constraints` (Inequality)
     - :math:`\mathbf{Z}^T(t) \bar{\mathbf{z}}(t) \leq \mathbf{z}_{0}(t)`
   * - :ref:`Simple Bound Constraints` (Upper)
     - :math:`\bar{\mathbf{z}}(t) \leq \mathbf{z}_{\max}(t)`
   * - :ref:`Simple Bound Constraints` (Lower)
     - :math:`\bar{\mathbf{z}}(t) \geq \mathbf{z}_{\min}(t)`


.. note:: Suggested future revision to notation:

   .. math::

      \underset{\bar{\mathbf{z}}(t)}{\text{minimize:}} \quad & o = \int_{t_0}^{t_f} L(t,\bar{\mathbf{z}}(t)) \mathrm{d}t + M(\mathbf{p}, \mathbf{x}_0, \mathbf{x}_f) \\
      \text{subject to:} \quad & \dot{\mathbf{x}}(t) = \mathbf{A}(t)\mathbf{x}(t) + \mathbf{B}(t)\mathbf{u}(t) + \mathbf{V}(t)\mathbf{p}(t) + \mathbf{d}(t) \\
      & \mathbf{H}^T(t) \bar{\mathbf{z}}(t) = \mathbf{h}(t) \\
      & \mathbf{G}^T(t) \bar{\mathbf{z}}(t) \leq \mathbf{g}(t) \\
      & \bar{\mathbf{z}}(t) \leq \mathbf{z}_{\max}(t) \\
      & \bar{\mathbf{z}}(t) \geq \mathbf{z}_{\min}(t) \\
      \text{where:} \quad & \mathbf{x}_0 = \mathbf{x}(t_0), \mathbf{x}_f = \mathbf{x}(t_f)

.. toctree::
   :hidden:

   /basics/lagrange
   /basics/mayer
   /basics/dynamics
   /basics/additional-constraints
   /basics/simple-bounds