Class Guides
=======================
Most classes in Yggdrasil are purely virtual. This allows for the interaction between class species (like physics classes
and integrators) to be abstracted and predictable.  

Physics Classes
--------------------
The primary class that you'll interact with in Yggdrasil is the ``Physics`` class. This class contains the prescription
for the State vector and all of the necessary logic for computing derivatives of the State vector. The general concept
for how physics classes work in Yggdrasil is that the ``NodeList`` holds immutable quantities that are valid for the beginning
of an integration step, while the State vector (usually a subset copy of some of the Fields on the ``NodeList``) may update
throughout a step. 
Yggdrasil's integrators expect assigned (derived) physics classes to override the physics base class method ``EvaluateDerivatives``
at least, and some may also override ``PrestepInitialize`` and ``FinalizeStep`` as well. 

Let's use ``constantGravity.cc`` in the ``src/Physics`` directory as an example of a derived class that overrides only the derivatives.

.. literalinclude:: ../../src/Physics/constantGravity.cc
   :language: c++
   :linenos:
   :lines: 4-27
   :lineno-start: 4

In the constructor for ``ConstantGravity``, we pass in the ``NodeList`` pointer, a ``PhysicalConstants`` object, and a ``Vector`` that
prescribes the direction and magnitude of the acceleration. Then we enroll the ``acceleration``, ``velocity``, and ``position`` Vector 
Fields using the type templated EnrollFields method from the Physics base class. 
This ensures that these fields exist on the NodeList, and if they don't, this 
method will create them. 

Next we enroll the ``position`` and ``velocity`` Fields to the ``State``.
In the case of 
the constant gravity package, the only derivatives are the position derivative and velocity derivative, and as a standard, we keep
Fields whose derivatives we want to iterate throughout an integration step
on the ``State`` object. For that reason, we do not assign ``acceleration`` to the ``State`` 
object for this physics class since it never changes.

.. literalinclude:: ../../src/Physics/constantGravity.cc
   :language: c++
   :linenos:
   :lines: 29-57
   :lineno-start: 36

``EvaluateDerivatives`` is used to compute the derivatives of the ``State`` object at each node. In this case, we are computing the 
change in velocity and position due to a constant acceleration, so we have a pair of ODEs:

.. math::
    \begin{aligned}
    \dot{x} &= v(t) + a \Delta t\\
    \dot{v} &= a\\
    \end{aligned}

The integrator invokes this method with intial and derivatives ``State`` vectors and a time step, and the physics class stores the value of those
derivatives back to the ``deriv`` ``State`` object.

.. note::
    The integrator will call this method multiple times per time step if the particular integrator in question is of high temporal order. 
    Thus, the value of ``dt`` here may be a partial step and the initial ``State`` object may change within a time step, 
    but in the case of a simple, forward Euler integration (:math:`\mathcal{O}(\Delta t)`), ``dt=0`` means that we are 
    computing the derivatives using ``State`` vector quantities evaluated at the current time. 

The final bit of code in ``EvaluateDerivatives`` merely calculates the minimum timestep based on a velocity condition.

Generally, the base class method for ``FinalizeStep`` should suffice for most physics
packages as the base class takes care of copying the final State after an integration step back to the ``NodeList``
on the derived physics class.  

.. note::
   If you don't want to recapitulate what is already in the base class ``FinalizeStep`` method, but you'd still like to
   tie up some loose ends at the end of an integration step (*e.g.* you may want to set floor limits on certain
   quantities), you can do so by overriding the ``FinalChecks()`` method. 
   By dafault, this method doesn't do anything, but is always called after ``FinalizeStep``.

Integrator Classes
------------------
Nihoggr's integrators follow a basic pattern of initializing the state of the physics objects at each step, applying boundary conditions,
evaluating derivatives, advancing the state of each physics object from those derivatives, and then finalizing each physics object's state.
An example of this pattern for simple forward Euler integration is shown below:

.. literalinclude:: ../../src/Integrators/integrator.cc
   :language: c++
   :lines: 6-66

.. note::
   While this class is a foward Euler integrator, it is also the base class for all integrators in Yggdrasil.

Equations of State
------------------
Equations of state are possibly the simplest classes in Yggdrasil. They consume Field objects (or referenced doubles) 
and set the values of other Fields according to their respective closure equations. Most of the logic for how they 
work is self-described by the base class interface file ``equationOfState.hh``.

.. literalinclude:: ../../src/EOS/equationOfState.hh
   :language: c++
   :linenos:
   :lines: 9-51
   :lineno-start: 9

.. note::
    Equations of state are not typically templated in Yggdrasil since they act only on scalar Fields. 