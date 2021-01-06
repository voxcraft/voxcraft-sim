VXA Tags
========

.. _vxa_tags:

.. note:: Here is a sample base.vxa file.

VXA.GPU
-------

Configuration related to GPU settings.

VXA.GPU.HeapSize
^^^^^^^^^^^^^^^^

The proportion of memories in GPUs that we want to use as Heap Memory. 
If we use the GPU both for display and computing, like on my laptop, we can choose `0.5`; if we use the GPU sololy for computing, like on a GPU server, we can use `1.0`.

CUDA manages all memories in GPUs. The default heap size is 8M bytes, which is too small if we want to use `malloc` in CUDA code.

`See CUDA Documentation <https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__DEVICE.html#group__CUDART__DEVICE_1g05956f16eaa47ef3a4efee84563ccb7d>`_

VXA.Simulator
-------------

Configuration related to the whole simulation settings.

VXA.Simulator.FitnessFunction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a MathTree (refer to MathTree)

The fitness function for a creature. The simulator will sort all creatures in one generation by their fitness scores and report them as an output file. (refer to Voxelyze3 parameter -o)

VXA.Simulator.Integration.DtFrac
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0.0 ~ 1.0

Simulator will compute a recommended time step for each simulation, but sometimes the simulation still exploded even using the recommended time step.

So we can define a fraction here to say that we want to use even smaller time step to make sure the simulation runs safely.

`1.0` if everything is fine, `0.5` if we want a safer run.


VXA.Simulator.Damping.BondDampingZ
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0.0 ~ 1.0

This is a name used in VX1. In VX2 it's called `InternalDamping` or `zetaInternal`.

When we compute the force for the link, we multiply the calculated force by a `dampingMultiplier()`, which is related to this internal damping variable.

In another words, if we want the voxels to jiggle a lot, we use `0.0`, if we want the voxels to calm down, we use `1.0`.

VXA.Simulator.Damping.ColDampingZ
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0.0 ~ 1.0

This is a name used in VX1. In VX2 it's called `CollisionDamping` or `zetaCollision`.

`0.0` if we want the voxel bouncing on the floor, `1.0` if we don't want it bouncing.

.. note:: **Why the robot still bounce even when this variable is set to zero?**

    Because the links between vertical voxels will produce a bouncing force. 
    If there's only one voxel, the bounce will be determined by this variable alone.

VXA.Simulator.Damping.SlowDampingZ
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0.0 ~ 1.0

This is a name used in VX1. In VX2 it's called `GlobalDamping` or `zeta`.

`0.0` is no global damping, `1.0` to damp everything in the simulation.

.. warning:: This is a strong damping term. The whole simulation will almost freeze if this variable is set to `0.1`.

VXA.Simulator.StopCondition.StopConditionFormula
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a MathTree (refer to MathTree)

Tell the simulation when to stop. We usually stop at t-4>0, which means stop at 4.0 sec.

.. code-block:: XML

    <mtSUB>
        <mtVAR>t</mtVAR>
        <mtCONST>4</mtCONST>
    </mtSUB>


VXA.Simulator.RecordHistory.RecordStepSize
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: an integer

`0` if we don't want record any `.history` file. (`.history` file can be played in VoxCAD)

For example, `100` if we want to record the history every 100 steps.

.. warning:: Recording history will cost much longer time! Now it is about 40x.

VXA.Simulator.RecordHistory.RecordVoxel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`0` if we don't want to record voxels, `1` if we do.

VXA.Simulator.RecordHistory.RecordLink
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`0` if we don't want to record links, `1` if we do.

Recording links is usually for debugging porposes. By setting the alpha value of the materials less than 1, we can see the links in the playback.

VXA.Simulator.AttachDetach.EnableCollision
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

This variable is controling the voxel-voxel collision. Not related to the collsion with floor.

`0` if we don't want voxel-voxel collision, `1` if we do.

.. note:: Collision detection takes `O(n^2)` time, so disable this feature can make simulation much faster.

VXA.Simulator.AttachDetach.EnableAttach
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When collision happens, we can enable attachment. Under certain condition (defined in `AttachCondition`), two voxels will stick together when collide.

VXA.Simulator.AttachDetach.AttachCondition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a set of MathTrees (refer to MathTree)

If we want attachment happens whenever collision happens, we can define `Condition_0`, `Condition_1`, up to `Condition_4`.

.. code-block:: XML

    <Condition_0>
        <mtCONST>1</mtCONST>
    </Condition_0>

VXA.Simulator.AttachDetach.SafetyGuard
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: an integer

when attachment happens, there will be a new link formed between two voxels. Sometimes the relative speed of two voxels is too large, the attachment will seem to be unrealistic.

This variable defines the number of steps in which there will be a special damping between two newly attached voxels.

.. note:: This is the number of steps, not in seconds, so it will change if step size changes.

VXA.Simulator.ForceField
^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: three MathTrees for x,y,z dimension (refer to MathTree)

If we want to apply an external force field to the simulation, we can define it here. We can define `x_forcefield`, `y_forcefield`, and `z_forcefield`.

Here is an example to define a force field that only has value on x direction.

`x_forcefield = 0-x`, which means everything will be pull to y axis.

.. code-block:: XML

    <x_forcefield>
        <mtSUB>
            <mtCONST>0</mtCONST>
            <mtCONST>x</mtCONST>
        </mtSUB>
    </x_forcefield>
    <y_forcefield>
        <mtCONST>0</mtCONST>
    </y_forcefield>

VXA.Simulator.EnableSignals
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`0` if we want to disable singals. `1` if we want to enable the singals.


VXA.Simulator.EnableCilia
^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`0` if we want to disable cilia. `1` if we want to enable the cilia.

VXA.Simulator.SavePositionOfAllVoxels
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`0` if we don't want the output report XML file contains the final positions of all voxels, `1` if we do.

VXA.Simulator.MaxDistInVoxelLengthsToCountAsPair
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a real number with no unit

Sometimes we need to count how many pairs of TARGET voxels are close to each other. By defining this variable, we can specify the threshold for counting.

`0` if we don't want to count close pairs.

.. note:: This quantity is the distance over average voxel length. For example, if the voxel length is 0.01 meter, then if we set this variable to 2 here, it means the distance is 2*0.01 meter.

**Value**: 0 or 1

`0` if we don't want to enable counting closeness

VXA.Environment
---------------

Configuration related to the whole virtual environment.

VXA.Environment.Gravity.GravEnabled
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`0` if we don't want gravity, `1` if we want gravity.

.. note:: If we disable the gravity here, we can still use the Force Field to define a downward force that is identical to gravity.

    However, if we enable force field will be a little bit slower than simply using this gravity variable.
    (refer to force field)

VXA.Environment.Gravity.GravAcc
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a real number in m/s^2

`-9.81` if we want to use the common gravity on earth. Negative means downward.

VXA.Environment.Gravity.FloorEnabled
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`1` if we want there to be a floor at `z=0`, so that thing won't fall forever.

.. note:: This variable is not related to whether to draw the floor in VoxCAD. That can be controled via a checkbox in VoxCAD.

VXA.Environment.Thermal.VaryTempEnabled
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`1` if we want enable temperature.

.. note:: **Why temperature?**

    One way to make the voxel actuate and let the robot move is to use a varying temperature.
    In this thermal expansion model, the rest length of the links between voxels will vary due to the temperature change.
    We can define different coefficient of thermal expansion (CTE) for different materials. 
    (refer to CTE)

VXA.Environment.Thermal.TempAmplitude
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a real number in degree Celsius

The amplitude of the temperature oscillation. If we want the voxels to actuate more, we can use larger amplitude.


VXA.Environment.Thermal.TempPeriod
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: Real number in second

The period of the temperature oscillation. If we want the actuate period to be longer, we can use larger period here.

.. note:: The temperature will oscillate as a sine wave.


VXA.VXC
-------

Configuration related to the creatures in the simulation.

VXA.VXC.Lattice.Lattice_Dim
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a real number in meter

The dimension (side length) of a voxel.

VXA.VXC.Palette
^^^^^^^^^^^^^^^

Palette means the materials.

VXA.VXC.Palette.Material
^^^^^^^^^^^^^^^^^^^^^^^^

One material.

VXA.VXC.Palette.Material.Name
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a string

The name of the material.

VXA.VXC.Palette.Material.Display
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The color of the material.

VXA.VXC.Palette.Material.Display.Red
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0.0 ~ 1.0

VXA.VXC.Palette.Material.Display.Green
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0.0 ~ 1.0

VXA.VXC.Palette.Material.Display.Blue
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0.0 ~ 1.0

VXA.VXC.Palette.Material.Display.Alpha
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0.0 ~ 1.0

`See RGBA color model <https://en.wikipedia.org/wiki/RGBA_color_model>`_

VXA.VXC.Palette.Material.Mechanical.isTarget
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`1` if we want voxels made of this material to be the target, `0` if we don't.

Target voxels will trigger special functionalities. For example, when a non-target voxel hit a target voxel, the former one will generate a signal; or, when we use `MeasureFitnessOfTargetMaterialOnly` tag, the fitness function will only take into account the target voxels instead of all voxels.

VXA.VXC.Palette.Material.Mechanical.isMeasured
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1 (default 1)

`1` if we want to measure voxels made by this material in all MathTree functions, especially the fitness function. `0` if we want to exclude this material.


VXA.VXC.Palette.Material.Mechanical.Fixed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`1` if we don't want this material to move at all.

The fixed voxels can serve as the environment, such as wall, steps, etc., or serve as a pin when we want to pin the robot down in space.

VXA.VXC.Palette.Material.Mechanical.sticky
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`1` if we want attachment can happen to this material, `0` if we don't.

VXA.VXC.Palette.Material.Mechanical.Cilia
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a real number

`1` if we want to enable cilia force for this material, `0` if we don't.

Other real number if you want to change the magnitude of the cilia force.

VXA.VXC.Palette.Material.Mechanical.isPaceMaker
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`1` if this material can generate periodic signals spontaneously, `0` if not.

VXA.VXC.Palette.Material.Mechanical.PaceMakerPeriod
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a real number in second

The pace maker can generate periodic signals spontaneously. The period between two signals is defined by this variable.

VXA.VXC.Palette.Material.Mechanical.signalValueDecay
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0.0 ~ 1.0

When the singal propagates to other part of the body, it has a decay ratio. This variable defines the ratio.

`0.0` means the signal cannot propagate at all, `1.0` means the signal never decay and can propagate to infinity.

VXA.VXC.Palette.Material.Mechanical.signalTimeDelay
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a real number in second

When the singal propagates to other part of the body, it has a travel speed. The signal may be delayed at every stop (in every voxel). This variable defines how much time it will delay in each voxel.

VXA.VXC.Palette.Material.Mechanical.inactivePeriod
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a real number in second

Inspired by the process of action potential in living cells, in which the cell will enter an inactive state for a while to prevent the signal traveling backward.
This variable defines the time period that a voxel stays inactive after sending out the signal.

`See Action Potential <https://en.wikipedia.org/wiki/Action_potential>`_

VXA.VXC.Palette.Material.Mechanical.MatModel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`0` for simple linear elastic model. The mechanical model for a perfectly elastic material is a simple spring.
`1` for linear elastic model that can fail. When materials are subjected to a large enough strain they fail by fracture. The voxels will detach.
(refer to AttachDetach)

.. note:: Fail by fracture model need `VXA.Simulator.AttachDetach.EnableDetach` to be 1.

VXA.VXC.Palette.Material.Mechanical.Elastic_Mod
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a real number in Pascal

The elastic modulus (a.k.a. Young's Modulus) describes the stiffness of a material. For soft robotics, we usually use 10^7 Pa like a rubber.

`See values for common materials <https://en.wikipedia.org/wiki/Young%27s_modulus#Approximate_values>`_

VXA.VXC.Palette.Material.Mechanical.Fail_Stress
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a real number in Pascal

When the stress is larger than this threshold, the material fail by fracture.

.. note:: This need `MatModel` to be `1` and `EnableDetach` to be 1.

VXA.VXC.Palette.Material.Mechanical.Density
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a real number in kg/m^3

For example, natural rubbber's density is about 1.5e+3 kg/m^3.

VXA.VXC.Palette.Material.Mechanical.Poissons_Ratio
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0.0 ~ 0.5

The ratio of the proportional decrease in a lateral measurement to the proportional increase in length in a sample of material that is elastically stretched.

For example, rubber has a ratio near 0.5, and cork is famous for a ratio near 0.0.

VXA.VXC.Palette.Material.Mechanical.CTE
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a small real number in 1/degree Celsius ( same as 1/K )

For example, plastics has CTE of about 10^(-4). To make the actuation more obvious, we can choose CTE to be 0.01. (in reality, we don't have such high CTE material.)



VXA.VXC.Palette.Material.Mechanical.uStatic
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0.0 ~ 5.0

This is a name in VX1. In VX2 it's called `StaticFriction`.

This is the static frictional coefficient. For example, rubber-rubber static friction coefficient is 1.16, and ice-ice static friction coefficient is 0.1.


VXA.VXC.Palette.Material.Mechanical.uDynamic
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0.0 ~ 1.0

This is a name in VX1. In VX2 it's called `KineticFriction`.

This is the kinetic frictional coefficient. For example, rubber-pavement kinetic friction coefficient is 0.8, and steel-ice kinetic friction coefficient is 0.01.

VXA.VXC.Palette.Material.Mechanical.Cilia
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: 0 or 1

`0` if this material does not exerting cilia force, `1` if it has.

VXA.VXC.Structure
^^^^^^^^^^^^^^^^^

The simulation start with a world of voxels, which are in a lattice structure.

We should always use the attribute `Compression="ASCII_READABLE"`, since by doing that we can see the data directly.

VXA.VXC.Structure.X_Voxels
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: an integer larger than zero

The x dimension of the voxel world.

VXA.VXC.Structure.Y_Voxels
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: an integer larger than zero

The y dimension of the voxel world.

VXA.VXC.Structure.Z_Voxels
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: an integer larger than zero

The z dimension of the voxel world.

VXA.VXC.Structure.Data
^^^^^^^^^^^^^^^^^^^^^^

This section defines material type for each position in the lattice.

The x and y dimension combined together form a layer. The first layer is z=1, the second layer is z=2, aranged from bottom up.

VXA.VXC.Structure.Data.Layer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a string of integers, which define the material for each voxel in one layer. 

The number is the index of material but in charactor type, `0` for nothing. For example, the first material corresponds to '1', the nine-th material corresponds to '9'.

If we have more than 9 materials, we can continue using ':',';','<'... following '9' which is in the ASCII order.

`See ASCII order here <https://www.ascii-code.com/>`_

The length of digits in one layer should equal X_Voxels * Y_Voxels.

The number of layers should equal Z_Voxels.

VXA.VXC.Structure.PhaseOffset
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One way to make the voxel actuate is using thermal expansion model, with `VaryTempEnabled` and `CTE`, the rest length of voxels will change due to temperature changing.

However, if all voxels actuate in phase, it is quite difficult to find any interesting behaviors.

Phase offset introduce more interesting behavior by allowing each voxel has its own phase.

Phase offset settings should have the same dimension as `Data`, and each value corresponds to the phase offset of that voxel.

VXA.VXC.Structure.PhaseOffset.Layer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a set of real number 0.0 ~ 1.0

VXA.VXC.Structure.BaseCiliaForce
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a set of real numbers, which define the cilia force in x,y,z dimension for each voxel in one layer.

.. note:: This feature works with `EnableCilia` = 1 and only apply to material with `Cilia` = 1.

VXA.VXC.Structure.ShiftCiliaForce
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Value**: a set of real numbers, which define the behavior shifting of the cilia force in x,y,z dimension for each voxel in one layer.

When a voxel has a signal larger than 0, there will be a shifting in behavior.

VXA.RawPrint
------------

**Value**: a string

What was passed here will be simply passed along to the history file (or standard output).
