# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
PYB11includes = ['"physics.hh"',
                '"kinematics.hh"',
                '"../Materials/material.hh"',
                '"../Materials/linearElastic.cc"',
                '"femLinearElasticity.cc"',
                '"constantForce.cc"',
                '"pointSourceGravity.cc"',
                '"nBodyGravity.cc"',
                '"waveEquation.cc"',
                '"hydro.hh"',
                '"simplePhysics.cc"',
                '"implicitPhysics.cc"',
                '"eulerHydro.cc"',
                '"gridHydroHLLE.cc"',
                '"gridHydroHLLC.cc"',
                '"gridHydroKT.cc"',
                '"kinetics.cc"',
                '"fem.cc"',
                '"thermalConduction.cc"',
                '"phaseCoupling.cc"',
                '"treeGravity.cc"',
                '"reactionDiffusion.cc"',
                '"complexWaveEquation.cc"',
                '"magField.cc"',
                '"cncPathPhysics.cc"',
                ]

from physics import *
from kinematics import *
from constantForce import *
from pointSourceGravity import *
from nBodyGravity import *
from waveEquation import *
from hydro import *
from simplePhysics import *
from implicitPhysics import *
from eulerHydro import *
from gridHydroHLLE import *
from gridHydroHLLC import *
from gridHydroKT import *
from kinetics import *
from fem import *
from femLinearElasticity import *
from thermalConduction import *
from phaseCoupling import *
from treeGravity import *
from reactionDiffusion import *
from complexWaveEquation import *
from magField import *
from cncPathPhysics import *