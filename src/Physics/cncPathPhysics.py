# Copyright (C) 2026  Cody Raskin

from PYB11Generator import *
from physics import *


@PYB11template("dim")
class CNCPathPhysics(Physics):

    def pyinit(self,
               nodeList="NodeList*",
               constants="PhysicalConstants&"):
        """CNC kinematic toolpath physics package.

        Parameters
        ----------
        nodeList : NodeList*
            NodeList containing the CNC tool node(s). Typically just 1 node.
        constants : PhysicalConstants&
            Physical constants object (passed through for consistency with other physics).
        """
        return

    def addLinearMove(self,
                      endPosition="const Lin::Vector<%(dim)s>&",
                      feed="double",
                      rapid="bool"):
        """Add a straight-line (G0/G1-style) move to the CNC path.

        Parameters
        ----------
        endPosition : Lin::Vector<dim> (const ref)
            Absolute end position of this move in machine coordinates.
        feed : double
            Linear speed along this segment (same units as your system, e.g. mm/s).
        rapid : bool
            True if this is a rapid move (G0). Currently treated the same as feed
            moves in the C++ implementation, but included for future differentiation.
        """
        return "void"

    def clearPath(self):
        """Clear all programmed CNC moves."""
        return "void"

    def pathComplete(self):
        """Return True if all programmed moves have been completed."""
        return "bool"


# Concrete template instantiations
CNCPathPhysics1d = PYB11TemplateClass(
    CNCPathPhysics,
    template_parameters=("1",),
    cppname="CNCPathPhysics<1>",
    pyname="CNCPathPhysics1d",
    docext=" (1D)."
)

CNCPathPhysics2d = PYB11TemplateClass(
    CNCPathPhysics,
    template_parameters=("2",),
    cppname="CNCPathPhysics<2>",
    pyname="CNCPathPhysics2d",
    docext=" (2D)."
)

CNCPathPhysics3d = PYB11TemplateClass(
    CNCPathPhysics,
    template_parameters=("3",),
    cppname="CNCPathPhysics<3>",
    pyname="CNCPathPhysics3d",
    docext=" (3D)."
)
