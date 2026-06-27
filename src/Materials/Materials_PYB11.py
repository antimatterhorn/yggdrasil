from PYB11Generator import *
PYB11includes = ['"material.hh"',
                 '"linearElastic.cc"']

PlaneCondition = PYB11enum(
    name="PlaneCondition",
    values=[
        ("Stress", "PlaneCondition::Stress"),
        ("Strain", "PlaneCondition::Strain"),
    ]
)

from material import *
from linearElastic import *
