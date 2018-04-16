# Import Standard Libraries
import os
import scipy as np
from configparser import ConfigParser

# Import Local Libraries
from dofspace import DofSpace
from solidModel import SolidModel
from constraints import Constraints
from properties import Properties, ElementType
from algebra import MatrixBuilder


# userInput =
# {
#     modules = ["mesh", "pbcgroups", "hardening", "loads"]

#     mesh =
#     {
#         type = "GmshInput"
#         file = "rve.msh"
#         doElemGroups = true
#     }

#     pbcgroups =
#     {
#         type = "PBCGroupInput"
#     }

#     hardening =
#     {
#         type = "Input"
#         file = "mises.data"
#     }

#     loads =
#     {
#         type = "Input"
#         file = "load.data"
#     }
# }

#   model       =
#   {
#     type   = "Multi";
#     models = [ "matrix", "fibers", "arclen" ];

#     matrix =
#     {
#       type     = "Solid";
#       elements = "gmsh1";

#       material =
#       {
#         type   = "Melro";
#         dim    = 2;
#         state  = "PLANE_STRAIN";   
        
#         young      = 3760.;
#         poisson    = 0.3;
#         poissonP   = 0.39;
#         rmTolerance = 1.e-10;
#         sigmaT = "st(x)";
#         sigmaC = "st(x)";
#       };
      
#       shape.type = "Triangle3";
#       shape.intScheme = "Gauss1";
#     };

#     fibers = 
#     {
#       type = "Solid";
#       elements = "gmsh0";

#       material = 
#       {
#         type = "Hooke";
#         dim = 2;
#         state = "PLANE_STRAIN";

#         young = 74000.;
#         poisson = 0.2;
#       };

#       shape.type = "Triangle3";
#       shape.intScheme = "Gauss1";
#     };

#     arclen = 
#     {
#       type = "DispArclen";
#      };

#       constraints = 
#       {
#         nodeGroups = [ "cornery","cornerx" ];
#         dofs = [ "dx","dy" ];
#         loaded = -1;
#         loads.loadTable = "load";
#       };
#     };
#   };
# };


# properties props()
# props.parseFile('input.pro')


# output:
# u, strain, sigma, f_int, K


