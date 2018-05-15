 {
    "input":
    {
        "modules" : ["mesh", "hardening", "loads"],

        "mesh" :
        {
            "type" : "Gmsh",
            "file" : "Examples/rve.msh",
            "rank" : 2,
            "doElemGroups" : true
        },

        "hardening" :
        {
            "type" : "Input",
            "file" : "mises.data"
        },

        "loads" :
        {
            "type" : "Input",
            "file" : "load.data"
        }
    },

    "model" :
    {
        "type"   : "Multi",
        "models" : [ "matrix", "fibers", "pbc" ],

        "matrix" :
        {
            "type"     : "Solid",
            "elements" : "gmsh1",
            "thickness" : 1,

            "material" :
            {
                "type"   : "Melro",
                "rank"    : 2,
                "state"  : "PLANE_STRAIN",
                "young"      : 3760,
                "poisson"    : 0.3,
                "poissonP"   : 0.39,
                "rmTolerance" : 0.0000000001,
                "sigmaT" : "st(x)",
                "sigmaC" : "st(x)"
            },
            
            "shape" :
            {
                "type" : "Tri3",
                "scheme" : "Gauss"
            }
        },

        "fibers" : 
        {
            "type" : "Solid",
            "elements" : "gmsh0",
            "thickness" : 1,

            "material" : 
            {
                "type" : "PlaneStrain",
                "young" : 74000,
                "poisson" : 0.2
            },

            "shape" :
            {
                "type" : "Tri3",
                "scheme" : "Gauss"
            }
        },
        
        "pbc" :
        {
            "type"     : "Periodic",
            "coarsenFactor" : 0.3,
            "imposedStrain" : [0.1, 0.2, 0.3]
        }
    }
 }