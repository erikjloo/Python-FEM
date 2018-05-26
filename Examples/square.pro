 {
    "input":
    {
        "modules" : ["mesh", "hardening", "loads"],

        "mesh" :
        {
            "type" : "Gmsh",
            "file" : "Examples/square.msh",
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
        "models" : [ "matrix", "pbc" ],

        "matrix" :
        {
            "type"     : "Solid",
            "elements" : "1",
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

        "pbc" :
        {
            "type"     : "Periodic",
            "coarsenFactor" : 0.8,
            "strainRate" : [0.0001, -0.0002, 0.003]
        }
    },
    "solver":
    {
        "type" : "numpy"
    }
 }