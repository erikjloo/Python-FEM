 {
    "userInput":
    {
        "modules" : ["mesh", "pbcgroups", "hardening", "loads"],

        "mesh" :
        {
            "type" : "GmshInput",
            "file" : "Examples/rve.msh",
            "doElemGroups" : true
        },

        "pbcgroups" :
        {
            "type" : "PBCGroupInput"
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
        "models" : [ "matrix", "fibers", "arclen" ],

        "matrix" :
        {
            "type"     : "Solid",
            "elements" : "gmsh1",

            "material" :
            {
                "type"   : "Melro",
                "dim"    : 2,
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
                "type" : "Triangle3",
                "intScheme" : "Gauss1"
            }
        },

        "fibers" : 
        {
            "type" : "Solid",
            "elements" : "gmsh0",

            "material" : 
            {
                "type" : "Hooke",
                "dim" : 2,
                "state" : "PLANE_STRAIN",

                "young" : 74000,
                "poisson" : 0.2
            },

            "shape" :
            {
                "type" : "Triangle3",
                "intScheme" : "Gauss1"
            }
        }
    }
 }