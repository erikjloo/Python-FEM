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
                "type" : "PlaneStrain",
                "young" : 5000,
                "poisson" : 0.3
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
            "strainRate" : [0.0001, -0.0002, 0.003],
            "shape" :
            {
                "type" : "Line2",
                "scheme" : "Gauss"
            }
        }
    },
    "nonlin":
    {
        "type" : "full",
        "niter" : 1,
        "tol" : 1e-5,
        "solver" :
        {
            "type" : "lstsq"
        }
    }
 }