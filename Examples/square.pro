 {
    "input":
    {
        "mesh" :
        {
            "type" : "Gmsh",
            "file" : "Examples/square.msh",
            "rank" : 2,
            "doElemGroups" : false
        },

        "load" :
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
            "thickness" : 1,

            "material" : 
            {
                "type" : "PlaneStrain",
                "young" : 500,
                "poisson" : 0.3
            },
            
            "shape" :
            {
                "type" : "Tri3"
            }
        },

        "pbc" :
        {
            "type"     : "Periodic",
            "coarsenFactor" : 1.0,
            "strainRate" : [0.0, 0.0, 0.005],
            "shape" :
            {
                "type" : "Line2"
            }
        }
    },
    "linsolve":
    {
        "solver" :
        {
            "type" : "lstsq"
        }
    },
    "nonlin":
    {   
        "niter" : 1,
        "solver" :
        {
            "type" : "lstsq"
        }
    }
 }