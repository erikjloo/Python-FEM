 {
    "Input":
    {

        "mesh" :
        {
            "type" : "Gmsh",
            "file" : "Examples/square.msh",
            "doElemGroups" : false
        },

        "loads" :
        {
            "type" : "Input",
            "file" : "load.data"
        }
    },

    "model" :
    {
        "type"     : "Solid",
        "elements" : "gmsh0",

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
            "type" : "Tri3",
            "scheme" : "Gauss"
        }
    }
 }