 {
    "input":
    {

        "mesh" :
        {
            "type" : "Gmsh",
            "file" : "Examples/2D_semicircle.msh",
            "rank" : 2,
            "doElemGroups" : false
        },

        "loads" :
        {
            "type" : "Input",
            "file" : "load.data"
        },

        "constraints" :
        {
            "type" : "Input",
            "file" : "constraints.data"
        }
    },

    "model" :
    {
        "type"     : "Solid",
        "elements" : "gmsh0",

        "material" :
        {
            "type"   : "PlaneStress",
            "young"  : 100000,
            "poisson": 0.2
        },
            
        "shape" :
        {
            "type" : "Tri3",
            "scheme" : "Gauss1"
        }
    }
 }