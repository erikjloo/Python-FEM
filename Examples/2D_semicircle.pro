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
            "file" : "Examples/2D_semicircle.xml"
        },

        "constraints" :
        {
            "type" : "Input",
            "file" : "Examples/2D_semicircle.xml"
        }
    },

    "model" :
    {
        "type"     : "Solid",
        "thickness" : 0.1,

        "material" :
        {
            "type"   : "PlaneStrain",
            "young"  : 100000,
            "poisson": 0.2
        },
            
        "shape" :
        {
            "type" : "Tri3",
            "scheme" : "Gauss1"
        }
    },
    "solver":
    {
        "type" : "numpy"
    }
 }