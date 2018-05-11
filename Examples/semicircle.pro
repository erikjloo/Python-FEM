 {
    "input":
    {

        "mesh" :
        {
            "type" : "Gmsh",
            "file" : "Examples/semicircle.msh",
            "rank" : 3,
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
            "type"   : "Hooke",
            "young"  : 100000,
            "poisson": 0.2
        },
            
        "shape" :
        {
            "type" : "Tetra4",
            "scheme" : "Gauss"
        }
    }
 }