 {
    "input":
    {

        "mesh" :
        {
            "type" : "Gmsh",
            "file" : "Examples/3D_semicircle.msh",
            "rank" : 3,
            "doElemGroups" : false
        },

        "loads" :
        {
            "type" : "Input",
            "file" : "Examples/3D_semicircle.xml"
        },

        "constraints" :
        {
            "type" : "Input",
            "file" : "Examples/3D_semicircle.xml"
        }
    },

    "model" :
    {
        "type"     : "Solid",

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
    },
    "linsolve":
    {
        "solver":
        {
            "type" : "lstsq"
        }
    }
 }