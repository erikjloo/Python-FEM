 {
    "input":
    {
        "modules" : ["mesh","load","cons"],
        "mesh" :
        {
            "type" : "Gmsh",
            "file" : "Examples/2D_semicircle.msh",
            "rank" : 2,
            "doElemGroups" : false
        },

        "load" :
        {
            "type" : "Loads",
            "file" : "Examples/2D_semicircle.xml"
        },

        "cons" :
        {
            "type" : "Constraints",
            "file" : "Examples/2D_semicircle.xml"
        }
    },
    "model" :
    {
        "type" : "Multi",
        "models" : ["model","load","cons"],
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
        "load" : 
        {
            "type" : "PointLoad",
            "loadTable" : "load"
        },
        "cons" :
        {
            "type" : "Constraints",
            "conTable" : "cons"
        }
    },
    "linsolve":
    {
        "solver": { "type" : "lstsq" }
    }
 }