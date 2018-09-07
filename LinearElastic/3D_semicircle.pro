 {
    "input":
    {
        "modules" : ["mesh","load","cons"],
        "mesh" :
        {
            "type" : "Gmsh",
            "file" : "LinearElastic/3D_semicircle.msh",
            "rank" : 3,
            "doElemGroups" : false
        },

        "load" :
        {
            "type" : "Loads",
            "file" : "LinearElastic/3D_semicircle.xml"
        },

        "cons" :
        {
            "type" : "Constraints",
            "file" : "LinearElastic/3D_semicircle.xml"
        }
    },
    "model" :
    {
        "type" : "Multi",
        "models" : ["model","load","cons"],
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