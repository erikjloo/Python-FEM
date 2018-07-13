 {
    "input":
    {
        "modules" : ["mesh", "load", "cons" ],

        "mesh" :
        {
            "type" : "XML",
            "file" : "Examples/uniaxial.xml",
            "rank" : 2,
            "doElemGroups" : false
        },

        "load" :
        {
            "type" : "Loads",
            "file" : "Examples/uniaxial.xml"
        },

        "cons" :
        {
            "type" : "Constraints",
            "file" : "Examples/uniaxial.xml"
        }
    },

    "model" :
    {
        "type" : "Multi",
        "models" : ["model","load","cons"],
        "model":
        {
            "type"     : "Solid",
            "thickness" : 1,

            "material" :
            {
                "type"   : "PlaneStrain",
                "young"      : 3760,
                "poisson"    : 0.3
            },        

            "shape" :
            {
                "type" : "Quad4",
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