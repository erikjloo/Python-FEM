 {
    "input":
    {
        "modules" : ["mesh", "loads", "constraints" ],

        "mesh" :
        {
            "type" : "XML",
            "file" : "Examples/uniaxial.xml",
            "rank" : 2,
            "doElemGroups" : false
        },

        "loads" :
        {
            "type" : "Input",
            "file" : "Examples/uniaxial.xml"
        },

        "constraints" :
        {
            "type" : "Input",
            "file" : "Examples/uniaxial.xml"
        }
    },

    "model" :
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
    "linsolve":
    {
        "solver" :
        {
            "type" : "lstsq"
        }
    }
 }