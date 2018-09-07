 {
    "control": { "nsteps" : 4 },
    "input":
    {
        "modules" : ["mesh","load","cons"],
        "mesh" :
        {
            "type" : "Gmsh",
            "file" : "Examples/square.msh",
            "rank" : 2,
            "doElemGroups" : false
        },
        "load" : 
        { 
            "type" : "Loads",
            "file" : "Examples/square.xml"
        },
        "cons":
        {
            "type" : "Constraints",
            "file" : "Examples/square.xml"
        }
    },

    "model" :
    {
        "type"   : "Multi",
        "models" : [ "matrix", "pbc" , "load" , "cons" ],

        "matrix" :
        {
            "type"     : "Solid",
            "thickness" : 1,
            "material" : 
            {"type" : "PlaneStrain", "young" : 900, "poisson" : 0.3 },
            "shape" : { "type" : "Tri3" }
        },

        "pbc" :
        {
            "type"     : "Periodic",
            "coarsenFactor" : 0.1,
            "shape" : { "type" : "Line2" }
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
        "solver" : { "type" : "solve" }
    },
    "nonlin":
    {   
        "niter" : 2,
        "solver" : { "type" : "solve" }
    },
    "sample":
    {
        "file": "Examples/lodi.dat",
        "dofs": [2,3]
    }
 }