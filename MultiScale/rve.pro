 {
    "control": { "nsteps" : 4 },
    "input":
    {
        "modules" : ["mesh" , "load" , "cons" ],
        "mesh" :
        {
            "type" : "Gmsh",
            "file" : "_mshes/rve26.msh",
            "rank" : 2,
            "doElemGroups" : true
        },
        "load" : 
        { 
            "type" : "Loads",
            "file" : "Examples/rve.xml"
        },
        "cons":
        {
            "type" : "Constraints",
            "file" : "Examples/rve.xml"
        }
    },

    "model" :
    {
        "type"   : "Multi",
        "models" : [ "matrix", "fibers", "pbc" , "cons" ],

        "matrix" :
        {
            "type"     : "Solid",
            "elements" : "gmsh1",
            "thickness" : 1,

            "material" : 
            {
                "type" : "PlaneStrain",
                "young" : 3760,
                "poisson" : 0.3
            },
            
            "shape" : { "type" : "Tri3" }
        },

        "fibers" : 
        {
            "type" : "Solid",
            "elements" : "gmsh0",
            "thickness" : 1,

            "material" : 
            {
                "type" : "PlaneStrain",
                "young" : 74000,
                "poisson" : 0.2
            },

            "shape" : { "type" : "Tri3" }
        },
        
        "pbc" :
        {
            "type"     : "Periodic",
            "coarsenFactor" : 0.7,
            "strain" : [0.1, -0.0428571428571429, 0.0],
            "strainRate" : [0.0, 0.0, 0.1],
            "shape": { "type" : "Line2" }
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
    "nonlin":
    {
        "type" : "NR",
        "niter" : 2,
        "tol" : 1e-4,
        "solver" : { "type" : "solve" }
    },
    "sample":
    {
        "file": "Examples/lodi.dat",
        "dofs": [3006,3007]
    }
 }