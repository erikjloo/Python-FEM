{
    "input": {
        "modules": [
            "mesh",
            "load",
            "cons"
        ],
        "mesh": {
            "file": "Examples/rve.msh",
            "type": "Gmsh",
            "rank": 2,
            "doElemGroups": true
        },
        "load": {
            "type": "Loads",
            "file": "Examples/rve.xml"
        },
        "cons": {
            "type": "Constraints",
            "file": "Examples/rve.xml"
        }
    },
    "model": {
        "type": "Multi",
        "models": [
            "matrix",
            "fibers",
            "pbc",
            "cons"
        ],
        "matrix": {
            "type": "Solid",
            "elements": "gmsh1",
            "thickness": 1,
            "shape": {
                "type": "Tri3",
                "scheme": "Gauss"
            },
            "material": {
                "young": 3760,
                "poisson": 0.3
            }
        },
        "fibers": {
            "type": "Solid",
            "elements": "gmsh0",
            "thickness": 1,
            "shape": {
                "type": "Tri3",
                "scheme": "Gauss"
            },
            "material": {
                "young": 74000,
                "poisson": 0.2
            }
        },
        "pbc": {
            "type": "Periodic",
            "strainRate": [
                0.0,
                0.0,
                0.1
            ],
            "coarsenFactor": 0.7,
            "shape": {
                "type": "Line2",
                "scheme": "Gauss"
            }
        },
        "cons": {
            "type": "Constraints",
            "conTable": "cons"
        }
    },
    "nonlin": {
        "type": "NR",
        "niter": 2,
        "tiny": 1e-10,
        "tol": 0.0001,
        "solver": {
            "type": "solve"
        }
    },
    "sample": {
        "file": "Examples/lodi.dat",
        "dofs": [
            3006,
            3007
        ]
    },
    "control": {
        "nsteps": 4
    }
}