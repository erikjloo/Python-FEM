{
    "input": {
        "mesh": {
            "file": "Examples/square.msh",
            "type": "Gmsh",
            "rank": 2,
            "doElemGroups": false
        }
    },
    "model": {
        "type": "Multi",
        "models": [
            "matrix",
            "pbc"
        ],
        "matrix": {
            "type": "Solid",
            "elements": "All",
            "thickness": 1,
            "shape": {
                "type": "Tri3",
                "scheme": "Gauss"
            },
            "material": {
                "young": 75000,
                "poisson": 0.2
            }
        },
        "pbc": {
            "type": "Periodic",
            "strainRate": [
                0.05,
                0.0,
                0.1
            ],
            "coarsenFactor": 0.2,
            "shape": {
                "type": "Line2",
                "scheme": "Gauss"
            }
        }
    },
    "nonlin": {
        "type": "full",
        "niter": 20,
        "tiny": 1e-20,
        "tol": 0.001,
        "solver": {
            "type": "rtfreechol"
        }
    },
    "sample": {
        "file": "Examples/lodi.dat",
        "dofs": [
            2,
            3
        ]
    }
}