{
    "input": {
        "modules": [
            "mesh",
            "load",
            "cons"
        ],
        "mesh": {
            "file": "Examples/2D_semicircle.msh",
            "type": "Gmsh",
            "rank": 2,
            "doElemGroups": false
        },
        "load": {
            "type": "Loads",
            "file": "Examples/2D_semicircle.xml"
        },
        "cons": {
            "type": "Constraints",
            "file": "Examples/2D_semicircle.xml"
        }
    },
    "model": {
        "type": "Multi",
        "models": [
            "model",
            "load",
            "cons"
        ],
        "model": {
            "type": "Solid",
            "elements": "All",
            "thickness": 0.1,
            "shape": {
                "type": "Tri3",
                "scheme": "Gauss1"
            },
            "material": {
                "young": 100000,
                "poisson": 0.2
            }
        },
        "load": {
            "type": "PointLoad",
            "loadTable": "load"
        },
        "cons": {
            "type": "Constraints",
            "conTable": "cons"
        }
    },
    "linsolve": {
        "solver": {
            "type": "lstsq"
        }
    }
}