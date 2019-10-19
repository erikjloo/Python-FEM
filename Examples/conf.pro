{
    "input": {
        "modules": [
            "mesh",
            "load",
            "cons"
        ],
        "mesh": {
            "file": "Examples/truss.xml",
            "type": "XML",
            "rank": 2,
            "doElemGroups": false
        },
        "load": {
            "type": "Loads",
            "file": "Examples/truss.xml"
        },
        "cons": {
            "type": "Constraints",
            "file": "Examples/truss.xml"
        }
    },
    "model": {
        "type": "Multi",
        "models": [
            "truss",
            "load",
            "cons"
        ],
        "truss": {
            "type": "Truss",
            "elements": "All",
            "area": 1,
            "young": 10000000
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
    "nonlin": {
        "type": "full",
        "niter": 6,
        "tiny": 1e-10,
        "tol": 0.0001,
        "solver": {
            "type": "solve"
        }
    },
    "sample": {
        "file": "Examples/truss.dat",
        "dofs": [
            2,
            3
        ]
    },
    "control": {
        "nsteps": 2
    }
}