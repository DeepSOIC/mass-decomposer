molecs = [
    ('NO', {
        'm29Int': 0.021892103,
        'm30Int': 0.072973677,
        'm31Int': 0.905134219,
    }),
    ('CO', {
        'm28Int': 0.436291522,
        'm29Int': 0.498264750,
        'm30Int': 0.058254834,
        'm31Int': 0.007188894,
    }),
    ('N2O', {
        'm28Int': 0.000321429,
        'm29Int': 0.012214286,
        'm30Int': 0.123892857,
        'm31Int': 0.149285714,
        'm44Int': 0.001785714,
        'm45Int': 0.067857143,
        'm46Int': 0.644642857,
    }),
    ('N2', {
        'm28Int': 0.0025,
        'm29Int': 0.095,
        'm30Int': 0.9025,
    }),    
    ('CO2', {
        'm44Int': 1.0,
    }),    
]

weights = {
    'm28Int': 0.1,
}

import MassDecomposer as MD
import sys
for arg in sys.argv[1:]:
    MD.processFile(arg, molecs, weights, decimal_separator= '.')