__doc__ = (
"""This file defines set of molecules with their mass peak magnitudes. It also invokes the processing.
It is recommended to put a (possibly, modified) copy of this file into the directory with 
your data files, then MassDecomposerApp.py will pick it up automatically.
""")

#edit to define peaks for each molecule. Be sure to define all contributions to every peak used.
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
        'm28Int': 0.150168492,
        'm29Int': 0.049915754,
        'm30Int': 0.125526537,
        'm31Int': 0.119208088,
        'm44Int': 0.01137321,
        'm45Int': 0.055181129,
        'm46Int': 0.621314238,
        'm47Int': 0.017481045,
    }),
    ('N2', {
        'm28Int': 0.0025,
        'm29Int': 0.095,
        'm30Int': 0.9025,
    }),    
    ('CO2nat', {
        'm44Int': 1.0,
    }),    
    ('CO2iso', {
        'm44Int': 0.359838915,
        'm45Int': 0.51355651 ,
        'm46Int': 0.004914336, #assuming low 18O due to T-induced isotopic exchange with ZnO
        'm47Int': 0.007013664, #assuming low 18O due to T-induced isotopic exchange with ZnO
    }),
]

# Set weights for peaks. Higher weight mean more importance of the corresponding equation 
# in the system. 
# Small or zero weights can be useful to exclude peaks that are noisy, drifty, or when not all 
# molecules that contribute to the peak are listed.
weights = {
    'm28Int': 0.1,
}

import MassDecomposer as MD
import sys
for arg in sys.argv[1:]:
    MD.processFile(arg, molecs, weights, decimal_separator= '.')