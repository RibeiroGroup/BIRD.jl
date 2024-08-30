import psi4
import numpy as np

LiH = psi4.geometry("""
0 1
Li
H             1    R""")

psi4.set_options({
    "basis" : "aug-cc-pvtz"
})

Rvals = np.arange(0.5, 10.0, 0.1)

for R in Rvals:
    LiH.R = R
    properties('ccsd', properties=['dipole'])
