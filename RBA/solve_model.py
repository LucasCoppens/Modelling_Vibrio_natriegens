"""Solve standard RBA problem."""

# python 2/3 compatibility
from __future__ import division, print_function

import os
import sys
sys.path.append("/Users/lucascoppens/Documents/Phd/Active/Vnat modelling/Vnat_v5/RBA/RBApy")

# package imports
import rba
import re
import copy
import numpy as np
import pandas as pd



def main():
    model = rba.RbaModel.from_xml("model")

    res = model.solve(bissection_tol = 0.01)

    res.write_fluxes('fluxes.csv', file_type="csv")
    print("Growth rate: ", res.mu_opt)

if __name__ == '__main__':
    main()
