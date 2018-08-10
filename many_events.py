'''
Program to update gravitational SME constraints with nine or more GW/EM counterparts
Author: Tim Mikulski
Date: 27 July 2018

This code derives new constraints on the gravitational coefficients of the Lorentz-
violating Standard Model Extension (SME) using linear programming. It requires that at
least nine events be present in data/events.csv.
'''

import sme_grav_utils as sgu
import numpy as np
import optparse
from scipy.optimize import linprog

bounds_order = ['s_00','s_10','Re s_11','Im s_11','s_20','Re s_21','Im s_21','Re s_22','Im s_22']
harmonics = [(0,0), (1,0), (1,1), (2,0), (2,1), (2,2)]
imaginary_idxs = [3, 6, 8]


# ============================
#         DEFINITIONS
# ============================


def parse():
    # Parse instructions from command line
    parser = optparse.OptionParser()

    parser.add_option('-d', '--data', default = 'data/',
    help = 'Directory containing prior events and constraints (events.csv and coeffs.csv)')
    parser.add_option('-w', '--writenew', action = 'store_true', default = False,
    help = 'Write these results to new copies of coeffs.csv and events.csv (default false)')
    parser.add_option('-u', '--update', action = 'store_true', default = False,
    help = 'Update coeffs.csv and events.csv with new data (default false)')
    
    opts = parser.parse_args()[0]

    return opts


def optimize(coeff, events):

    cmax = cmin = np.zeros(9)
    cmax[coeff] = -1
    cmin[coeff] = 1

    a = np.zeros(9)
    bmax = bmin = np.array([])
    for event in events:

        ra = event[0]
        dec = event[1]
        dvmax = event[2]
        dvmin = event[3]
        ph = ra * 360/24 * np.pi / 180
        th = (90 - dec) * np.pi / 180
        temp = np.array([])

        for ii in harmonics:
            if ii[1] == 0:
                temp = np.append(temp, .5 * (-1 ** ii[0]) * sph_harm(ii[1], ii[0], ph, th))
            else:
                temp = np.append(temp, .5 * (-1 ** ii[0]) * np.real(sph_harm(ii[1], ii[0], ph, th)))
                temp = np.append(temp, .5 * (-1 ** ii[0]) * np.imag(sph_harm(ii[1], ii[0], ph, th)))
        
        a = np.append(a, temp, axis=0)
        bmax = np.append(bmax,dvmax)
        bmin = np.append(bmin, dvmin)
    
    a = np.delete(a, 0, 0)
    dvmax_upper = linprog(cmax, a, bmax).x
    dvmax_lower = linprog(cmin, a, bmax).x
    dvmin_upper = linprog(cmax, a, bmin).x
    dvmin_lower = linprog(cmin, a, bmin).x

    bounds = sgu.combine_constraints([dvmax_lower, dvmax_upper], [dvmin_lower, dvmin_upper])
    return bounds


# ============================
#            MAIN
# ============================

opts = parse()
old = sgu.get_old_constraints(opts.data + 'coeffs.csv')
events = sgu.read_events(opts.data)

bounds = []
for i in len(bounds_order):
    if len(bounds) == 0:
        bounds.append(optimize(i, events))
    else:
        temp_bounds = optimize(i, event)
        bounds = sgu.combine_constraints(bounds, temp_bounds)

sgu.display_constraints(bounds)

if opts.writenew:
    sgu.write_constraints(bounds, opts)
    
if opts.update:
    sgu.update_constraints(bounds, opts)