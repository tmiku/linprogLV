'''
Program to update gravitational SME constraints with a single GW/EM counterpart
Author: Tim Mikulski
Date: 27 July 2018

This code derives new constraints on the gravitational coefficients of the Lorentz-
violating Standard Model Extension (SME). Input values are  distance to the event,
fractional velocity difference between gravitational and EM signals, and sky position
(RA and declination) of the event.
'''

from __future__ import print_function
from scipy.special import sph_harm
import sme_grav_utils as sgu
import numpy as np
import optparse


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

    parser.add_option('--dvmin', help = 'Minimum fractional GW/EM velocity difference', action = 'store', type = 'float')
    parser.add_option('--dvmax', help = 'Maximum fractional GW/EM velocity difference', action = 'store', type = 'float')
    parser.add_option('--ra', help = 'Right ascension of event, use hours', action = 'store', type = 'float')
    parser.add_option('--dec', help = 'Declination of event, use degrees', action = 'store', type = 'float')

    (opts, args) = parser.parse_args()

    if len(args) != 0:
        print("Warning: some options not read.")

    return opts


def new_constraints(opts):
    ph = opts.ra * 360/24 * np.pi / 180 # azimuthal angle in radians
    th = (90 - opts.dec) * np.pi / 180 # polar angle in radians

    mins = []
    maxes = []
    negs = []
    neg_indices = {(0,0):0, (1,0):1, (1,1,'re'):2, (1,1,'im'):3, (2,0):4, (2,1,'re'):5, (2,1,'im'):6, (2,2,'re'):7, (2,2,'im'):8}
    d = 4

    for j in range(d-1):
        m = 0
        maxes.append(opts.dvmax/(.5*(-1)**j*np.real(sph_harm(m,j,ph,th))))
        mins.append(opts.dvmin/(.5*(-1)**j*np.real(sph_harm(m,j,ph,th))))
        if (sph_harm(m,j,ph,th) > 0 and j == 1) or (sph_harm(m,j,ph,th) < 0 and j != 1):
            negs.append((j,m))

        for m in range(1,j+1):
            maxes.append(opts.dvmax/((-1)**j*np.real(sph_harm(m,j,ph,th))))
            mins.append(opts.dvmin/((-1)**j*np.real(sph_harm(m,j,ph,th))))
            maxes.append(-opts.dvmax/((-1)**j*np.imag(sph_harm(m,j,ph,th))))
            mins.append(-opts.dvmin/((-1)**j*np.imag(sph_harm(m,j,ph,th))))

            if (np.real(sph_harm(m,j,ph,th)) > 0 and j == 1) or (np.real(sph_harm(m,j,ph,th)) < 0 and j != 1):
                negs.append((j,m,'re'))

            if (np.imag(sph_harm(m,j,ph,th)) < 0 and j == 1) or (np.imag(sph_harm(m,j,ph,th)) > 0 and j != 1):
                negs.append((j,m,'im'))

    mins = -abs(np.array(mins))
    maxes = abs(np.array(maxes))

    for item in negs:
        (mins[neg_indices[item]], maxes[neg_indices[item]]) = (-maxes[neg_indices[item]], -mins[neg_indices[item]])

    return np.array([mins,maxes]).round(decimals=16)

# ============================
#             MAIN
# ============================

opts = parse()

old = sgu.get_old_constraints(opts.data + 'coeffs.csv')
new = new_constraints(opts)
bounds = sgu.combine_constraints(old, new)
sgu.display_constraints(bounds)

if opts.writenew:
    sgu.write_constraints(bounds, opts, 'single')
    
if opts.update:
    sgu.update_constraints(bounds, opts, 'many')