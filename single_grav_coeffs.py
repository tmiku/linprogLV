'''
Program to update gravitational SME constraints with GW/EM counterpart data
Author: Tim Mikulski
Date: In progress, last updated 2 July 2018

This code derives new constraints on the gravitational coefficients of the Lorentz-
violating Standard Model Extension (SME). Input values are  distance to the event, time
difference between gravitational and EM signals, the RA and declination of the event, and
associated error with these values. This raw data form past events will be read in from a
file, and used to fully optimize the new coefficients.

Most of the new_constraints function is directly from Jay Tasson's jupyter notebook 
calc_sme_constraints90v3 used in LIGO-P1700308.
'''

import scipy.special
import numpy as np
import texttable as tt
import optparse

def parse():
    # Parse instructions from command line
    parser = optparse.OptionParser()

    parser.add_option('-d', '--data', help = 'Data file for prior constraints')
    parser.add_option('-w', '--writenew', help = 'Update data file with new constraints (default false)', default = store_true)

    parser.add_option('--dvmin', help = 'Minimum fractional difference in velocity between GW and EM signal')
    parser.add_option('--dvmax', help = 'Maximum fractional difference in velocity between GW and EM signal')
    parser.add_option('--ra')
    parser.add_option('--dec')

    opts = parser.parse_args()

    return opts

def old_constraints(loc):
    data_file = open(loc)
    lines = data_file.readlines()
    mins = lines[1].split(',')
    maxes = lines[2].split(',')
    cons = [mins,maxes]
    return cons

def new_constraints(dvmin, dvmax, ra, dec):
    th = ra * np.pi / 180
    ph = (90 - dec) * np.pi / 180

    mins = []
    maxes = []
    d = 4

    for j in range(d-1):
        m = 0
        mins.append(dvmin/(scipy.special.sph_harm(m, j, th, ph)))

# ============================
#             MAIN
# ============================

opts = parse()
old = old_constraints(opts.data)
new = new_constraints(opts.dvmin, opts.dvmax, opts.ra, opts.dec)
comb = [[],[]]

for i in range(len(new[0])):

    if new[0][i] > old[0][i]:
        comb[0].append(new[0][i])
    else:
        comb[0].append(old[0][i])
    
    if new[1][i] < old[1][i]:
        comb[1].append(new[1][i])
    else:
        comb[1].append(old[1][i])

upper