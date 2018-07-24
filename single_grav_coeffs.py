'''
Program to update gravitational SME constraints with GW/EM counterpart data
Author: Tim Mikulski
Date: In progress, last updated 24 July 2018

This code derives new constraints on the gravitational coefficients of the Lorentz-
violating Standard Model Extension (SME). Input values are  distance to the event, time
difference between gravitational and EM signals, the RA and declination of the event, and
associated error with these values. This raw data form past events will be read in from a
file, and used to fully optimize the new coefficients.

Most of the new_constraints function is directly from Jay Tasson's jupyter notebook 
calc_sme_constraints90v3.ipynb used in LIGO-P1700308.
'''

import numpy as np
import texttable as tt
import optparse
from scipy.special import sph_harm

def parse():
    # Parse instructions from command line
    parser = optparse.OptionParser()

    parser.add_option('-d', '--data', help = 'Data file for prior constraints')
    parser.add_option('-w', '--writenew', help = 'Write these results to coeffs.csv and events.csv (default false)', default = store_true)

    parser.add_option('--dvmin', help = 'Minimum fractional difference in velocity between GW and EM signal')
    parser.add_option('--dvmax', help = 'Maximum fractional difference in velocity between GW and EM signal')
    parser.add_option('--ra', help = 'Right ascension of event, use hours')
    parser.add_option('--dec', help = 'Declination of event, use degrees')

    opts = parser.parse_args()

    return opts

# ============================
#         DEFINITIONS
# ============================

def get_old_constraints(loc):
    data_file = open(loc)
    lines = data_file.readlines()
    mins = lines[1].split(',')
    maxes = lines[2].split(',')
    cons = [mins,maxes]
    return cons

def new_constraints(dvmin, dvmax, ra, dec):
    ph = ra * 360/24 * np.pi / 180 # azimuthal angle in radians
    th = (90 - dec) * np.pi / 180 # polar angle in radians

    mins = []
    maxes = []
    d = 4


    for j in range(d-1):
        # Need to talk to jay about this for loop -- why -1^j and not -1^(1+j)?
        m = 0
        maxes.append(dvmin/(.5*(-1)**j*np.real(sph_harm(m,j,ph,th))))
        mins.append(dvmax/(.5*(-1)**j*np.real(sph_harm(m,j,ph,th))))

        for m in range(1,j+1):
            maxes.append(dvmin/(.5*(-1)**j*np.real(sph_harm(m,j,ph,th))))
            mins.append(dvmax/(.5*(-1)**j*np.real(sph_harm(m,j,ph,th))))
            maxes.append(dvmin/(.5*(-1)**j*np.imag(sph_harm(m,j,ph,th))))
            mins.append(dvmax/(.5*(-1)**j*np.imag(sph_harm(m,j,ph,th))))

    return [mins,maxes]


# ============================
#             MAIN
# ============================

opts = parse()
old = get_old_constraints(opts.data)
new = new_constraints(opts.dvmin, opts.dvmax, opts.ra, opts.dec)
bounds = [[],[]]

for i in range(len(new[0])):

    if new[0][i] > old[0][i]:
        bounds[0].append(new[0][i])
    else:
        bounds[0].append(old[0][i])
    
    if new[1][i] < old[1][i]:
        bounds[1].append(new[1][i])
    else:
        bounds[1].append(old[1][i])

bounds_order = ['s_00','s_10','-Re s_11','Im s11','-s_20','-Re s_21','Im s_21','Re s_22','-Im s_22']
tab = tt.Texttable()
tab.header(['Lower','s','upper'])

for i in range(len(new)):
    tab.add_row([bounds[0][i], bounds_order[i], bounds[1][i]])

a = tab.draw()
print a

if opts.writenew:
    