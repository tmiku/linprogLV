'''
Program to update gravitational SME constraints with a single GW/EM counterpart data
Author: Tim Mikulski
Date: 27 July 2018

This code derives new constraints on the gravitational coefficients of the Lorentz-
violating Standard Model Extension (SME). Input values are  distance to the event,
fractional velocity difference between gravitational and EM signals, and sky position
(RA and declination) of the event.
'''

import sme_grav_utils as sgu
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


# ============================
#             MAIN
# ============================

opts = parse()

old = sgu.get_old_constraints(opts.data + 'coeffs.csv')
new = sgu.new_constraints(opts)
bounds = sgu.combine_constraints(old, new)
sgu.display_constraints(bounds)

if opts.writenew:
    sgu.write_constraints(bounds, opts)
    
if opts.update:
    sgu.update_constraints(bounds, opts)