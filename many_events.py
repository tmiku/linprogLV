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
import cvxopt
from scipy.optimize import linprog
from scipy.special import sph_harm

bounds_order = ['s_00','s_10','Re s_11','Im s_11','s_20','Re s_21','Im s_21','Re s_22','Im s_22']
harmonics = [(0,0), (1,0), (1,1), (2,0), (2,1), (2,2)]


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
    parser.add_option('-p', '--package', default = 'scipy',
    help = 'Package to use for solving linear program. Default = scipy, or choose cvxopt')
    parser.add_option('-s', '--solver', default= 'simplex',
    help = 'Method of solving linear program. Default = simplex for scipy/built-in for cvxopt')

    opts = parser.parse_args()[0]

    return opts


def optimize(coeff, events, package='scipy', solver='simplex'):

    c_upper = np.zeros(9)
    c_lower = np.zeros(9)
    c_upper[coeff] = -1
    c_lower[coeff] = 1

    a = np.zeros((2 * len(events), 9))
    b = np.array([])
    n_bound = 0

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
                temp = np.append(temp, .5 * (-1 ** ii[0]) * np.real(sph_harm(ii[1], ii[0], ph, th)))
            else:
                temp = np.append(temp, (-1 ** ii[0]) * np.real(sph_harm(ii[1], ii[0], ph, th)))
                temp = np.append(temp, -1 * (-1 ** ii[0]) * np.imag(sph_harm(ii[1], ii[0], ph, th)))
        
        a[n_bound] = temp
        n_bound += 1
        b = np.append(b, dvmax)
        a[n_bound] = -temp
        n_bound += 1
        b = np.append(b, -dvmin)

    if package == 'scipy':
        min_sol = linprog(c_lower, a, b, bounds=(None,None), method=solver)
        max_sol = linprog(c_upper, a, b, bounds=(None,None), method=solver)
        lower = min_sol.x[coeff]
        upper = max_sol.x[coeff]
    
    elif package == 'cvxopt':
        (a, b_max, b_min, c_upper, c_lower) = map(cvxopt.matrix, (a, b_max, b_min, c_upper, c_lower))
        upper = cvxopt.solvers.lp(c_upper, a, b, solver='glpk' if solver=='glpk' else None)['x'][coeff]
        lower = cvxopt.solvers.lp(c_lower, a, b, solver='glpk' if solver=='glpk' else None)['x'][coeff]
    
    bound = [lower, upper]
    

    return np.array(bound)


# ============================
#            MAIN
# ============================

if __name__ == '__main__':
    opts = parse()
    old = sgu.get_old_constraints(opts.data + 'coeffs.csv')
    events = sgu.read_events(opts.data)

    bounds = np.array([[0,0]])
    for i in range(len(bounds_order)):
        bounds = np.append(bounds, np.array([optimize(i, events)]),axis=0)
    new = np.delete(bounds,0,0)
    new = np.transpose(new)
    new = new.tolist()

    bounds = sgu.combine_constraints(old, new)
    sgu.display_constraints(bounds)

    if opts.writenew:
        sgu.write_constraints(bounds, opts, 'many')
        
    if opts.update:
        sgu.update_constraints(bounds, opts, 'many')