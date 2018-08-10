'''
General-use functions for the sme_grav_constraints codebase
Author: Tim Mikulski
Date: 27 July 2018

This code contains functions used by multiple other scripts in the SME_grav_constraints
codebase that refines constratints on the coefficients for Lorentz violation in the
gravitational sector of the Standard Model Extension (SME).

The new_constraints function is derived from code in Jay Tasson's jupyter notebook 
calc_sme_constraints90v3.ipynb used for these coefficients in LIGO-P1700308.
'''

from __future__ import print_function
import numpy as np
import os
import texttable as tt
from scipy.special import sph_harm
from shutil import copyfile

bounds_order = ['s_00','s_10','Re s_11','Im s_11','s_20','Re s_21','Im s_21','Re s_22','Im s_22']


def get_old_constraints(loc):
    with open(loc, 'r') as old_file:
        lines = old_file.readlines()
    mins = map(float,lines[1].split(','))
    maxes = map(float,lines[2].split(','))
    cons = np.array([mins,maxes])
    return cons


def combine_constraints(old, new):
    bounds = [[],[]]
    flags=[[],[]]


    for i in range(len(new[0])):
    
        if new[0][i] > old[0][i]:
            bounds[0].append(new[0][i])
            flags[0].append(i)
        else:
            bounds[0].append(old[0][i])
        
        if new[1][i] < old[1][i]:
            bounds[1].append(new[1][i])
            flags[1].append(i)
        else:
            bounds[1].append(old[1][i])
    
    if len(flags[0]) > 0 or len(flags[1]) > 0:
        print("Constraints improved:")
        for i in flags[0]:
            print(bounds_order[i] + ' lower bound')
        for i in flags[1]:
            print(bounds_order[i] + ' upper bound')
    else:
        print("No constraints improved.")

    return bounds


def display_constraints(bounds):
    tab = tt.Texttable()
    tab.header(['Lower','s','upper'])

    for i in range(len(bounds[0])):
        tab.add_row([bounds[0][i], bounds_order[i], bounds[1][i]])

    a = tab.draw()
    print(a)


def update_constraints(bounds, opts, loc = 'data/'):
    with open(loc + 'events.csv','a') as events:
        events.write('%s,%s,%s,%s\n' % (opts.ra, opts.dec, opts.dvmax, opts.dvmin))

    header = ','.join(bounds_order) + '\n'
    mins_str = ','.join([str(b) for b in bounds[0]]) + '\n'
    maxes_str = ','.join([str(b) for b in bounds[1]]) + '\n'

    with open(loc + 'coeffs.csv','w') as coeffs:
        coeffs.write(header + mins_str + maxes_str)


def write_constraints(bounds, opts, loc = 'data/'):
    if 'outputs' not in os.listdir('.'):
        os.mkdir('outputs')
        
    copyfile(loc + 'events.csv', 'outputs/events.csv')
    with open('outputs/events.csv','a') as events:
        events.write('%s,%s,%s,%s\n' % (opts.ra, opts.dec, opts.dvmax, opts.dvmin))
    
    header = ','.join(bounds_order) + '\n'
    mins_str = ','.join([str(b) for b in bounds[0]]) + '\n'
    maxes_str = ','.join([str(b) for b in bounds[1]]) + '\n'

    with open('outputs/coeffs.csv','w') as coeffs:
        coeffs.write(header + mins_str + maxes_str)


def read_events(data):
    events_file = open(data + 'events.csv')
    events = []
    first_line = True
    for line in events_file:
        if first_line:
            first_line = False
        else:
            events.append(map(float,line.rstrip().split(',')))
    return np.array(events)