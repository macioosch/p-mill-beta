#!/usr/bin/env python3
from math import sqrt
from sys import argv
import csv

g = 9.80665
m = 3.2672563597333851e-02
Er = 1.0666666666666667e+11
radius = 0.01
h_static = -(9/16.)**(1/3.0) * ( g*m/(Er*sqrt(radius)) )**(2/3.0)

old_v = 0.0

def energy(h, v):
    return 0.5*m*v**2 + m*g*h

with open(argv[1], newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    found = 0
    for row in spamreader:
        if len(row) == 7 and not row[0].startswith("#"):
            new_v = float(row[4])
            if old_v >= 0 and new_v < 0:
                v = new_v
                h = float(row[2]) - h_static
                E = energy(h, v)
                t = float(row[0])
                print(found, E, t, sep='\t')
                found += 1
            old_v = new_v
