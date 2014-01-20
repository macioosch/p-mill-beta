#!/usr/bin/env python3
from math import sqrt
from sys import argv
import csv

g = 9.80665
m = 4.084070449666731e-3
Er = 1e+11
radius = 0.005

old_v = 0.0

def energy(h, v):
    return m*(0.5*v**2 + g*h)

with open(argv[1], newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    found = 0
    for row in spamreader:
        if not row[0].startswith("#"):
            new_v = float(row[4])
            if old_v >= 0 and new_v < 0:
                v = new_v
                h = float(row[2])
                E = energy(h, v)
                t = float(row[0])
                print(found, E, t, sep='\t')
                found += 1
            old_v = new_v
