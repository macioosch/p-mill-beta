#!/usr/bin/env python3
from math import sqrt
import csv

g = 9.80665
m = 0.03246312408709453
Er = 1.098901098901099e8
Rr = 0.005
h_static = -(9/16.)**(1/3.0) * ( g*m/(Er*sqrt(Rr)) )**(2/3.0)

old_v = 0.0

def energy(h, v):
    return 0.5*m*v**2 + m*g*h

with open('output/steel-test.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter='\t')
    found = 0
    for row in spamreader:
        if len(row) == 3:
            new_v = float(row[2])
            if old_v >= 0 and new_v < 0:
                v = new_v
                h = float(row[1]) - h_static
                E = energy(h, v)
                print(found, E, sep='\t')
                found += 1
            old_v = new_v
