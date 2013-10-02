#!/usr/bin/python

#A program for a simulated annealing algorithm
import sys
import re
import random
import math

## GLOBALS ##
# input #
constants = []
num_iterations = 0

# vary these #
temp = 400.0 #starting effective temperature
cool = .7 #scaling factor for temp
start_x = 0

# output #
best = 0
best_energy = 0

# print an error message #
def error(str):
    print("ERROR: "+str)

# pick initial x point #
def pick_initial_x():
    x = start_x
    return x
    
# move from x to a new point.  x can move either left or right #
def move_x(x):
    global temp
    x1 = x + 0.5
    x2 = x - 0.5
    e = get_energy(x)
    e1 = get_energy(x1)
    e2 = get_energy(x2)
    # if x1 is better than or equal to x and x2, return x1 #
    if (e1 == min(e, e2, e1)):
        return x1
    # if x2 is better than or equal to x and x2, return x2 #
    if (e2 == min(e, e1, e2)):
        return x2
    # if x is better than x1 and x2, make a judgement call #
    prob1 = math.exp(-e1/temp)
    prob2 = math.exp(-e2/temp)
    r = random.random()
    # allow some bad choices to happen at higher temperatures #
    if (prob1 == max(r, prob1, prob2)):
        return x1
    if (prob2 == max(r, prob2, prob1)):
        return x2
    return x
    
# cooling schedule function #
def cooling():
    global temp
    temp *= cool
    
# calculate energy from current x using energy function #
def get_energy(x):
    i = 0
    y = 0
    for c in constants:
        y += c*(x**i)
        i+=1
    return y

def simulate_annealing():
    x = pick_initial_x()
    global best
    best = x
    global best_energy
    best_energy = get_energy(x)
    next_x = x
    last_x = x
    e = get_energy(x)
    for iter in range (0, num_iterations):
        e = get_energy(x)
        # remember the best position if it was encountered. #
        if best_energy > e:
            best_energy = e
            best = x
        next_x = move_x(x)
        # in case x has stopped moving #
        # (prevent math.exp from overloading) #
        if (x == next_x or next_x == last_x):
            break
        last_x = x
        x = next_x
        cooling()
    return (x,e)
    
## MAIN FUNCTION ##

# get file and extract each line #
file = open(sys.argv[1])
lines = file.readlines()
if (len(lines) != 2):
    error("file has incorrect number of lines")
    sys.exit()
for l in range(0,len(lines)):
    lines[l]=lines[l].rstrip()

# get input variables #
str_constants = re.split(r'\s', lines[0])
for s in str_constants:
    constants.append(float(s))
num_iterations = int(lines[1])

# simulate the algorithm #
last = simulate_annealing()
print("Algorithm returns:\tx = "+str(best)+"\tenergy = "+str(best_energy))
print("Algorithm's last x was:\tx = "+str(last[0])+"\tenergy = "+str(last[1]))