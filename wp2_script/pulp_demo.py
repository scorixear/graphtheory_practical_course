#!/usr/bin/env python3

import sys
import pulp

model = pulp.LpProblem("Maximising_Problem", pulp.LpMaximize)

V = [ pulp.LpVariable('R'+str(i), lowBound=0, upBound=50, cat='Integer') for i in range(2) ]

V.append( pulp.LpVariable('R3', lowBound=-10, upBound=10, cat='Continuous'))

# Objective function
model += 20 * V[0] + 10 * V[1] + 15 * V[2] , "Profit"

# Constraints
model += 3 * V[0] + 2 * V[1] + 5 * V[2] <= 55
model += 2 * V[0] + 1 * V[1] + 1 * V[2] <= 26
model += 1 * V[0] + 1 * V[1] + 3 * V[2] <= 30


# model += 5 * V[0] + 2 * V[1] + 4 * V[2] <= 57
listconstraints = []
for index, coef in [ (0, 5), (1, 2), (2, 4) ]:
   listconstraints.append(coef * V[index])

model += pulp.lpSum(listconstraints) <= 57

model.solve()
pulp.LpStatus[model.status]


for v in V:
   print(v.name,":", v.varValue)
   
print(pulp.value(model.objective))
