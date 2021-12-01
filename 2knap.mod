# This is the 2-knapsack formulation for switched shunts
# 28 October 2020
# Daniel Bienstock

option knitroampl_auxfiles rc;

param Asize >= 0 integer; #number of integer variables

param n{j in 1..Asize}; #upper bound
param b{j in 1..Asize}; #
param binit;

var x{j in 1..Asize} integer >= 0, <= n[j];
var s{j in 1..2} >= 0;

minimize difference:
  s[1] + s[2];

subject to Def:
binit = s[1] - s[2] + sum{j in 1..Asize} b[j]*x[j];

subject to Overkill:
s[1] + s[2] = 0

# impose this directly in variable definition above
#subject to upbound{j in 1..Asize}: #constraint (61) upper bound
#x[j] <= n[j];
