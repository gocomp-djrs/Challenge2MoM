# Model used to solve base case ACOPF for the GO Competition 2 from Team NU-Columbia-Artelys.
# 21 Mar 2019
# Richard Waltz


param I >= 0 integer; #numbuses
param G >= 0 integer; #numgens (active only)
#param E >= 0 integer; #active non-transformer branches (not used)
#param F >= 0 integer; #active transformer branches (not used)
param LC >= 0 integer; #linecount=2*(E+F)
param Hg{1..G} >= 0; #sample pts for cost gen interpolation

#param pi := 3.14;
#param sbase; # not used

# bounds and other bus data
param nvlo{1..I};     #lobnd on voltage magnitudes
param nvhi{1..I};     #upbnd on voltage magnitudes
param nvalo{1..I};    #lobnd on voltage angles
param nvahi{1..I};    #upbnd on voltage angles
param bcslopu{1..I};  #lobnd on shunts
param bcshipu{1..I};  #upbnd on shunts
param area{1..I};  #bus area

# generator bounds
param minpg{1..G};    #lobnd on real power generation
param maxpg{1..G};    #upbnd on real power generation
param minqg{1..G};    #lobnd on reactive power generation
param maxqg{1..G};    #upbnd on reactive power generation

# Balance equation parameters
param loadP{1..I};
param fsglpu{1..I};
param loadQ{1..I};
param fsblpu{1..I};

# Info on generators and lines at each bus
param gencount{1..I} integer; # number of gens at each bus
set gensatbus {i in 1..I: gencount[i]>0}; # set of generators at each bus
param linecount{1..I} integer; # number of lines connected to each bus
set linesatbus {i in 1..I: linecount[i]>0}; # set of lines connected to each bus
#param genisactive{1..G} integer; # may be set to 0 by contingency
#param lineisactive{1..LC} integer; # may be set to 0 by contingency

# provided initval parameters

# line flow coefficients / values
param Gff{1..LC}; # ge
param Gft{1..LC}; # -ge
param Bff{1..LC}; # (be + be^CH / 2)
param Bft{1..LC}; # -be
param transformer{1..LC};
param linebusi{1..LC};
param linebusj{1..LC};
param shiftangle{1..LC};
param rating{1..LC};

# Primary Variables
var busVmags{i in 1..I}, >=nvlo[i], <=nvhi[i];
var busVangles{i in 1..I}, >=nvalo[i], <=nvahi[i];
var busCshunts{i in 1..I}, >=bcslopu[i], <=bcshipu[i];
var genP{i in 1..G}, >=minpg[i], <=maxpg[i];
var genQ{i in 1..G}, >=minqg[i], <=maxqg[i];

# Variables to construct linear gencost objective from interpolation
#param cgh{g in 1..G, h in 1..Hg[g]}, >= 0.0; # gen cost coefficients in (2)
#param pgh{g in 1..G, h in 1..Hg[g]}, >= 0.0; # power values in (3)
param cgh{g in 1..G, h in 1..Hg[g]}; # gen cost coefficients in (2); can be negative?
param pgh{g in 1..G, h in 1..Hg[g]}; # power values in (3); can be negative?
var tgh{g in 1..G, h in 1..Hg[g]}, >=0, <=1; # interpolation coefficients

# Defined variables for gencost (2)
var gencost{g in 1..G} = sum{h in 1..Hg[g]} cgh[g,h]*tgh[g,h];

# Defined variables for line flow definitions
var lineP{i in 1..LC} = Gff[i]*busVmags[linebusi[i]]^2 +
	  	      ( Gft[i]*cos(busVangles[linebusi[i]] - busVangles[linebusj[i]] - shiftangle[i]) +
		        Bft[i]*sin(busVangles[linebusi[i]] - busVangles[linebusj[i]] - shiftangle[i]) ) * (busVmags[linebusi[i]]*busVmags[linebusj[i]]) ;
var lineQ{i in 1..LC} = -Bff[i]*busVmags[linebusi[i]]^2 +
	  	      (-Bft[i]*cos(busVangles[linebusi[i]] - busVangles[linebusj[i]] - shiftangle[i]) +
		        Gft[i]*sin(busVangles[linebusi[i]] - busVangles[linebusj[i]] - shiftangle[i]) ) * (busVmags[linebusi[i]]*busVmags[linebusj[i]]) ;

# Objective function (generation cost)
# minimize bus_cost - load_benefit + line_cost + transformer_cost + gen_cost
minimize f:
	#sum{g in 1..G} genP[g];
        sum{g in 1..G} gencost[g];

# Power interpolation constraints (3)
subject to PowerInterp{g in 1..G}:
        genP[g] = sum{h in 1..Hg[g]} pgh[g,h]*tgh[g,h];

# Interpolation coefficient normalization (5)
subject to InterpNorm{g in 1..G: Hg[g]>0}:
        sum{h in 1..Hg[g]} tgh[g,h] = 1;

# Real power balance constraints
subject to PBalance{i in 1..I}:
        sum{g in gensatbus[i]: gencount[i]>0} genP[g] - loadP[i] - fsglpu[i]*busVmags[i]^2 -
        sum{l in linesatbus[i]: linecount[i]>0} lineP[l] = 0;

# Reactive power balance constraints
subject to QBalance{i in 1..I}:
        sum{g in gensatbus[i]: gencount[i]>0} genQ[g] - loadQ[i] - (-fsblpu[i] - busCshunts[i])*(busVmags[i]^2) -
        sum{l in linesatbus[i]: linecount[i]>0} lineQ[l] = 0;

# Line current constraints (nontransformers); squared form
subject to LineCurrentNT{i in 1..LC: transformer[i]==0}:
        lineP[i]^2 + lineQ[i]^2 - (rating[i]*busVmags[linebusi[i]])^2 <= 0;


# Line current constraints (transformers); squared form
subject to LineCurrentT{i in 1..LC: transformer[i]>0}:
        lineP[i]^2 + lineQ[i]^2 - rating[i]^2 <= 0;

#data;
