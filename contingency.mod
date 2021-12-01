# This is a constrained NLP model formulation for the GO Competition from Team djrs.
# This one uses an L1 relaxation of balance and line current constraints.
# 21 Mar 2019
# Richard Waltz


param I >= 0 integer; #numbuses
param G >= 0 integer; #numgens (active only)
#param E >= 0 integer; #active non-transformer branches (not used)
#param F >= 0 integer; #active transformer branches (not used)
param LC >= 0 integer; #linecount=2*(E+F)
param Hg{1..G} >= 0; #sample pts for cost gen interpolation
#param J >= 0 integer; #numloads (need this? load defined for each bus?)

#param pi := 3.14;
#param sbase := 100; # not used
#param maxsoftviol:= 0.019; #hard limit on soft constraint violation
#param maxsoftviol:= 0.49; #hard limit on soft constraint violation
param maxviolBalance >=0; #hard limit on soft constraint violation for balance equations
param maxviolLC >=0;      #hard limit on soft constraint violation for line current constraints

# bounds and other bus data
param nvlo{1..I};     #lobnd on voltage magnitudes
param nvhi{1..I};     #upbnd on voltage magnitudes
param nvalo{1..I};    #lobnd on voltage angles
param nvahi{1..I};    #upbnd on voltage angles
param bcslopu{1..I};  #lobnd on shunts
param bcshipu{1..I};  #upbnd on shunts
param area{1..I};     #bus area
param freeVmags{1..I}; # = 1 if bus can change voltage, = 0 otherwise

# Base solution as data
param baseVmags{1..I};    # used in model and also may be used to define initial pt
param baseVangles{1..I};  # may be used to define initial pt
param baseCshunts{1..I};  # may be used to define initial pt
param basegenP{1..G};     # used in model and also may be used to define initial pt
param basegenQ{1..G};     # may be used to define initial pt

#lower and upper bounds on Delta_k
param minDelta_Scen;
param maxDelta_Scen;

# generator bounds
param minpg{1..G};    #lobnd on real power generation
param maxpg{1..G};    #upbnd on real power generation
param minqg{1..G};    #lobnd on reactive power generation
param maxqg{1..G};    #upbnd on reactive power generation

# some upper bounds on comp variables
#param maxVmagdrop{1..I};
#param maxVmagrise{1..I};
param maxshortp{1..G};

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
param alphapart{1..G};  # = alpha when participating, = 0 if not
param yesnopart{1..G};  # = 1 when participating, = 0 if not
param movinggen{1..G};  # = 1 when allowed to vary with affine rule for participation
param notoperating{1..G}; # = 0 when out in the contingency
param lowDelta;

# Primary Variables
var busVmags{i in 1..I}, >=nvlo[i], <=nvhi[i];
#var busVmagdrop{i in 1..I}, >=0, <=maxVmagdrop[i];
#var busVmagrise{i in 1..I}, >=0, <=maxVmagrise[i];
#var busVmagdrop{i in 1..I}, >=0, <=nvhi[i];
#var busVmagrise{i in 1..I}, >=0, <=nvhi[i];
var busVmagdrop{i in 1..I}, >=0, <=nvhi[i]-nvlo[i];
var busVmagrise{i in 1..I}, >=0, <=nvhi[i]-nvlo[i];
#var busVmagdrop{i in 1..I}, >=0, <=baseVmags[i]-nvlo[i];
#var busVmagrise{i in 1..I}, >=0, <=baseVmags[i]-nvlo[i];
var busVangles{i in 1..I}, >=nvalo[i], <=nvahi[i];
var busCshunts{i in 1..I}, >=bcslopu[i], <=bcshipu[i];
var genP{i in 1..G}, >=minpg[i], <=maxpg[i];
var genPdrop{i in 1..G}, >=0, <=maxpg[i] - minpg[i];
var genPrise{i in 1..G}, >=0, <=maxpg[i] - minpg[i];
var shortageP{i in 1..G}, >=0, <=maxpg[i] - minpg[i];
#var shortageP{i in 1..G}, >=0, <=maxshortp[i];
var excessP{i in 1..G}, >=0, <=maxpg[i] - minpg[i];
var genQ{i in 1..G}, >=minqg[i], <=maxqg[i];
var shortageQ{i in 1..G}, >=0, <=maxqg[i] - minqg[i];
var busSumshortageQ{i in 1..I}, >=0;
var excessQ{i in 1..G}, >=0, <=maxqg[i] - minqg[i];
var busSumexcessQ{i in 1..I}, >=0;
#var Delta_scen, >=minDelta_Scen, <=maxDelta_Scen;
var Delta_scen, >= lowDelta;

# Variables to construct linear gencost objective from interpolation
#param cgh{g in 1..G, h in 1..Hg[g]}, >= 0.0; # gen cost coefficients in (2)
#param pgh{g in 1..G, h in 1..Hg[g]}, >= 0.0; # power values in (3)
param cgh{g in 1..G, h in 1..Hg[g]}; # gen cost coefficients in (2); can be negative?
param pgh{g in 1..G, h in 1..Hg[g]}; # power values in (3); can be negative?
var tgh{g in 1..G, h in 1..Hg[g]}, >=0, <=1; # interpolation coefficients

# Variables to construct linear loadben objective from interpolation
param clh{g in 1..G, h in 1..Hg[g]}; # load benefit coefficients in (15/16); can be negative?
param plh{g in 1..G, h in 1..Hg[g]}; # load values in (13); can be negative?
var tlh{g in 1..G, h in 1..Hg[g]}, >=0, <=1; # interpolation coefficients

# Relaxation variables
var relaxPplus{i in 1..I}, >=0;
var relaxPminus{i in 1..I}, >=0;
var relaxQplus{i in 1..I}, >=0;
var relaxQminus{i in 1..I}, >=0;
#var relaxLineNT{i in 1..LC: transformer[i]==0}, >=0;
#var relaxLineT{i in 1..LC: transformer[i]>0}, >=0;
var relaxLineNT{i in 1..LC: transformer[i]==0}, >=0, <=maxviolLC;
var relaxLineT{i in 1..LC: transformer[i]>0}, >=0 <=maxviolLC;

# Defined variables for gencost (2)
var gencost{g in 1..G} = sum{h in 1..Hg[g]} cgh[g,h]*tgh[g,h];

#TBD defined variables for loadP/loadQ based on the fraction that has cleared, deterined by variable loadFrac

# Defined variables for line flow definitions
var lineP{i in 1..LC} = Gff[i]*busVmags[linebusi[i]]^2 +
	  	      ( Gft[i]*cos(busVangles[linebusi[i]] - busVangles[linebusj[i]] - shiftangle[i]) +
		        Bft[i]*sin(busVangles[linebusi[i]] - busVangles[linebusj[i]] - shiftangle[i]) ) * (busVmags[linebusi[i]]*busVmags[linebusj[i]]) ;
var lineQ{i in 1..LC} = -Bff[i]*busVmags[linebusi[i]]^2 +
	  	      (-Bft[i]*cos(busVangles[linebusi[i]] - busVangles[linebusj[i]] - shiftangle[i]) +
		        Gft[i]*sin(busVangles[linebusi[i]] - busVangles[linebusj[i]] - shiftangle[i]) ) * (busVmags[linebusi[i]]*busVmags[linebusj[i]]) ;



# Objective function (L1 violation of balance and line current constraints)
minimize f:
	sum{i in 1..I} (relaxPplus[i] + relaxPminus[i] + relaxQplus[i] + relaxQminus[i]) + # bus cost
        sum{j in 1..LC: transformer[j]==0} (relaxLineNT[j]) +                              # line cost
        sum{k in 1..LC: transformer[k]>0} (relaxLineT[k]) +                                # transformer cost
        sum{g in 1..G} (gencost[g]);                                                       # gen cost (see base case)
        # - sum{i in 1..I} (ben[i]*]loadP[i]);                                             # TBD: load benefit


# Real power balance constraints
subject to PBalance{i in 1..I}:
        sum{g in gensatbus[i]: gencount[i]>0} genP[g] - loadP[i] - fsglpu[i]*busVmags[i]^2 -
        sum{l in linesatbus[i]: linecount[i]>0} lineP[l] = relaxPplus[i] - relaxPminus[i];

#enforce hard upper bound on violation
subject to PBalanceViolLimit{i in 1..I}: relaxPplus[i] + relaxPminus[i] <= maxviolBalance;

# Reactive power balance constraints
subject to QBalance{i in 1..I}:
        sum{g in gensatbus[i]: gencount[i]>0} genQ[g] - loadQ[i] - (-fsblpu[i] - busCshunts[i])*(busVmags[i]^2) -
        sum{l in linesatbus[i]: linecount[i]>0} lineQ[l] = relaxQplus[i] - relaxQminus[i];

#enforce hard upper bound on violation
subject to QBalanceViolLimit{i in 1..I}: relaxQplus[i] + relaxQminus[i] <= maxviolBalance;

# Line current constraints (nontransformers); squared form
subject to LineCurrentNT{i in 1..LC: transformer[i]==0}:
        lineP[i]^2 + lineQ[i]^2 - (rating[i]*busVmags[linebusi[i]])^2 <= relaxLineNT[i];


# Line current constraints (transformers); squared form
subject to LineCurrentT{i in 1..LC: transformer[i]>0}:
       lineP[i]^2 + lineQ[i]^2 - rating[i]^2 <= relaxLineT[i];

# Settling real power disjunctions. This is done in constraint Pset, explained as follows
# Suppose first that yesnopart[i] = 1.
#         if movinggen[i] = 1 generator moves according to affine rule, and if movinggen[i] = 0 then generator is at upper bound
# Suppose next that yesnopart[i] = 0.
#         if notoperating[i] = 0, i.e. generator is operating, then generator stays fixed at base value.  If notoperating[i] = 1, then generator is at zero


subject to Pset{i in 1..G}:
        genP[i] = yesnopart[i]*( movinggen[i]*(basegenP[i] + alphapart[i]*Delta_scen) + (1 - movinggen[i])*maxpg[i] )   +  (1 - yesnopart[i])*(1 - notoperating[i])*basegenP[i];

# Reactive power disjunctions
subject to vdrop{i in 1..I}:
        busVmagdrop[i] + busVmags[i] >= baseVmags[i];

subject to vrise{i in 1..I}:
        busVmagrise[i] - busVmags[i] >= -baseVmags[i];

subject to Shortage{g in 1..G}:
        shortageQ[g] + genQ[g] = maxqg[g];

subject to SumshortageQ{i in 1..I: gencount[i]>0}:
#subject to SumshortageQ{i in 1..I}:
        - busSumshortageQ[i] + sum{g in gensatbus[i]: gencount[i]>0} shortageQ[g] = 0;

subject to Qshortage0{i in 1..I: gencount[i] > 0 and freeVmags[i] > 0}:
#subject to Qshortage0{i in 1..I}:
        0 <= busSumshortageQ[i] complements busVmagdrop[i] >= 0;

subject to Excess{g in 1..G}:
        excessQ[g] - genQ[g] = - minqg[g];

subject to SumexcessQ{i in 1..I: gencount[i]>0}:
#subject to SumexcessQ{i in 1..I}:
        - busSumexcessQ[i] + sum{g in gensatbus[i]: gencount[i]>0} excessQ[g] = 0;

subject to Qexcess0{i in 1..I: gencount[i] > 0 and freeVmags[i] > 0}:
#subject to Qexcess0{i in 1..I}:
        0 <= busSumexcessQ[i] complements busVmagrise[i] >= 0;

#data;
