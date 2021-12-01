# This is a constrained NLP model formulation used for the contingencies
# for the GO2 Competition from Team ACN.
# 30 Dec 2020
# Daniel Bienstock and Richard Waltz

var count := 0; #dummy counting variable just used in the scripts

param Deltact >=0;
param deltarctg >= 0;
param I >= 0 integer; #numbuses
param G >= 0 integer; #numgens
param L >= 0 integer; #numloads
param NTRcnt >= 0 integer; # 2x number of non-transformer branches
param NTRcnt_original; # should be exactly half of the above
param TRcnt >= 0 integer; # 2x nmber transformer branches
param TRcnt_original; # should be exactly half of the above
param SWSHcnt >= 0 integer; #number of switchable shunts

# param used to exclude dropped generator or line from model
param genout_index integer;
param ntrout_index integer;
param trout_index integer;

#param maxviolBalance >=0; #hard limit on soft constraint violation for balance equations

param numcblocks{1..G} >= 0; #number of cblocks
param pgcostcblock{g in 1..G, n in 1..numcblocks[g]}; # gen block cost coefficient
param pgmaxcblock{g in 1..G, n in 1..numcblocks[g]}; # gen cblock max
param Oncost{g in 1..G};
param Sucost{g in 1..G};
param Sdcost{g in 1..G};
param evlo{1..I};     #emergency lobnd on voltage magnitudes
param evhi{1..I};     #emergency upbnd on voltage magnitudes
param nvalo{1..I};     #lobnd on voltage angles (crude)
param nvahi{1..I};     #upbnd on voltage angles (crude)
set nontransatbus {i in 1..I}; # set of nontrs  connected to each bus
set transatbus {i in 1..I}; # set of trs  connected to each bus
param Gffntr{1..NTRcnt}; # ge
param Gftntr{1..NTRcnt}; # -ge
param Bffntr{1..NTRcnt}; # (be + be^CH / 2)
param Bftntr{1..NTRcnt}; # -be
param ntrbusi{1..NTRcnt}; # i bus
param ntrbusj{1..NTRcnt}; # i bus
param ntrp{1..NTRcnt}; # recall that an ntr is always part of a pair, i.e. for a line (i,j) there is an i-version and a j-version.  this is the count of the common parent


param trbusi{1..TRcnt}; # i bus
param trbusj{1..TRcnt}; # i bus
param tr_o_or_r{1..TRcnt};
param trp{1..TRcnt}; # recall that an tr is always part of a pair, i.e. for a line (i,j) there is an i-version and a j-version.  this is the count of the common parent

param gmf{f in 1..TRcnt_original};
param bmf{f in 1..TRcnt_original};
param barxstf{f in 1..TRcnt_original};
param maxxstf{f in 1..TRcnt_original};
param minxstf{f in 1..TRcnt_original};
param brevetauf{f in 1..TRcnt_original};
param taustf{f in 1..TRcnt_original};
param tau0_fixed{f in 1..TRcnt_original};
param brevethetaf{f in 1..TRcnt_original};
param thetastf{f in 1..TRcnt_original};
param theta0f{f in 1..TRcnt_original};
param Ftauvalid{f in 1..TRcnt_original};
param Fthetavalid{f in 1..TRcnt_original};
param Fnuvalid{f in 1..TRcnt_original};
param numnux{f in 1..TRcnt_original};

param trueg{f in 1..TRcnt_original, m in minxstf[f]..maxxstf[f]};
param trueb{f in 1..TRcnt_original, m in minxstf[f]..maxxstf[f]};

#Info on generators and lines at each bus
param gencount{1..I} integer; # number of gens at each bus
param gencountat0{1..I} integer; # number of gens at each bus that are on
set gensatbus{i in 1..I: gencount[i]>0}; # set of generators at each bus
set gensatbusat0{i in 1..I: gencountat0[i]>0}; # set of generators at each bus that are on
param maxpg{1..G};
param minpg{1..G};
param maxqg{1..G};
param minqg{1..G};

#Info on loads
param loadnumcblocks{1..L} >= 0; #number of load cblocks
param loadcostcblock{l in 1..L, n in 1..loadnumcblocks[l]}; # load block cost coefficients
param loadmaxcblock{l in 1..L, n in 1..loadnumcblocks[l]}; # load cblock max
param loadcount{1..I} integer; # number of loads at each bus
param loadisactive{1..L} integer; #1 or 0 according to active
param activeloadcountat0{1..I} integer; # number of loads at each bus that are on
set loadsatbus{i in 1..I: loadcount[i]>0}; # set of loads at each bus
set activeloadsatbusat0{i in 1..I: activeloadcountat0[i]>0}; # set of loads at each bus that are active


param tmin{l in 1..L}; #= 0 if not active, say with tmax
param tmax{l in 1..L};
param pexist{l in 1..L};  #p0 in documentation, but we are using 0 to refer to base case, so 'exist' refers to the existing amount
param qexist{l in 1..L};  #q0 in documentation, but we are using 0 to refer to base case
param loadP0{l in 1..L};  #final load from base case solution
param loadprdctg{l in 1..L}; #load ramp rate down limit for contingency
param loadpructg{l in 1..L}; #load ramp rate up limit for contingency
param fsglpu{1..I};
param fsblpu{1..I};
set activeswshatbus{i in 1..I}; # set of active switched shunts at each bus
param numactiveswshatbus{i in 1..I}; # number of active switched shunts at each bus
param Pnumcblocks >= 0; #number of bus real power cblocks
param Qnumcblocks >= 0; #number of bus reactive power cblocks
param Snumcblocks >= 0; #number of line/tr apparent power cblocks
param Pcostcblock{n in 1..Pnumcblocks}; # bus block real power cost coefficients
param Qcostcblock{n in 1..Qnumcblocks}; # bus block reactive power cost coefficients
param Scostcblock{n in 1..Snumcblocks}; # line/tr apparent power cost coefficients
param Pmaxcblock{n in 1..Pnumcblocks}; # bus cblock real power max
param Qmaxcblock{n in 1..Qnumcblocks}; # bus cblock reactive power max
param Smaxcblock{n in 1..Snumcblocks}; # line/tr cblock apparent power max

param swshactive{1..SWSHcnt}; #=1 if the switched shunt is active
param bhasize{h in 1..SWSHcnt};
param bha{h in 1..SWSHcnt, a in 1..bhasize[h]}; #  'b' numbers associated with an active shunt
param nha{h in 1..SWSHcnt, a in 1..bhasize[h]}; #  'n' numbers associated with an active shunt

param genonexist{g in 1..G}; #genon0 from base solution: status of generator g at base solution
param genpgexist{g in 1..G}; #genP0 from base solution: status of generator real power pg at base solution
param genprdctg{g in 1..G}; #generator ramp rate down limit for contingencies
param genpructg{g in 1..G}; #generator ramp rate up limit for contingencies
param gensdqualctg{g in 1..G}; #=1 if qualified to shut down
param gensuqualctg{g in 1..G}; #=1 if qualified to start up
param genstartup0{g in 1..G} binary; #>= 0, <= 1;
param genshutdown0{g in 1..G} binary; #>= 0, <= 1;

param ntrqualsw{i in 1..NTRcnt_original}; #=1 if qualified to switch
param ntrcostsw{i in 1..NTRcnt_original}; #cost to switch
param ntrstat{i in 1..NTRcnt_original}; #status of branch in existing solution (ntrON0 from base solution))
param ntrratingct{i in 1..NTRcnt_original}; #current rating
param ntroriginal{i in 1..NTRcnt}; #original or reversed
param trqualsw{i in 1..TRcnt_original}; #=1 if qualified to switch
param trcostsw{i in 1..TRcnt_original}; #cost to switch
param trstat{i in 1..TRcnt_original}; #status of transformer in existing solution (trON0 from base solution))
param trratingct{i in 1..TRcnt_original}; #transformer current rating
param troriginal{i in 1..TRcnt}; #original or reversed

#var xha{h in 1..SWSHcnt, a in 1..bhasize[h]: swshactive[h]*nha[h,a]>0} integer >= 0; #initialize 0???
var xha{h in 1..SWSHcnt, a in 1..bhasize[h]: swshactive[h]*nha[h,a]>0} >= 0; #relax and round???
var bhavar{h in 1..SWSHcnt: swshactive[h]>0}; #initialize 0???

var pgblock0{g in 1..G, n in 1..numcblocks[g]} >= 0, <= pgmaxcblock[g,n];
var gencost0{g in 1..G} >= 0;

var loadblock0{l in 1..L, n in 1..loadnumcblocks[l]} >= 0, <= loadmaxcblock[l,n];  # constraint (14)

var Pplusblock0{i in 1..I, n in 1..Pnumcblocks} >= 0, <= Pmaxcblock[n];  # constraint (4)
var Pminusblock0{i in 1..I, n in 1..Pnumcblocks} >= 0, <= Pmaxcblock[n];  # constraint (6)
var Qplusblock0{i in 1..I, n in 1..Qnumcblocks} >= 0, <= Qmaxcblock[n];  # constraint (9)
var Qminusblock0{i in 1..I, n in 1..Qnumcblocks} >= 0, <= Qmaxcblock[n];  # constraint (10)

#var relaxPplus{i in 1..I}, >=0;
#var relaxPminus{i in 1..I}, >=0;
#var relaxQplus{i in 1..I}, >=0;
#var relaxQminus{i in 1..I}, >=0;
# imposing hard upper bound to avoid very large penalties
var relaxPplus{i in 1..I}, >=0, <= Pcostcblock[Pnumcblocks];
var relaxPminus{i in 1..I}, >=0, <= Pcostcblock[Pnumcblocks];
var relaxQplus{i in 1..I}, >=0, <= Qcostcblock[Qnumcblocks];
var relaxQminus{i in 1..I}, >=0, <= Qcostcblock[Qnumcblocks];

var genon{g in 1..G} binary; # xon_g0, commitment variable;
#var genon{g in 1..G} binary, := genonexist[g]; # initialize to existing?
#var genon{g in 1..G} binary, := 1; # initialize to all ON?
#var genon{g in 1..G} binary, := 0; # initialize to all OFF?

var genstartup{g in 1..G} binary; #>= 0, <= 1;
var genshutdown{g in 1..G} binary; #>= 0, <= 1;
#var genP{g in 1..G} >= minpg[g], <= maxpg[g];
#var genQ{g in 1..G} >= minqg[g], <= maxqg[g];
var genP{g in 1..G};
var genQ{g in 1..G};
var busVmags{i in 1..I}, >=evlo[i], <=evhi[i];     #initialize 1.0 for flat start
var busVangles{i in 1..I}, >=nvalo[i], <=nvahi[i]; #initialize 0.0 for flat start

var ntrON{i in 1..NTRcnt_original} binary; #controls if branch is on or off;
#var ntrON{i in 1..NTRcnt_original} binary, := ntrstat[i]; #initialize to existing?
#var ntrON{i in 1..NTRcnt_original} binary, := 1; #initialize all ON?
#var ntrON{i in 1..NTRcnt_original} binary, := 0; #initialize all OFF?

var trON{f in 1..TRcnt_original} binary; #controls if transformer is on or off; variable (59), x^sw_fk for k = 0
#var trON{f in 1..TRcnt_original} binary, := trstat[f]; #initialize to existing?
#var trON{f in 1..TRcnt_original} binary, := 1; #initialize all ON?
#var trON{f in 1..TRcnt_original} binary, := 0; #initialize all OFF?

#var xstf{f in 1..TRcnt_original} integer; #position of transformer i, i.e. variable (61), x^st_fk, initialize 0?
var xstf{f in 1..TRcnt_original}; #relax and round?
var binxstf{f in 1..TRcnt_original, m in minxstf[f]..maxxstf[f]} >= 0, <= 1; #used to define xstf in parallel (relaxed version)
#var binxstf{f in 1..TRcnt_original, m in minxstf[f]..maxxstf[f]} binary; #>= 0, <= 1; #should be binary

var tauf{f in 1..TRcnt_original}, := 1.0; #variable (62), tau_fk  # how to initialize and keep non-zero?
var thetaf{f in 1..TRcnt_original}; #variable (64), theta_fk
var gf{f in 1..TRcnt_original}; #variable (66) or (68)
#var gf_actual{f in 1..TRcnt}; #like gf0, but accounts for reversal
var Gf{f in 1..TRcnt}; #at origin bus equals gf/tauf^2 + gmf, at destination just gf
var bf{f in 1..TRcnt_original}; #variable (67) or (69)
#var bf_actual{f in 1..TRcnt}; #like bf0, but accounts for reversal
var Bf{f in 1..TRcnt}; #at origin bus equals bf/tauf^2 + bmf, at destination just gf

var ntrP{i in 1..NTRcnt} = ntrON[ntrp[i]]*(Gffntr[i]*busVmags[ntrbusi[i]]^2 +
                      ( Gftntr[i]*cos(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) +
                        Bftntr[i]*sin(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) ) * (busVmags[ntrbusi[i]]*busVmags[ntrbusj[i]])) ; # constraint (50)

var ntrQ{i in 1..NTRcnt} = ntrON[ntrp[i]]*(-Bffntr[i]*busVmags[ntrbusi[i]]^2 +
                      (-Bftntr[i]*cos(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) +
                       Gftntr[i]*sin(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) ) * (busVmags[ntrbusi[i]]*busVmags[ntrbusj[i]])) ;  # constraint (51)


var trP{f in 1..TRcnt} = trON[trp[f]]*(Gf[f]*busVmags[trbusi[f]]^2
                     - ( gf[trp[f]]*cos(busVangles[trbusi[f]] - busVangles[trbusj[f]] - tr_o_or_r[f]*thetaf[trp[f]]) +
                         bf[trp[f]]*sin(busVangles[trbusi[f]] - busVangles[trbusj[f]] - tr_o_or_r[f]*thetaf[trp[f]]) )  * (busVmags[trbusi[f]]*busVmags[trbusj[f]])/tauf[trp[f]] ); # constraints (72,74)


var trQ{f in 1..TRcnt} = trON[trp[f]]*(-Bf[f]*busVmags[trbusi[f]]^2
                     + ( bf[trp[f]]*cos(busVangles[trbusi[f]] - busVangles[trbusj[f]] - tr_o_or_r[f]*thetaf[trp[f]])
                         -gf[trp[f]]*sin(busVangles[trbusi[f]] - busVangles[trbusj[f]] - tr_o_or_r[f]*thetaf[trp[f]]) )  * (busVmags[trbusi[f]]*busVmags[trbusj[f]])/tauf[trp[f]] ); # constraints (73, 75)

var t0{l in 1..L}, >=tmin[l], <=tmax[l]; #constraint (39). Note that As per equation (204), Section B11, active status never changes.  initialize to 1.0?

var loadP{l in 1..L} = pexist[l]*t0[l]; #constraint (40). Recall that t0[l] = 0 if not active
var loadQ{l in 1..L} = qexist[l]*t0[l]; #constraint (41).

var buscost0{i in 1..I} = sum{n in 1..Pnumcblocks} Pcostcblock[n]*(Pplusblock0[i,n] + Pminusblock0[i,n]) +
                          sum{n in 1..Qnumcblocks} Qcostcblock[n]*(Qplusblock0[i,n] + Qminusblock0[i,n]);  # constraint (11), without the delta multiplier

var loadcost0{l in 1..L} = sum{n in 1..loadnumcblocks[l]} loadcostcblock[l,n]*loadblock0[l,n];  # constraint (15), without the delta multiplier

maximize ourobjective:
sum{g in 1..G: g != genout_index} Deltact*(-gencost0[g] - Oncost[g]*genon[g]) - sum{g in 1..G: g != genout_index} (Sucost[g]*genstartup[g] + Sdcost[g]*genshutdown[g]) #total gen cost
- sum{i in 1..I} Deltact*buscost0[i]                                                               #total bus power imbalance cost
+ sum{l in 1..L: loadisactive[l] > 0} Deltact*loadcost0[l]                                         #total load benefit
- sum{e in 1..NTRcnt_original: ntrqualsw[e]>0 and ntrstat[e]=0 and e != ntrout_index} ntrcostsw[e]*ntrON[e]              #penalize line switch ON
- sum{e in 1..NTRcnt_original: ntrqualsw[e]>0 and ntrstat[e]>0 and e != ntrout_index} ntrcostsw[e]*(ntrstat[e]-ntrON[e]) #penalize line switch OFF
- sum{f in 1..TRcnt_original: trqualsw[f]>0 and trstat[f]=0 and f != trout_index} trcostsw[f]*trON[f]                   #penalize transformer switch ON
- sum{f in 1..TRcnt_original: trqualsw[f]>0 and trstat[f]>0 and f != trout_index} trcostsw[f]*(trstat[f]-trON[f]);      #penalize transformer switch OFF


# bus imbalance objective components

subject to totalPplus{i in 1..I}: #constraint (3)
relaxPplus[i] = sum{n in 1..Pnumcblocks} Pplusblock0[i,n];

subject to totalPminus{i in 1..I}: #constraint (5)
relaxPminus[i] = sum{n in 1..Pnumcblocks} Pminusblock0[i,n];

#enforce hard upper bound on violation
#subject to PBalanceViolLimit{i in 1..I}: relaxPplus[i] + relaxPminus[i] <= Pcostcblock[Pnumcblocks];

subject to totalQplus{i in 1..I}: #constraint (7)
relaxQplus[i] = sum{n in 1..Qnumcblocks} Qplusblock0[i,n];

subject to totalQminus{i in 1..I}: #constraint (8)
relaxQminus[i] = sum{n in 1..Qnumcblocks} Qminusblock0[i,n];

#enforce hard upper bound on violation
#subject to QBalanceViolLimit{i in 1..I}: relaxQplus[i] + relaxQminus[i] <= Qcostcblock[Qnumcblocks];

# load benefit objective components

subject to totalload{l in 1..L}: #constraint (13)
loadP[l] = sum{n in 1..loadnumcblocks[l]} loadblock0[l,n];

# generator cost objective components

subject to totalgen{g in 1..G: g != genout_index}: #constraint (27)
genP[g] = sum{n in 1..numcblocks[g]} pgblock0[g,n];

subject to costgen{g in 1..G: g != genout_index}: #first term in (29), without the delta multiplier
gencost0[g] = sum{n in 1..numcblocks[g]} pgcostcblock[g,n]*pgblock0[g,n];

subject to gencommit{g in 1..G: g != genout_index}: #constraint (82)
genon[g] - genonexist[g] = genstartup[g] - genshutdown[g];

# If geonexist=0, force genshutdown=0 (can't shutdown if aready off)
subject to forcegensdoff{g in 1..G: g != genout_index && genonexist[g] == 0}:
genshutdown[g] = 0;

# If geonexist=1, force genstartup=0 (can't startup if aready on)
subject to forcegensuoff{g in 1..G: g != genout_index && genonexist[g] > 0}:
genstartup[g] = 0;

# Fix variables related to genout (probably not necessary, but doesn't hurt)
subject to fixgenoff{g in 1..G: g == genout_index}: #if gen out; force off
genon[g] = 0;
subject to fixgensuoff{g in 1..G: g == genout_index}: #if gen out; force startup off
genstartup[g] = 0;
subject to fixgensdon{g in 1..G: g == genout_index && genonexist[g] > 0}: #if gen out and was on; force shutdown
genshutdown[g] = 1;
subject to fixgensdoff{g in 1..G: g == genout_index && genonexist[g] == 0}: #if gen out and was off; force shutdown off
genshutdown[g] = 0;

subject to gendecide0{g in 1..G: g != genout_index}: #constraint (84)
genstartup[g] + genshutdown[g] <= 1;

#subject to defntrP{i in 1..NTRcnt}:
# ntrP[i] = ntrON[ntrp[i]]*(Gffntr[i]*busVmags[ntrbusi[i]]^2 +
#                      ( Gftntr[i]*cos(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) +
#                        Bftntr[i]*sin(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) ) * (busVmags[ntrbusi[i]]*busVmags[ntrbusj[i]])) ; # constraint (50)

#subject to defntrQ{i in 1..NTRcnt}:
# ntrQ[i] = ntrON[ntrp[i]]*(-Bffntr[i]*busVmags[ntrbusi[i]]^2 +
#                      (-Bftntr[i]*cos(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) +
#                       Gftntr[i]*sin(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) ) * (busVmags[ntrbusi[i]]*busVmags[ntrbusj[i]])) ;  # constraint (51)

#subject to deftrP{f in 1..TRcnt}:
# trP[f] = trON[trp[f]]*(Gf[f]*busVmags[trbusi[f]]^2
#                     - ( gf[trp[f]]*cos(busVangles[trbusi[f]] - busVangles[trbusj[f]] - thetaf[trp[f]]) +
#                         bf[trp[f]]*sin(busVangles[trbusi[f]] - busVangles[trbusj[f]] - thetaf[trp[f]]) )  * (busVmags[trbusi[f]]*busVmags[trbusj[f]])/tauf[trp[f]] ); # constraints (72,74)


#subject to deftrQ{f in 1..TRcnt}:
# trQ[f] = trON[trp[f]]*(-Bf[f]*busVmags[trbusi[f]]^2
#                     + ( bf[trp[f]]*cos(busVangles[trbusi[f]] - busVangles[trbusj[f]] - thetaf[trp[f]]) -
#                         gf[trp[f]]*sin(busVangles[trbusi[f]] - busVangles[trbusj[f]] - thetaf[trp[f]]) )  * (busVmags[trbusi[f]]*busVmags[trbusj[f]])/tauf[trp[f]] ); # constraints (73, 75)

# Active power balance constraints (35)
subject to PBalance{i in 1..I}:
 sum{g in gensatbus[i]: gencount[i]>0 and g != genout_index} genP[g]
 - sum{l in activeloadsatbusat0[i]: activeloadcountat0[i]>0}loadP[l] - fsglpu[i]*busVmags[i]^2  #note: activeloadsatbusat0[i] is the set of STATUS=1 loads at bus i
 -  sum{j in nontransatbus[i]} ntrP[j] - sum{j in transatbus[i]} trP[j] = relaxPplus[i] - relaxPminus[i];

# Reactive power balance constraints (38)
subject to QBalance{i in 1..I}:
 sum{g in gensatbus[i]: gencount[i]>0 and g != genout_index} genQ[g]
 - sum{l in activeloadsatbusat0[i]: activeloadcountat0[i]>0}loadQ[l] + fsblpu[i]*busVmags[i]^2  + (sum{s in activeswshatbus[i]: numactiveswshatbus[i]>0}bhavar[s])*busVmags[i]^2
 -  sum{j in nontransatbus[i]} ntrQ[j] - sum{j in transatbus[i]} trQ[j] = relaxQplus[i] - relaxQminus[i];

# Line current constraints (nontransformers); squared form (55,56)
subject to LineCurrentNTR{i in 1..NTRcnt}:
        ntrP[i]^2 + ntrQ[i]^2 - (ntrratingct[ntrp[i]]*busVmags[ntrbusi[i]])^2 <= 0;

# Line current constraints (transformers); squared form
subject to LineCurrentTR{i in 1..TRcnt}:
        trP[i]^2 + trQ[i]^2 - trratingct[trp[i]]^2 <= 0;

# load ramping constraints (44,45))
subject to ramping{l in 1..L: loadisactive[l]>0}:
loadP0[l] - loadprdctg[l]*deltarctg <= loadP[l] <= loadP0[l] + loadpructg[l]*deltarctg;   #constraints (44-45).  Note that loadP[l] = 0 if load l is not active


subject to fixntrON{i in 1..NTRcnt_original: ntrqualsw[i]=0 and i != ntrout_index}: #constraint (49)
ntrON[i] = ntrstat[i];

subject to fixntrOFF{i in 1..NTRcnt_original: i == ntrout_index}: #if line out; force off
ntrON[i] = 0;

subject to fixtrON{f in 1..TRcnt_original: trqualsw[f]=0 && f != trout_index}: #constraint (60)
trON[f] = trstat[f];

subject to fixtrOFF{i in 1..NTRcnt_original: i == trout_index}: #if transformer line out; force off
trON[i] = 0;

subject to upboundxstf{f in 1..TRcnt_original: f != trout_index}: #constraint (61) upper bound
xstf[f] <= maxxstf[f];

subject to loboundxstf{f in 1..TRcnt_original: f != trout_index}: #constraint (61) lower bound
xstf[f] >= minxstf[f];

#subject to seqxstf{f in 1..TRcnt_original, m in 1..barxstf[f]-1: barxstf[f] - 1 > 0}:
#binxstf[f,m+1] <= binxstf[f,m];

subject to pickxstf{f in 1..TRcnt_original: f != trout_index}:
sum{m in minxstf[f]..maxxstf[f]} binxstf[f,m] = 1;

subject to defxstf{f in 1..TRcnt_original: f != trout_index}:
xstf[f] = sum{m in minxstf[f]..maxxstf[f]} m*binxstf[f,m];

subject to defgf{f in 1..TRcnt_original: Fnuvalid[f] > 0 and f != trout_index}:
  gf[f] = sum{m in minxstf[f]..maxxstf[f]} trueg[f,m]*binxstf[f,m];

subject to fixeddefgf{f in 1..TRcnt_original: Fnuvalid[f] == 0 and f != trout_index}:
  gf[f] = trueg[f,0];

subject to defGfo{f in 1..TRcnt: tr_o_or_r[f] == 1}:
  Gf[f] = gf[trp[f]]/(tauf[trp[f]]^2) +  gmf[trp[f]];

subject to defGfd{f in 1..TRcnt: tr_o_or_r[f] == -1}:
  Gf[f] = gf[trp[f]];

subject to defBfo{f in 1..TRcnt: tr_o_or_r[f] == 1}:
  Bf[f] = bf[trp[f]]/(tauf[trp[f]]^2) +  bmf[trp[f]];

subject to defBfd{f in 1..TRcnt: tr_o_or_r[f] == -1}:
  Bf[f] = bf[trp[f]];

subject to fixeddefbf{f in 1..TRcnt_original: Fnuvalid[f] == 0 and f != trout_index}:
  bf[f] = trueb[f,0];

subject to defbf{f in 1..TRcnt_original: Fnuvalid[f] > 0 and f != trout_index}:
  bf[f] = sum{m in minxstf[f]..maxxstf[f]} trueb[f,m]*binxstf[f,m];

#subject to deftheta0{f in 1..TRcnt_original: Fthetavalid[f] > 0 and f != trout_index}:
#  thetaf[f] = sum{m in minxstf[f]..maxxstf[f]} truetheta[f,m]*binxstf[f,m];



subject to tapratioactive0{f in 1..TRcnt_original: Ftauvalid[f] > 0 and f != trout_index}: #constraint (62)
tauf[f] = brevetauf[f] + taustf[f]*xstf[f];

subject to tapratioinactive0{f in 1..TRcnt_original: Ftauvalid[f] == 0 and f != trout_index}: #constraint (63)
tauf[f] = tau0_fixed[f];

subject to phaseshiftactive0{f in 1..TRcnt_original: Fthetavalid[f] > 0 and f != trout_index}: #constraint (64)
thetaf[f] = brevethetaf[f] + thetastf[f]*xstf[f];

subject to phaseshiftinactive0{f in 1..TRcnt_original: Fthetavalid[f] == 0 and f != trout_index}: #constraint (65)
thetaf[f] = theta0f[f];


subject to upboundxha{h in 1..SWSHcnt, a in 1..bhasize[h]: swshactive[h]*nha[h,a]>0}: #constraint (46)
xha[h,a] <= nha[h,a];

subject to susceptance0{h in 1..SWSHcnt: swshactive[h]>0}: #constraint (47)
bhavar[h] = sum{a in 1..bhasize[h]:nha[h,a]>0} bha[h,a]*xha[h,a];



###### generators on/off, check for 'active' lines or something like that

# WRITING AS A DOUBLE-SIDED CONSTRAINT WAS BUGGY!
#subject to genPbounds{g in 1..G: g != genout_index}: #constraint (85)
#minpg[g]*genon[g] <= genP[g] <= maxpg[g]*genon[g];

subject to genPboundsLO{g in 1..G: g != genout_index}: #constraint (85)
genP[g] >= minpg[g]*genon[g];

subject to genPboundsUP{g in 1..G: g != genout_index}: #constraint (85)
genP[g] <= maxpg[g]*genon[g];

# WRITING AS A DOUBLE-SIDED CONSTRAINT WAS BUGGY!
#subject to genQbounds{g in 1..G: g != genout_index}: #constraint (86)
#minqg[g]*genon[g] <= genQ[g] <= maxqg[g]*genon[g];

subject to genQboundsLO{g in 1..G: g != genout_index}: #constraint (86)
genQ[g] >= minqg[g]*genon[g];

subject to genQboundsUP{g in 1..G: g != genout_index}: #constraint (86)
genQ[g] <= maxqg[g]*genon[g];

subject to genrampingup0{g in 1..G: g != genout_index}: #constraint (89)
genP[g] <= (genpgexist[g] + genpructg[g]*deltarctg)*(genon[g] - genstartup[g]) + (minpg[g] + genpructg[g]*deltarctg)*genstartup[g];

subject to genrampingdown0{g in 1..G: g != genout_index}: #constraint (90)
genP[g] >= (genpgexist[g] - genprdctg[g]*deltarctg)*(genon[g] - genstartup[g]);

subject to gennostartup{g in 1..G: gensuqualctg[g] == 0 and g != genout_index}: #constraint (92)
genstartup[g] = 0;

subject to gennoshutdown{g in 1..G: gensdqualctg[g] == 0 and g != genout_index}: #constraint (94)
genshutdown[g] = 0;

#Generator started in base cannot shutdown
subject to gennoshutdown2{g in 1..G: genstartup0[g] > 0 and g != genout_index}: #constraint (95)
genshutdown[g] = 0;

#Generator shut down in base cannot be started
subject to gennostartup2{g in 1..G: genshutdown0[g] > 0 and g != genout_index}: #constraint (96)
genstartup[g] = 0;



# Make these defined variables
#var ntrP{i in 1..NTRcnt}; #active power flow injected into nontrans i (recall that each line appears twice, once at each end bus)
#var ntrQ{i in 1..NTRcnt}; #reactive power flow injected into nontrans i
#var trP{i in 1..TRcnt};  #active power flow injected into trans i
#var trQ{i in 1..TRcnt};  #reactive power flow injected into trans i
