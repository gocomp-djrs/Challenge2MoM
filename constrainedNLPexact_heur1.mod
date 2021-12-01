# This is a constrained NLP model formulation for the GO2 Competition from Team djrs.
# 28 October 2020
# Daniel Bienstock

#option knitroampl_auxfiles rc;

var count := 0; #dummy counting variable just used in the scripts

param Delta >=0;
param deltar >= 0;
param I >= 0 integer; #numbuses
param G >= 0 integer; #numgens
param L >= 0 integer; #numloads
param NTRcnt >= 0 integer; # 2x number of non-transformer branches
param NTRcnt_original; # should be exactly half of the above
param TRcnt >= 0 integer; # 2x nmber transformer branches
param TRcnt_original; # should be exactly half of the above
param SWSHcnt >= 0 integer; #number of switchable shunts

#param maxviolBalance >=0; #hard limit on soft constraint violation for balance equations

param numcblocks{1..G} >= 0; #number of cblocks in base case
param pgcostcblock{g in 1..G, n in 1..numcblocks[g]}; # gen block cost coefficient base case
param pgmaxcblock{g in 1..G, n in 1..numcblocks[g]}; # gen cblock max  base case
param Oncost{g in 1..G};
param Sucost{g in 1..G};
param Sdcost{g in 1..G};
param nvlo{1..I};     #lobnd on voltage magnitudes
param nvhi{1..I};     #upbnd on voltage magnitudes
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
param heurmaxxstf{f in 1..TRcnt_original};
param heurminxstf{f in 1..TRcnt_original};
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
param gencountat0{1..I} integer; # number of gens at each bus that are on in base case
set gensatbus{i in 1..I: gencount[i]>0}; # set of generators at each bus
set gensatbusat0{i in 1..I: gencountat0[i]>0}; # set of generators at each bus that are on in base case
param maxpg{1..G};
param minpg{1..G};
param maxqg{1..G};
param minqg{1..G};

#Info on loads
param loadnumcblocks{1..L} >= 0; #number of load cblocks in base case
param loadcostcblock{l in 1..L, n in 1..loadnumcblocks[l]}; # load block cost coefficient base case
param loadmaxcblock{l in 1..L, n in 1..loadnumcblocks[l]}; # load cblock max base case
param loadcount{1..I} integer; # number of loads at each bus
param loadisactive{1..L} integer; #1 or 0 according to active
param activeloadcountat0{1..I} integer; # number of loads at each bus that are on in base case
set loadsatbus{i in 1..I: loadcount[i]>0}; # set of loads at each bus
set activeloadsatbusat0{i in 1..I: activeloadcountat0[i]>0}; # set of loads at each bus that are active in base case


param tmin{l in 1..L}; #= 0 if not active, say with tmax
param tmax{l in 1..L};
param pexist{l in 1..L};  #p0 in documentation, but we are using 0 to refer to base case, so 'exist' refers to the existing amount
param qexist{l in 1..L};  #q0 in documentation, but we are using 0 to refer to base case
param prd{l in 1..L}; #load ramp rate down limit
param pru{l in 1..L}; #load ramp rate up limit
param fsglpu{1..I};
param fsblpu{1..I};
set activeswshatbus{i in 1..I}; # set of active switched shunts at each bus
param numactiveswshatbus{i in 1..I}; # number of active switched shunts at each bus
param Pnumcblocks >= 0; #number of bus real power cblocks in base case
param Qnumcblocks >= 0; #number of bus reactive power cblocks in base case
param Snumcblocks >= 0; #number of line/tr apparent power cblocks
param Pcostcblock{n in 1..Pnumcblocks}; # bus block real power cost coefficient base case
param Qcostcblock{n in 1..Qnumcblocks}; # bus block reactive power cost coefficient base case
param Scostcblock{n in 1..Snumcblocks}; # line/tr apparent power cost coefficients
param Pmaxcblock{n in 1..Pnumcblocks}; # bus cblock real power max base case
param Qmaxcblock{n in 1..Qnumcblocks}; # bus cblock reactive power max base case
param Smaxcblock{n in 1..Snumcblocks}; # line/tr cblock apparent power max

param swshactive{1..SWSHcnt}; #=1 if the switched shunt is active
param bhasize{h in 1..SWSHcnt};
param bha{h in 1..SWSHcnt, a in 1..bhasize[h]}; #  'b' numbers associated with an active shunt
param nha{h in 1..SWSHcnt, a in 1..bhasize[h]}; #  'n' numbers associated with an active shunt

param genonexist{g in 1..G}; #what they call xon,0_g, i.e. status of generator g in existing solution
param genpgexist{g in 1..G}; #p_g^0, i.e. status of generator real power pg in existing solution
param genprd{g in 1..G}; #generator ramp rate down limit
param genpru{g in 1..G}; #generator ramp rate up limit
param gensdqual{g in 1..G}; #=1 if qualified to shut down
param gensuqual{g in 1..G}; #=1 if qualified to start up

param ntrqualsw{i in 1..NTRcnt_original}; #=1 if qualified to switch
param ntrcostsw{i in 1..NTRcnt_original}; #cost to switch
param ntrstat{i in 1..NTRcnt_original}; #status of branch in existing solution
param UPntr0{i in 1..NTRcnt_original}; #upper bound on binary variable NTRON0 used to fix it
param LOntr0{i in 1..NTRcnt_original}; #lower bound on binary variable NTRON0 used to fix it

param UPtr0{i in 1..TRcnt_original}; #upper bound on binary variable TRON0 used to fix it
					    param LOtr0{i in 1..TRcnt_original}; #upper bound on binary variable TRON0 used to fix it

param UPgen0{g in 1..G}; #upper bound on binary variable genon0 (used to fix)
param LOgen0{g in 1..G}; #lower bound on binary variable genon0 (used to fix)
param UPgenstartup0{g in 1..G}; #upper bound on binary variable genon0 (used to fix)
param LOgenstartup0{g in 1..G}; #lower bound on binary variable genon0 (used to fix)
param UPgenshutdown0{g in 1..G}; #upper bound on binary variable genon0 (used to fix)
param LOgenshutdown0{g in 1..G}; #lower bound on binary variable genon0 (used to fix)

param UPxha0{h in 1..SWSHcnt, a in 1..bhasize[h]: swshactive[h]*nha[h,a]>0};
param LOxha0{h in 1..SWSHcnt, a in 1..bhasize[h]: swshactive[h]*nha[h,a]>0};
param ntroriginal{i in 1..NTRcnt}; #original or reversed
param ntrrating0{i in 1..NTRcnt}; #current rating
param trqualsw{i in 1..TRcnt_original}; #=1 if qualified to switch
param trcostsw{i in 1..TRcnt_original}; #cost to switch
param trstat{i in 1..TRcnt_original}; #status of transformer in existing solution
param troriginal{i in 1..TRcnt}; #original or reversed
param trrating_original{i in 1..TRcnt_original}; #transformer current rating
param trrating0{i in 1..TRcnt}; #transformer current rating



var xha0{h in 1..SWSHcnt, a in 1..bhasize[h]: swshactive[h]*nha[h,a]>0} integer >= 0; #initialize 0???
#var xha0{h in 1..SWSHcnt, a in 1..bhasize[h]: swshactive[h]*nha[h,a]>0} >= 0; #relax and round???
var bhavar0{h in 1..SWSHcnt: swshactive[h]>0}; #initialize 0???

var pgblock0{g in 1..G, n in 1..numcblocks[g]} >= 0, <= pgmaxcblock[g,n];
var gencost0{g in 1..G} >= 0;

var loadblock0{l in 1..L, n in 1..loadnumcblocks[l]} >= 0, <= loadmaxcblock[l,n];  # constraint (14)

var Pplusblock0{i in 1..I, n in 1..Pnumcblocks} >= 0, <= Pmaxcblock[n];  # constraint (4)
var Pminusblock0{i in 1..I, n in 1..Pnumcblocks} >= 0, <= Pmaxcblock[n];  # constraint (6)
var Qplusblock0{i in 1..I, n in 1..Qnumcblocks} >= 0, <= Qmaxcblock[n];  # constraint (9)
var Qminusblock0{i in 1..I, n in 1..Qnumcblocks} >= 0, <= Qmaxcblock[n];  # constraint (10)

var Strblock{f in 1..TRcnt_original, n in 1..Snumcblocks} >= 0, <= Smaxcblock[n]*trrating_original[f];  # constraint (25)
#var relaxStr{f in 1..TRcnt_original} = sum{n in 1..Snumcblocks} Strblock[f,n]; #constraint (24)
var relaxStr{f in 1..TRcnt_original}, >=0; # slacks for transformer line current constraints
#var relaxStr{f in 1..TRcnt_original}, >=0, <= Scostcblock[Snumcblocks]; # try hard upper bound?

#var relaxPplus0{i in 1..I}, >=0;
#var relaxPminus0{i in 1..I}, >=0;
#var relaxQplus0{i in 1..I}, >=0;
#var relaxQminus0{i in 1..I}, >=0;
# imposing hard upper bound to avoid very large penalties
var relaxPplus0{i in 1..I}, >=0, <= Pcostcblock[Pnumcblocks];
var relaxPminus0{i in 1..I}, >=0, <= Pcostcblock[Pnumcblocks];
var relaxQplus0{i in 1..I}, >=0, <= Qcostcblock[Qnumcblocks];
var relaxQminus0{i in 1..I}, >=0, <= Qcostcblock[Qnumcblocks];

var genon0{g in 1..G} binary; # xon_g0, commitment variable;
#var genon0{g in 1..G} binary, := genonexist[g]; # initialize to existing?
#var genon0{g in 1..G} binary, := 1; # initialize to all ON?
#var genon0{g in 1..G} binary, := 0; # initialize to all OFF?

var genstartup0{g in 1..G} binary; #>= 0, <= 1;
var genshutdown0{g in 1..G} binary; #>= 0, <= 1;
#var genP0{g in 1..G} >= minpg[g], <= maxpg[g];
#var genQ0{g in 1..G} >= minqg[g], <= maxqg[g];
var genP0{g in 1..G};
var genQ0{g in 1..G};
var busVmags{i in 1..I}, >=nvlo[i], <=nvhi[i];     #initialize 1.0 for flat start
var busVangles{i in 1..I}, >=nvalo[i], <=nvahi[i]; #initialize 0.0 for flat start

var ntrON0{i in 1..NTRcnt_original} binary; #controls if branch is on or off;
#var ntrON0{i in 1..NTRcnt_original} binary, := ntrstat[i]; #initialize to existing?
#var ntrON0{i in 1..NTRcnt_original} binary, := 1; #initialize all ON?
#var ntrON0{i in 1..NTRcnt_original} binary, := 0; #initialize all OFF?

var trON0{f in 1..TRcnt_original} binary; #controls if transformer is on or off; variable (59), x^sw_fk for k = 0
#var trON0{f in 1..TRcnt_original} binary, := trstat[f]; #initialize to existing?
#var trON0{f in 1..TRcnt_original} binary, := 1; #initialize all ON?
#var trON0{f in 1..TRcnt_original} binary, := 0; #initialize all OFF?

#var xstf0{f in 1..TRcnt_original} integer; #position of transformer i, i.e. variable (61), x^st_fk for k = 0, initialize 0?
var xstf0{f in 1..TRcnt_original}; #relax and round?
var binxstf0{f in 1..TRcnt_original, m in minxstf[f]..maxxstf[f]} >= 0, <= 1; #used to define xstf0 in parallel (relaxed version)
set fixed0binxstf0{f in 1..TRcnt_original}; # set of m indices such that binstf0[f,m] = 0 is forced

#var binxstf0{f in 1..TRcnt_original, m in minxstf[f]..maxxstf[f]} binary; #>= 0, <= 1; #should be binary

var tauf0{f in 1..TRcnt_original}, := 1.0; #variable (62), tau_fk for k = 0  # how to initialize and keep non-zero?
var thetaf0{f in 1..TRcnt_original}; #variable (64), theta_fk for k = 0
var gf0{f in 1..TRcnt_original}; #variable (66) or (68) for k = 0
var gf0_actual{f in 1..TRcnt}; #like gf0, but accounts for reversal
var Gf0{f in 1..TRcnt}; #at origin bus equals gf0/tauf0^2 + gmf, at destination just gf0
var bf0{f in 1..TRcnt_original}; #variable (67) or (69) for k = 0
var bf0_actual{f in 1..TRcnt}; #like bf0, but accounts for reversal
var Bf0{f in 1..TRcnt}; #at origin bus equals bf0/tauf0^2 + bmf, at destination just gf0

var ntrP0{i in 1..NTRcnt} = ntrON0[ntrp[i]]*(Gffntr[i]*busVmags[ntrbusi[i]]^2 +
                      ( Gftntr[i]*cos(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) +
                        Bftntr[i]*sin(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) ) * (busVmags[ntrbusi[i]]*busVmags[ntrbusj[i]])) ; # constraint (50)

var ntrQ0{i in 1..NTRcnt} = ntrON0[ntrp[i]]*(-Bffntr[i]*busVmags[ntrbusi[i]]^2 +
                      (-Bftntr[i]*cos(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) +
                       Gftntr[i]*sin(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) ) * (busVmags[ntrbusi[i]]*busVmags[ntrbusj[i]])) ;  # constraint (51)


var trP0{f in 1..TRcnt} = trON0[trp[f]]*(Gf0[f]*busVmags[trbusi[f]]^2
                     - (  gf0[trp[f]]*cos(busVangles[trbusi[f]] - busVangles[trbusj[f]] - tr_o_or_r[f]*thetaf0[trp[f]]) +
                           bf0[trp[f]]*sin(busVangles[trbusi[f]] - busVangles[trbusj[f]] - tr_o_or_r[f]*thetaf0[trp[f]]) )  * (busVmags[trbusi[f]]*busVmags[trbusj[f]])/tauf0[trp[f]] ); # constraints (72,74)


var trQ0{f in 1..TRcnt} = trON0[trp[f]]*(-Bf0[f]*busVmags[trbusi[f]]^2
                     + ( bf0[trp[f]]*cos(busVangles[trbusi[f]] - busVangles[trbusj[f]] -  tr_o_or_r[f]*thetaf0[trp[f]])
                         -gf0[trp[f]]*sin(busVangles[trbusi[f]] - busVangles[trbusj[f]]  - tr_o_or_r[f]*thetaf0[trp[f]]) )  * (busVmags[trbusi[f]]*busVmags[trbusj[f]])/tauf0[trp[f]] ); # constraints (73, 75)

var t0{l in 1..L}, >=tmin[l], <=tmax[l]; #constraint (39). Note that As per equation (204), Section B11, active status never changes.  initialize to 1.0?

var loadP0{l in 1..L} = pexist[l]*t0[l]; #constraint (40). Recall that t0[l] = 0 if not active
var loadQ0{l in 1..L} = qexist[l]*t0[l]; #constraint (41).


var buscost0{i in 1..I} = sum{n in 1..Pnumcblocks} Pcostcblock[n]*(Pplusblock0[i,n] + Pminusblock0[i,n]) +
                          sum{n in 1..Qnumcblocks} Qcostcblock[n]*(Qplusblock0[i,n] + Qminusblock0[i,n]);  # constraint (11), without the delta multiplier
var busPcost0{i in 1..I} = sum{n in 1..Pnumcblocks} Delta*Pcostcblock[n]*(Pplusblock0[i,n] + Pminusblock0[i,n]);
var totalbusPcost0 = sum{i in 1..I}busPcost0[i];
var busQcost0{i in 1..I} = sum{n in 1..Qnumcblocks} Delta*Qcostcblock[n]*(Qplusblock0[i,n] + Qminusblock0[i,n]);
#var busQcost0{i in 1..I} >=0;
var totalbusQcost0 = sum{i in 1..I}busQcost0[i];

var loadcost0{l in 1..L} = sum{n in 1..loadnumcblocks[l]} loadcostcblock[l,n]*loadblock0[l,n];  # constraint (15), without the delta multiplier
var Qfsh{i in 1..I} = -fsblpu[i]*busVmags[i]^2;
var Qswsh{i in 1..I} = -(sum{s in activeswshatbus[i]: numactiveswshatbus[i]>0}bhavar0[s])*busVmags[i]^2;
var fullgencost{g in 1..G} = Delta*(gencost0[g] + Oncost[g]*genon0[g]) + (Sucost[g]*genstartup0[g] + Sdcost[g]*genshutdown0[g]); #total gen cost at 0
var totalgencost = 	sum{g in 1..G} fullgencost[g];
var totalbuscost = totalbusPcost0 + totalbusQcost0;
var genPatbus{i in 1..I} = sum{g in gensatbus[i]: gencount[i]>0} genP0[g];
var genQatbus{i in 1..I} = sum{g in gensatbus[i]: gencount[i]>0} genQ0[g];
var loadQatbus{i in 1..I} = sum{l in activeloadsatbusat0[i]: activeloadcountat0[i]>0}loadQ0[l];

var trexceed{f in 1..TRcnt_original} = sum{n in 1..Snumcblocks} Scostcblock[n]*Strblock[f,n];   #2nd part of (22) w/out delta multiplier

maximize ourobjective:
-totalgencost - totalbuscost
+ sum{l in 1..L: loadisactive[l] > 0} Delta*loadcost0[l]                                         #total load benefit at 0
- sum{e in 1..NTRcnt_original: ntrqualsw[e]>0 and ntrstat[e]=0} ntrcostsw[e]*ntrON0[e]              #penalize line switch ON
- sum{e in 1..NTRcnt_original: ntrqualsw[e]>0 and ntrstat[e]>0} ntrcostsw[e]*(ntrstat[e]-ntrON0[e]) #penalize line switch OFF
- sum{f in 1..TRcnt_original: trqualsw[f]>0 and trstat[f]=0} trcostsw[f]*trON0[f]                   #penalize transformer switch ON
- sum{f in 1..TRcnt_original: trqualsw[f]>0 and trstat[f]>0} trcostsw[f]*(trstat[f]-trON0[f])       #penalize transformer switch OFF
- sum{f in 1..TRcnt_original} (Delta*trexceed[f]);           #penalize transformer exceedance

#- sum{e in 1..NTRcnt_original: ntrstat[e]=0} ntrcostsw[e]*ntrON0[e]               #penalize line switch ON
#- sum{e in 1..NTRcnt_original: ntrstat[e]>0} ntrcostsw[e]*(ntrstat[e]-ntrON0[e])  #penalize line switch OFF
#- sum{f in 1..TRcnt_original: trstat[f]=0} trcostsw[f]*trON0[f]                   #penalize transformer switch ON
#- sum{f in 1..TRcnt_original: trstat[f]>0} trcostsw[f]*(trstat[f]-trON0[f]);      #penalize transformer switch OFF

# bus imbalance objective components

subject to totalPplus0{i in 1..I}: #constraint (3)
relaxPplus0[i] = sum{n in 1..Pnumcblocks} Pplusblock0[i,n];

subject to totalPminus0{i in 1..I}: #constraint (5)
relaxPminus0[i] = sum{n in 1..Pnumcblocks} Pminusblock0[i,n];

#enforce hard upper bound on violation
#subject to PBalanceViolLimit{i in 1..I}: relaxPplus0[i] + relaxPminus0[i] <= Pcostcblock[Pnumcblocks];

subject to totalQplus0{i in 1..I}: #constraint (7)
relaxQplus0[i] = sum{n in 1..Qnumcblocks} Qplusblock0[i,n];

subject to totalQminus0{i in 1..I}: #constraint (8)
relaxQminus0[i] = sum{n in 1..Qnumcblocks} Qminusblock0[i,n];

#subject to busQcostdef{i in 1..I}:
#busQcost0[i] = Delta*(sum{n in 1..Qnumcblocks} Qcostcblock[n]*(Qplusblock0[i,n] + Qminusblock0[i,n]));

#enforce hard upper bound on violation
#subject to QBalanceViolLimit{i in 1..I}: relaxQplus0[i] + relaxQminus0[i] <= Qcostcblock[Qnumcblocks];

# load benefit objective components

subject to totalload0{l in 1..L}: #constraint (13)
loadP0[l] = sum{n in 1..loadnumcblocks[l]} loadblock0[l,n];

# generator cost objective components

subject to totalgen0{g in 1..G}: #constraint (27)
genP0[g] = sum{n in 1..numcblocks[g]} pgblock0[g,n];

subject to costgen0{g in 1..G}: #first term in (29), without the delta multiplier
gencost0[g] = sum{n in 1..numcblocks[g]} pgcostcblock[g,n]*pgblock0[g,n];

subject to gencommit0{g in 1..G}: #constraint (82)
genon0[g] - genonexist[g] = genstartup0[g] - genshutdown0[g];

subject to gendecide0{g in 1..G}: #constraint (84)
genstartup0[g] + genshutdown0[g] <= 1;

# If geonexist=0, force genshutdown=0 (can't shutdown if aready off)
subject to forcegensdoff0{g in 1..G: genonexist[g] == 0}:
genshutdown0[g] = 0;

# If geonexist=1, force genstartup=0 (can't startup if aready on)
subject to forcegensuoff0{g in 1..G: genonexist[g] > 0}:
genstartup0[g] = 0;

#subject to defntrP0{i in 1..NTRcnt}:
# ntrP0[i] = ntrON0[ntrp[i]]*(Gffntr[i]*busVmags[ntrbusi[i]]^2 +
#                      ( Gftntr[i]*cos(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) +
#                        Bftntr[i]*sin(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) ) * (busVmags[ntrbusi[i]]*busVmags[ntrbusj[i]])) ; # constraint (50)

#subject to defntrQ0{i in 1..NTRcnt}:
# ntrQ0[i] = ntrON0[ntrp[i]]*(-Bffntr[i]*busVmags[ntrbusi[i]]^2 +
#                      (-Bftntr[i]*cos(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) +
#                       Gftntr[i]*sin(busVangles[ntrbusi[i]] - busVangles[ntrbusj[i]]) ) * (busVmags[ntrbusi[i]]*busVmags[ntrbusj[i]])) ;  # constraint (51)

#subject to deftrP0{f in 1..TRcnt}:
# trP0[f] = trON0[trp[f]]*(Gf0[f]*busVmags[trbusi[f]]^2
#                     - ( gf0[trp[f]]*cos(busVangles[trbusi[f]] - busVangles[trbusj[f]] - thetaf0[trp[f]]) +
#                         bf0[trp[f]]*sin(busVangles[trbusi[f]] - busVangles[trbusj[f]] - thetaf0[trp[f]]) )  * (busVmags[trbusi[f]]*busVmags[trbusj[f]])/tauf0[trp[f]] ); # constraints (72,74)


#subject to deftrQ0{f in 1..TRcnt}:
# trQ0[f] = trON0[trp[f]]*(-Bf0[f]*busVmags[trbusi[f]]^2
#                     + ( bf0[trp[f]]*cos(busVangles[trbusi[f]] - busVangles[trbusj[f]] - thetaf0[trp[f]]) -
#                         gf0[trp[f]]*sin(busVangles[trbusi[f]] - busVangles[trbusj[f]] - thetaf0[trp[f]]) )  * (busVmags[trbusi[f]]*busVmags[trbusj[f]])/tauf0[trp[f]] ); # constraints (73, 75)

# Active power balance constraints (35)
subject to PBalanceat0{i in 1..I}:
genPatbus[i]
 - sum{l in activeloadsatbusat0[i]: activeloadcountat0[i]>0}loadP0[l] - fsglpu[i]*busVmags[i]^2  #note: activeloadsatbusat0[i] is the set of STATUS=1 loads at bus i
 -  sum{j in nontransatbus[i]} ntrP0[j] - sum{j in transatbus[i]} trP0[j] = relaxPplus0[i] - relaxPminus0[i];

# Reactive power balance constraints (38)
subject to QBalanceat0{i in 1..I}:
genQatbus[i]
 - loadQatbus[i] - Qfsh[i] - Qswsh[i]
 -  sum{j in nontransatbus[i]} ntrQ0[j] - sum{j in transatbus[i]} trQ0[j] = relaxQplus0[i] - relaxQminus0[i];

# Line current constraints (nontransformers); squared form (55,56)
subject to LineCurrentNTR{i in 1..NTRcnt}:
        ntrP0[i]^2 + ntrQ0[i]^2 - (ntrrating0[i]*busVmags[ntrbusi[i]])^2 <= 0;

# Line current constraints (transformers); squared form (77,78)
subject to LineCurrentTR{i in 1..TRcnt}:
        trP0[i]^2 + trQ0[i]^2 <= (trrating_original[trp[i]] + relaxStr[trp[i]])^2;

subject to totalStr{i in 1..TRcnt_original}: #constraint (24)
relaxStr[i] = sum{n in 1..Snumcblocks} Strblock[i,n];

#subject to relaxStrpos{i in 1..TRcnt_original}: #relaxation variable non-negative
#relaxStr[i] >= 0;

subject to fixStr{i in 1..TRcnt_original}: #avoid largest penalty
if (Snumcblocks>1 and Scostcblock[Snumcblocks]>=1e6) then Strblock[i,Snumcblocks] = 0;


subject to ramping0{l in 1..L: loadisactive[l]>0}:
pexist[l] - prd[l]*deltar <= loadP0[l] <= pexist[l] + pru[l]*deltar;   #constraints (42-43).  Note that loadP0[l] = 0 if load l is not active


subject to fixntrON0{i in 1..NTRcnt_original: ntrqualsw[i]=0}: #constraint (49)
ntrON0[i] = ntrstat[i];

subject to fixtrON0{f in 1..TRcnt_original: trqualsw[f]=0}: #constraint (60)
trON0[f] = trstat[f];

subject to upboundxstf0{f in 1..TRcnt_original}: #constraint (61) upper bound
xstf0[f] <= heurmaxxstf[f];

subject to loboundxstf0{f in 1..TRcnt_original}: #constraint (61) lower bound
xstf0[f] >= heurminxstf[f];

#subject to seqxstf0{f in 1..TRcnt_original, m in 1..barxstf[f]-1: barxstf[f] - 1 > 0}:
#binxstf0[f,m+1] <= binxstf0[f,m];

subject to pickxstf0{f in 1..TRcnt_original}:
sum{m in minxstf[f]..maxxstf[f]} binxstf0[f,m] = 1;

subject to defxstf0{f in 1..TRcnt_original}:
xstf0[f] = sum{m in minxstf[f]..maxxstf[f]} m*binxstf0[f,m];

subject to defgf0{f in 1..TRcnt_original: Fnuvalid[f] > 0}:
  gf0[f] = sum{m in minxstf[f]..maxxstf[f]} trueg[f,m]*binxstf0[f,m];

subject to fixeddefgf0{f in 1..TRcnt_original: Fnuvalid[f] == 0}:
  gf0[f] = trueg[f,0];

subject to defGf0o{f in 1..TRcnt: tr_o_or_r[f] == 1}:
  Gf0[f] = gf0[trp[f]]/(tauf0[trp[f]]^2) +  gmf[trp[f]];

subject to defGf0d{f in 1..TRcnt: tr_o_or_r[f] == -1}:
  Gf0[f] = gf0[trp[f]];

subject to defBf0o{f in 1..TRcnt: tr_o_or_r[f] == 1}:
  Bf0[f] = bf0[trp[f]]/(tauf0[trp[f]]^2) +  bmf[trp[f]];

subject to defBf0d{f in 1..TRcnt: tr_o_or_r[f] == -1}:
  Bf0[f] = bf0[trp[f]];

subject to fixeddefbf0{f in 1..TRcnt_original: Fnuvalid[f] == 0}:
  bf0[f] = trueb[f,0];

subject to defbf0{f in 1..TRcnt_original: Fnuvalid[f] > 0}:
  bf0[f] = sum{m in minxstf[f]..maxxstf[f]} trueb[f,m]*binxstf0[f,m];

#subject to deftheta0{f in 1..TRcnt_original: Fthetavalid[f] > 0}:
#  thetaf0[f] = sum{m in minxstf[f]..maxxstf[f]} truetheta[f,m]*binxstf0[f,m];



subject to tapratioactive0{f in 1..TRcnt_original: Ftauvalid[f] > 0}: #constraint (62)
tauf0[f] = brevetauf[f] + taustf[f]*xstf0[f];

subject to tapratioinactive0{f in 1..TRcnt_original: Ftauvalid[f] == 0}: #constraint (63)
tauf0[f] = tau0_fixed[f];

subject to phaseshiftactive0{f in 1..TRcnt_original: Fthetavalid[f] > 0}: #constraint (64)
thetaf0[f] = brevethetaf[f] + thetastf[f]*xstf0[f];

subject to phaseshiftinactive0{f in 1..TRcnt_original: Fthetavalid[f] == 0}: #constraint (65)
thetaf0[f] = theta0f[f];


subject to upboundxha0{h in 1..SWSHcnt, a in 1..bhasize[h]: swshactive[h]*nha[h,a]>0}: #constraint (46)
xha0[h,a] <= nha[h,a];

subject to susceptance0{h in 1..SWSHcnt: swshactive[h]>0}: #constraint (47)
bhavar0[h] = sum{a in 1..bhasize[h]:nha[h,a]>0} bha[h,a]*xha0[h,a];



###### generators on/off, check for 'active' lines or something like that

# WRITING AS A DOUBLE-SIDED CONSTRAINT WAS BUGGY!
#subject to genP0bounds{g in 1..G}: #constraint (85)
#minpg[g]*genon0[g] <= genP0[g] <= maxpg[g]*genon0[g];

subject to genP0boundsLO{g in 1..G}: #constraint (85)
genP0[g] >= minpg[g]*genon0[g];

subject to genP0boundsUP{g in 1..G}: #constraint (85)
genP0[g] <= maxpg[g]*genon0[g];

# WRITING AS A DOUBLE-SIDED CONSTRAINT WAS BUGGY!
#subject to genQ0bounds{g in 1..G}: #constraint (86)
#minqg[g]*genon0[g] <= genQ0[g] <= maxqg[g]*genon0[g];

subject to genQ0boundsLO{g in 1..G}: #constraint (86)
genQ0[g] >= minqg[g]*genon0[g];

subject to genQ0boundsUP{g in 1..G}: #constraint (86)
genQ0[g] <= maxqg[g]*genon0[g];

subject to genrampingup0{g in 1..G}: #constraint (87)
genP0[g] <= (genpgexist[g] + genpru[g]*deltar)*(genon0[g] - genstartup0[g]) + (minpg[g] + genpru[g]*deltar)*genstartup0[g];

subject to genrampingdown0{g in 1..G}: #constraint (88)
genP0[g] >= (genpgexist[g] - genprd[g]*deltar)*(genon0[g] - genstartup0[g]);

subject to gennostartup0{g in 1..G: gensuqual[g] == 0}: #constraint (91)
genstartup0[g] = 0;

subject to gennoshutdown0{g in 1..G: gensdqual[g] == 0}: #constraint (93)
genshutdown0[g] = 0;

subject to fixthebins{f in 1..TRcnt_original}:
sum{m in fixed0binxstf0[f]} binxstf0[f,m] = 0;

# Make these defined variables
#var ntrP0{i in 1..NTRcnt}; #active power flow injected into nontrans i (recall that each line appears twice, once at each end bus)
#var ntrQ0{i in 1..NTRcnt}; #reactive power flow injected into nontrans i
#var trP0{i in 1..TRcnt};  #active power flow injected into trans i
#var trQ0{i in 1..TRcnt};  #reactive power flow injected into trans i


#sum{g in 1..G} Delta*(-gencost0[g] - Oncost[g]*genon0[g]) - sum{g in 1..G} (Sucost[g]*genstartup0[g] + Sdcost[g]*genshutdown0[g]) #total gen cost at 0

subject to upbndntrON0{i in 1..NTRcnt_original}:
ntrON0[i] <= UPntr0[i];

subject to lobndntrON0{i in 1..NTRcnt_original}:
ntrON0[i] >= LOntr0[i];

subject to upbndtrON0{i in 1..TRcnt_original}:
trON0[i] <= UPtr0[i];

subject to lobndtrON0{i in 1..TRcnt_original}:
trON0[i] >= LOtr0[i];

subject to upbndgenON0{g in 1..G}:
genon0[g] <= UPgen0[g];

subject to lobndgenON0{g in 1..G}:
genon0[g] >= LOgen0[g];

subject to upbndgenstartup0{g in 1..G}:
genstartup0[g] <= UPgenstartup0[g];

subject to lobndgenstartup0{g in 1..G}:
genstartup0[g] >= LOgenstartup0[g];

subject to upbndgenshutdown0{g in 1..G}:
genshutdown0[g] <= UPgenshutdown0[g];

subject to lobndgenshutdown0{g in 1..G}:
genshutdown0[g] >= LOgenshutdown0[g];




subject to upxha0{h in 1..SWSHcnt, a in 1..bhasize[h]: swshactive[h]*nha[h,a]>0}:
xha0[h,a] <= UPxha0[h,a];

subject to loxha0{h in 1..SWSHcnt, a in 1..bhasize[h]: swshactive[h]*nha[h,a]>0}:
xha0[h,a] >= LOxha0[h,a];