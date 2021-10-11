function generators_subproblem(sol::PFsol,dualsol::Dualsol)::Tuple{Dict{Tuple{String,Int64},Float64},PFsol}
	m=Model(optimizer_with_attributes(Mosek.Optimizer,
	"MSK_IPAR_PRESOLVE_USE"=>1,
	"MSK_DPAR_INTPNT_CO_TOL_DFEAS" => 1e-7,
	"MSK_DPAR_INTPNT_TOL_DFEAS" => 1e-7,
	"MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-7,
	"MSK_DPAR_INTPNT_TOL_PFEAS" => 1e-7,
	"MSK_DPAR_INTPNT_CO_TOL_MU_RED" => 1e-7,
	"MSK_DPAR_INTPNT_TOL_MU_RED" => 1e-7,
	"MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1e-7,
	"MSK_DPAR_INTPNT_TOL_REL_GAP" => 1e-7,
	"MSK_DPAR_MIO_TOL_ABS_RELAX_INT" => 1e-7,
	"MSK_DPAR_MIO_TOL_FEAS" => 1e-7,
	"MSK_DPAR_INTPNT_CO_TOL_NEAR_REL" => 1,
	))
	@variable(m, pg[seg_bid])
	@variable(m, -1e3 <= qg[dist_seg_bid] <= 1e3)
	@expression(m, net_pg_bid[bd=bid], sum(pg[sb] for sb=bid_2_seg_bids[bd]))
	@expression(m, net_qg_bid[bd=dist_bid], sum(qg[sb] for sb=bid_2_seg_bids[bd]))
	### Binary variables ###
	@variable(m, sa[seg_bid], Bin)
	@variable(m, qa[bid], Bin)
	@variable(m, qta[qt_bid], Bin)
	@variable(m, alpha[bid], Bin)
	@variable(m, omega[bid], Bin)

	@constraint(m, PDown[sb=seg_bid], -pg[sb] <= -sa[sb]*p_min[sb])
	@constraint(m, PUp[sb=seg_bid], pg[sb] <= sa[sb]*p_max[sb])

	@constraint(m, ExQtBidsConstr[i=ex_qt_bid_id],
	sum(qta[ExQt] for ExQt=ex_qt_bid[i]) <= 1)

	@constraint(m, QtaDefLow[qt=qt_bid],
	qta[qt] <= sum(qa[bd] for bd=bid if bd[2]==qt))

	@constraint(m, QtaDefHigh[bd=bid],
	qa[bd] <= qta[bd[2]])

	@constraint(m, DefQa[sb=seg_bid],
	sa[sb] <= qa[(sb[1],sb[2],sb[3],sb[5])])

	@constraint(m, RampMinIncr[rcI=RampIndexMinIncr],
	+ sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
	- sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
	- get(RCRealPowerIncr,(rcI[1],rcI[2],rcI[4]),0)
	>= 0)

	@constraint(m, RampMaxIncr[rcI=RampIndexMaxIncr],
	+ sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
	+ get(RCRealPowerIncr,(rcI[1],rcI[2],rcI[4]),0)
	- sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
	>= 0)

	@constraint(m, DiscQtoP[BQpd=QPDiscConstrIndex],
	[get(QPDiscMax,(BQpd[2],BQpd[3],BQpd[5]),0),
	net_pg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])],
	net_qg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]]
	in SecondOrderCone()
	)

	@constraint(m, ConstrainToUpFromLineQtoP[BHpc=QPHPConstrIndex1],
	net_qg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	-QPHPOffset[(BHpc[5],BHpc[2])]
	-QPHPSlope[(BHpc[5],BHpc[2])]*net_pg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	>= 0)

	@constraint(m, ConstrainToDownFromLineQtoP[BHpc=QPHPConstrIndex0],
	QPHPOffset[(BHpc[5],BHpc[2])]
	+QPHPSlope[(BHpc[5],BHpc[2])]*net_pg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	-net_qg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])] >= 0)

	@constraint(m, AlphaOmegaDef[B=AlphaOmegaLate],
	qa[B]
	- sum(qa[bd] for bd=bid if bd[1]==B[1] && bd[1]==B[2] && bd[4]==B[4]-1)
	- alpha[B]
	+ omega[B]
	== 0)

	@constraint(m, MutexAlphaOmega[B=AlphaOmegaIndex],
	alpha[B] + omega[B] <= 1)

	@constraint(m, AlphaInit[B=AlphaOmegaEarly],
	qa[B] - alpha[B] == 0)

	@constraint(m, OmegaInit[B=AlphaOmegaEarly],
	omega[B] == 0)

	@constraint(m, MinDurationConstr[mdp=MinDurIndex],
	+ sum(qa[bd] for bd=bid if bd[2]==mdp[1] && bd[3]==mdp[3] && bd[4]==mdp[5])
 	- sum(alpha[bd] for bd=bid if bd[2]==mdp[1] && bd[3]==mdp[2] && bd[4]==mdp[4])
	>= 0)

	@constraint(m, NoNewAct[B=NoNewActIndex],
	alpha[B] == 0)

	# @expression(m, obj, +sum((dualsol.priceP[sb[1],sb[5]]-cost[sb][1])*pg[sb] for sb=seg_bid)
	# +sum(dualsol.priceQ[sb[1],sb[5]]*qg[sb] for sb=dist_seg_bid)
	# -sum((dualsol.priceP[sb[1],sb[5]]-cost[sb][1])*sol.pg[sb] for sb=seg_bid)
	# -sum(dualsol.priceQ[sb[1],sb[5]]*sol.qg[sb] for sb=dist_seg_bid))
	#
	# @constraint(m, ObjPos, obj >= 0)

	### OBJECTIVE ###
	@objective(m, Max, +sum((dualsol.priceP[sb[1],sb[5]]-cost[sb][1])*pg[sb] for sb=seg_bid)
	+sum(dualsol.priceQ[sb[1],sb[5]]*qg[sb] for sb=dist_seg_bid)
	-sum((dualsol.priceP[sb[1],sb[5]]-cost[sb][1])*sol.pg[sb] for sb=seg_bid)
	-sum(dualsol.priceQ[sb[1],sb[5]]*sol.qg[sb] for sb=dist_seg_bid)
	)

	# for sb=seg_bid if abs(dualsol.priceP[sb[1],sb[5]]-cost[sb][1]) < 1e-6
	# 	#fix(pg[sb],sol.pg[sb];force=true)
	# 	@constraint(m, pg[sb]==sol.pg[sb])
	# end end
	# for sb=dist_seg_bid if abs(dualsol.priceQ[sb[1],sb[5]]) < 1e-6
	# 	#fix(qg[sb],sol.qg[sb];force=true)
	# 	@constraint(m, qg[sb]==sol.qg[sb])
	# end end

	#println(m)
	optimize!(m)
	println(termination_status(m))
	println(primal_status(m))

	pgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(sb=> value(pg[sb]) for sb=seg_bid)
	qgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(sb=> value(qg[sb]) for sb=dist_seg_bid)

	locg=Dict{Tuple{String,Int64},Float64}((j,t) =>
		+mapreduce(sum,+,(dualsol.priceP[sb[1],sb[5]]-cost[sb][1])*pgHist[sb] for sb=seg_bid if sb[1] in DNs[j] && t==sb[5];init=0)
		+mapreduce(sum,+,dualsol.priceQ[sb[1],sb[5]]*qgHist[sb] for sb=dist_seg_bid if sb[1] in DNs[j] && t==sb[5];init=0)
		-mapreduce(sum,+,(dualsol.priceP[sb[1],sb[5]]-cost[sb][1])*sol.pg[sb] for sb=seg_bid if sb[1] in DNs[j] && t==sb[5];init=0)
		-mapreduce(sum,+,dualsol.priceQ[sb[1],sb[5]]*sol.qg[sb] for sb=dist_seg_bid if sb[1] in DNs[j] && t==sb[5];init=0)
	for j=keys(DNs) for t=times)
	merge!(locg,Dict{Tuple{String,Int64},Float64}(("TN",t)=>
	+sum((dualsol.priceP[sb[1],sb[5]]-cost[sb][1])*pgHist[sb] for sb=tr_seg_bid if t==sb[5])
	-sum((dualsol.priceP[sb[1],sb[5]]-cost[sb][1])*sol.pg[sb] for sb=tr_seg_bid if t==sb[5]) for t=times))
	return (locg,
	PFsol(false,pgHist,qgHist,sol.c,sol.s,sol.theta,sol.fp,sol.fq,sol.fpTo,sol.fqTo,sol.slackPUp,sol.slackPDown,sol.slackQUp,sol.slackQDown))
end

function gen_subproblem(node::Int64,sol::PFsol,dualsol::Dualsol)
	m=Model(optimizer_with_attributes(Mosek.Optimizer,
	"MSK_IPAR_PRESOLVE_USE"=>1,"MSK_IPAR_LOG" => 0,
	"MSK_DPAR_INTPNT_CO_TOL_DFEAS" => 1e-7,
	"MSK_DPAR_INTPNT_TOL_DFEAS" => 1e-7,
	"MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-7,
	"MSK_DPAR_INTPNT_TOL_PFEAS" => 1e-7,
	"MSK_DPAR_INTPNT_CO_TOL_MU_RED" => 1e-7,
	"MSK_DPAR_INTPNT_TOL_MU_RED" => 1e-7,
	"MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1e-7,
	"MSK_DPAR_INTPNT_TOL_REL_GAP" => 1e-7,
	"MSK_DPAR_MIO_TOL_ABS_RELAX_INT" => 1e-7,
	"MSK_DPAR_MIO_TOL_FEAS" => 1e-7,
	"MSK_DPAR_INTPNT_CO_TOL_NEAR_REL" => 1,
	))
	seg_bid_i=[sb for sb=seg_bid if sb[1]==node]
	bid_i=[bd for bd=bid if bd[1]==node]
	qt_bid_i=Set([bd[2] for bd=bid_i])
	@variable(m, pg[seg_bid_i])
	@variable(m, -1e3 <= qg[seg_bid_i] <= 1e3)
	@expression(m, net_pg_bid[bd=bid_i], sum(pg[sb] for sb=bid_2_seg_bids[bd]))
	@expression(m, net_qg_bid[bd=bid_i], sum(qg[sb] for sb=bid_2_seg_bids[bd]))
	### Binary variables ###
	@variable(m, sa[seg_bid_i], Bin)
	@variable(m, qa[bid_i], Bin)
	@variable(m, qta[qt_bid_i], Bin)
	@variable(m, alpha[bid_i], Bin)
	@variable(m, omega[bid_i], Bin)

	@constraint(m, PDown[sb=seg_bid_i], -pg[sb] <= -sa[sb]*p_min[sb])
	@constraint(m, PUp[sb=seg_bid_i], pg[sb] <= sa[sb]*p_max[sb])

	@constraint(m, ExQtBidsConstr[i=ex_qt_bid_id],
	sum(qta[ExQt] for ExQt=ex_qt_bid[i] if ExQt in qt_bid_i) <= 1)

	@constraint(m, QtaDefLow[qt=qt_bid_i],
	qta[qt] <= sum(qa[bd] for bd=bid_i if bd[2]==qt))

	@constraint(m, QtaDefHigh[bd=bid_i],
	qa[bd] <= qta[bd[2]])

	@constraint(m, DefQa[sb=seg_bid_i],
	sa[sb] <= qa[(sb[1],sb[2],sb[3],sb[5])])

	@constraint(m, RampMinIncr[rcI=[rimi for rimi=RampIndexMinIncr if rimi[1] in qt_bid_i]],
	+ sum(net_pg_bid[bd] for bd=bid_i if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
	- sum(net_pg_bid[bd] for bd=bid_i if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
	- get(RCRealPowerIncr,(rcI[1],rcI[2],rcI[4]),0)
	>= 0)

	@constraint(m, RampMaxIncr[rcI=[rimi for rimi=RampIndexMaxIncr if rimi[1] in qt_bid_i]],
	+ sum(net_pg_bid[bd] for bd=bid_i if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
	+ get(RCRealPowerIncr,(rcI[1],rcI[2],rcI[4]),0)
	- sum(net_pg_bid[bd] for bd=bid_i if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
	>= 0)

	@constraint(m, DiscQtoP[BQpd=[qpdc for qpdc=QPDiscConstrIndex if qpdc[1]==node]],
	[get(QPDiscMax,(BQpd[2],BQpd[3],BQpd[5]),0),
	net_pg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])],
	net_qg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]]
	in SecondOrderCone()
	)

	@constraint(m, ConstrainToUpFromLineQtoP[BHpc=[qphp for qphp=QPHPConstrIndex1 if qphp[1]==node]],
	net_qg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	-QPHPOffset[(BHpc[5],BHpc[2])]
	-QPHPSlope[(BHpc[5],BHpc[2])]*net_pg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	>= 0)

	@constraint(m, ConstrainToDownFromLineQtoP[BHpc=[qphp for qphp=QPHPConstrIndex0 if qphp[1]==node]],
	QPHPOffset[(BHpc[5],BHpc[2])]
	+QPHPSlope[(BHpc[5],BHpc[2])]*net_pg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	-net_qg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])] >= 0)

	@constraint(m, AlphaOmegaDef[B=[bd for bd=AlphaOmegaLate if bd in bid_i]],
	qa[B]
	- sum(qa[bd] for bd=bid_i if bd[1]==B[1] && bd[1]==B[2] && bd[4]==B[4]-1)
	- alpha[B]
	+ omega[B]
	== 0)

	@constraint(m, MutexAlphaOmega[B=[bd for bd=AlphaOmegaIndex if bd[1]==node]],
	alpha[B] + omega[B] <= 1)

	@constraint(m, AlphaInit[B=[bd for bd=AlphaOmegaEarly if bd[1]==node]],
	qa[B] - alpha[B] == 0)

	@constraint(m, OmegaInit[B=[bd for bd=AlphaOmegaEarly if bd[1]==node]],
	omega[B] == 0)

	@constraint(m, MinDurationConstr[mdp=[mdi for mdi=MinDurIndex if mdi[1] in qt_bid_i]],
	+ sum(qa[bd] for bd=bid_i if bd[2]==mdp[1] && bd[3]==mdp[3] && bd[4]==mdp[5])
 	- sum(alpha[bd] for bd=bid_i if bd[2]==mdp[1] && bd[3]==mdp[2] && bd[4]==mdp[4])
	>= 0)

	@constraint(m, NoNewAct[B=[nnai for nnai=NoNewActIndex if nnai[1]==node]],
	alpha[B] == 0)

	# @expression(m, obj, +sum((dualsol.priceP[sb[1],sb[5]]-cost[sb][1])*pg[sb] for sb=seg_bid)
	# +sum(dualsol.priceQ[sb[1],sb[5]]*qg[sb] for sb=dist_seg_bid)
	# -sum((dualsol.priceP[sb[1],sb[5]]-cost[sb][1])*sol.pg[sb] for sb=seg_bid)
	# -sum(dualsol.priceQ[sb[1],sb[5]]*sol.qg[sb] for sb=dist_seg_bid))
	#
	# @constraint(m, ObjPos, obj >= 0)

	### OBJECTIVE ###
	@objective(m, Max, +sum((dualsol.priceP[sb[1],sb[5]]-cost[sb][1])*pg[sb] for sb=seg_bid_i)
	+sum(dualsol.priceQ[sb[1],sb[5]]*qg[sb] for sb=seg_bid_i if node in dist_nodes)
	-sum((dualsol.priceP[sb[1],sb[5]]-cost[sb][1])*sol.pg[sb] for sb=seg_bid_i)
	-sum(dualsol.priceQ[sb[1],sb[5]]*sol.qg[sb] for sb=seg_bid_i if node in dist_nodes)
	)

	# for sb=seg_bid if abs(dualsol.priceP[sb[1],sb[5]]-cost[sb][1]) < 1e-6
	# 	#fix(pg[sb],sol.pg[sb];force=true)
	# 	@constraint(m, pg[sb]==sol.pg[sb])
	# end end
	# for sb=dist_seg_bid if abs(dualsol.priceQ[sb[1],sb[5]]) < 1e-6
	# 	#fix(qg[sb],sol.qg[sb];force=true)
	# 	@constraint(m, qg[sb]==sol.qg[sb])
	# end end

	#println(m)
	optimize!(m)
	#println(termination_status(m))
	#println(primal_status(m))
	@assert(termination_status(m)==MOI.OPTIMAL)
	return objective_value(m)
end

function locg_decentralized(sol::PFsol,dualsol::Dualsol)
	locg=Dict{Int64,Float64}()
	counting=1
	for i=nodes
		locg[i]=gen_subproblem(i,sol,dualsol)
		#println(">$counting. SOLVED node $i")
		counting+=1
	end
	println("\n>>> LOCG = $(sum(values(locg)))")
	return locg
end

function operator_subproblem_relax(sol::PFsol,dualsol::Dualsol,preprocess=0)::Tuple{Dict{Int64,Float64},PFsol}
	m=Model(optimizer_with_attributes(Mosek.Optimizer,"MSK_IPAR_PRESOLVE_USE"=>1,
	"MSK_DPAR_INTPNT_CO_TOL_DFEAS" => 1e-7,
	"MSK_DPAR_INTPNT_TOL_DFEAS" => 1e-7,
	"MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-7,
	"MSK_DPAR_INTPNT_TOL_PFEAS" => 1e-7,
	"MSK_DPAR_INTPNT_CO_TOL_MU_RED" => 1e-7,
	"MSK_DPAR_INTPNT_TOL_MU_RED" => 1e-7,
	"MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1e-7,
	"MSK_DPAR_INTPNT_TOL_REL_GAP" => 1e-7,
	"MSK_DPAR_MIO_TOL_ABS_RELAX_INT" => 1e-7,
	"MSK_DPAR_MIO_TOL_FEAS" => 1e-7,
	"MSK_DPAR_MIO_TOL_REL_GAP" => 1e-2,
	"MSK_DPAR_INTPNT_CO_TOL_NEAR_REL" => 1,
	))

	### VARIABLES ###
	@variable(m, 0 <= c[c_index,times])
	@variable(m, s[s_index,times])
	@variable(m, -pi <= theta[nodes,times] <= pi)

	### FLOW DEFINITION ###
	# @variable(m, fp[e=edges,t=times])
	# @variable(m, fpTo[e=edges,t=times])
	# @variable(m, fq[e=dist_edges,t=times])
	# @variable(m, fqTo[e=dist_edges,t=times])
	max_diff=1e3
	@variable(m, sol.fp[e,t] - max_diff <= fp[e=edges,t=times] <= sol.fp[e,t] + max_diff)
	@variable(m, sol.fpTo[e,t] - max_diff <= fpTo[e=edges,t=times] <= sol.fpTo[e,t] + max_diff)
	@variable(m, sol.fq[e,t] - max_diff <= fq[e=dist_edges,t=times] <= sol.fq[e,t] + max_diff)
	@variable(m, sol.fqTo[e,t] - max_diff <= fqTo[e=dist_edges,t=times] <= sol.fqTo[e,t] + max_diff)

	@constraint(m, DefFpTr[e=tr_edges,t=times],
	X_pu[e]*fp[e,t]== (theta[e_fr[e],t]-theta[e_to[e],t]))
	@constraint(m, DefFpX0[e=tr_edges,t=times],
	fp[e,t] + fpTo[e,t] == 0)

	DefFpDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	DefFpToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	for l=dist_edges for t=times
		if R_pu[l] == 0.0
			get!(DefFpDist,(l,t), @constraint(m, fp[l,t] + fpTo[l,t] == 0))
			@constraint(m, c[(e_fr[l],e_to[l]),t] ==  c[(e_to[l],e_to[l]),t])
			@constraint(m, c[(e_fr[l],e_fr[l]),t] ==  c[(e_fr[l],e_to[l]),t])
			@constraint(m, s[(e_fr[l],e_to[l]),t] == 0.0)
			@constraint(m, theta[e_fr[l],t] == theta[e_to[l],t])
		else
			get!(DefFpDist,(l,t),
			@constraint(m,
				fp[l,t] == -Gbus[e_fr[l],e_to[l]]*c[(e_fr[l],e_fr[l]),t]
				+ Gbus[e_fr[l],e_to[l]]*c[(e_fr[l],e_to[l]),t]
				- Bbus[e_fr[l],e_to[l]]*s[(e_fr[l],e_to[l]),t])
			)
			get!(DefFpToDist,(l,t),
			@constraint(m,
				# fpTo[l,t] == -Gbus[e_to[l],e_fr[l]]*c[(e_to[l],e_to[l]),t]
				# + Gbus[e_to[l],e_fr[l]]*c[(e_fr[l],e_to[l]),t]
				# + Bbus[e_to[l],e_fr[l]]*s[(e_fr[l],e_to[l]),t])
				fpTo[l,t]
				+ fp[l,t]
				+ Gbus[e_to[l],e_fr[l]]*(c[(e_fr[l],e_fr[l]),t] - 2*c[(e_fr[l],e_to[l]),t] + c[(e_to[l],e_to[l]),t])
				== 0)
			)
		end
	end end

	DefFqDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	DefFqToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	for l=dist_edges for t=times
		if X_pu[l] == 0.0
			get!(DefFqDist,(l,t), @constraint(m, fq[l,t] + fqTo[l,t] == 0))
		else
			get!(DefFqDist,(l,t),
			@constraint(m,
				fq[l,t] == Bbus[e_fr[l],e_to[l]]*c[(e_fr[l],e_fr[l]),t]
				- Bbus[e_fr[l],e_to[l]]*c[(e_fr[l],e_to[l]),t]
				- Gbus[e_fr[l],e_to[l]]*s[(e_fr[l],e_to[l]),t])
			)
			get!(DefFqToDist,(l,t),
			@constraint(m,
				# fqTo[l,t] == Bbus[e_to[l],e_fr[l]]*c[(e_to[l],e_to[l]),t]
				# - Bbus[e_to[l],e_fr[l]]*c[(e_fr[l],e_to[l]),t]
				# + Gbus[e_to[l],e_fr[l]]*s[(e_fr[l],e_to[l]),t])
				fq[l,t]
				+ fqTo[l,t]
				- Bbus[e_to[l],e_fr[l]]*(c[(e_fr[l],e_fr[l]),t] - 2*c[(e_fr[l],e_to[l]),t] + c[(e_to[l],e_to[l]),t])
				== 0)
			)
		end
	end end

	### CONSTRAINTS ###
	@constraint(m, TrPowerFlowLimitUp[e=intersect(line_limit,tr_edges),t=times],
	fp[e,t] <= rate_a[e])

	@constraint(m, TrPowerFlowLimitDown[e=intersect(line_limit,tr_edges),t=times],
	-fp[e,t] <= rate_a[e])

	@constraint(m, VoltageLimitsDown[i=v_limit,t=times],
	-c[(i,i),t] <= -vsq_min[i])

	@constraint(m, VoltageLimitsUp[i=v_limit,t=times],
	c[(i,i),t] <= vsq_max[i])

	@constraint(m, LinePowerLimit1[e=intersect(line_limit,dist_edges),t=times],
	[rate_a[e],fp[e,t],fq[e,t]] in SecondOrderCone())

	@constraint(m, LinePowerLimit2[e=intersect(line_limit,dist_edges),t=times],
	[rate_a[e],fpTo[e,t],fqTo[e,t]] in SecondOrderCone())

	mod_s_ind=[(e_fr[e],e_to[e]) for e=dist_edges if R_pu[e]!=0.0]
	@constraint(m, SOCPconstraint[e=mod_s_ind,t=times],
	[0.5*c[(e[1],e[1]),t],c[(e[2],e[2]),t],c[e,t],s[e,t]] in RotatedSecondOrderCone())

	for i=slack_bus for t=times @constraint(m,theta[i,t]==0) end end

	if preprocess!=0
		for e=preprocess for t=times
			@constraint(m, fp[e,t]==sol.fp[e,t])
			@constraint(m, fpTo[e,t]==sol.fpTo[e,t])
			@constraint(m, fq[e,t]==sol.fq[e,t])
			@constraint(m, fqTo[e,t]==sol.fqTo[e,t])
			@constraint(m, c[(e_fr[e],e_to[e]),t]==sol.c[(e_fr[e],e_to[e]),t])
			@constraint(m, s[(e_fr[e],e_to[e]),t]==sol.s[(e_fr[e],e_to[e]),t])
		end end
	end

	@expression(m, obj,
	+ sum(-dualsol.priceP[e_fr[e],t]*(fp[e,t]-sol.fp[e,t])
	- dualsol.priceP[e_to[e],t]*(fpTo[e,t]-sol.fpTo[e,t]) for e=edges for t=times)
	+ sum(-dualsol.priceQ[e_fr[e],t]*(fq[e,t]-sol.fq[e,t])
	- dualsol.priceQ[e_to[e],t]*(fqTo[e,t]-sol.fqTo[e,t]) for e=dist_edges for t=times)
	)

	### OBJECTIVE ###
	@objective(m, Max, obj)

	optimize!(m)
	println(primal_status(m))
	println(termination_status(m))

	cHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((i,t) => value(c[i,t]) for i=c_index for t=times)
	sHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((l,t) => value(s[l,t]) for l=s_index for t=times)
	thetaHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(theta[i,t]) for i=nodes for t=times)
	fpHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fp[l,t]) for l=edges for t=times)
	fqHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fq[l,t]) for l=dist_edges for t=times)
	fpToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fpTo[l,t]) for l=edges for t=times)
	fqToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fqTo[l,t]) for l=dist_edges for t=times)
	return (Dict{Int64,Float64}(t =>
	+ mapreduce(sum,+,init=0,-dualsol.priceP[e_fr[e],t]*(value(fp[e,t])-sol.fp[e,t])
	- dualsol.priceP[e_to[e],t]*(value(fpTo[e,t])-sol.fpTo[e,t]) for e=edges)
	+ mapreduce(sum,+,init=0,-dualsol.priceQ[e_fr[e],t]*(value(fq[e,t])-sol.fq[e,t])
	- dualsol.priceQ[e_to[e],t]*(value(fqTo[e,t])-sol.fqTo[e,t]) for e=dist_edges) for t=times),
	PFsol(false,sol.pg,sol.qg,cHist,sHist,thetaHist,fpHist,fqHist,fpToHist,fqToHist,sol.slackPUp,sol.slackPDown,sol.slackQUp,sol.slackQDown))
end

function operator_subproblem(sol::PFsol,dualsol::Dualsol)::Tuple{Dict{Int64,Float64},PFsol}
 	m=Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>"yes",
	"bound_frac" => 1e-5,
	"bound_push" => 1e-5,
	"slack_bound_frac" => 1e-5,
	"slack_bound_push" => 1e-5,
	"bound_mult_init_method" => "mu-based",
	"mu_init" => 1e-4,
	))
	### VARIABLES ###
	@variable(m, 0 <= c[c_index,times])
	@variable(m, s[s_index,times])
	@variable(m, -pi <= theta[nodes,times] <= pi)

	### FLOW DEFINITION ###
	@variable(m, fp[edges,times])
	@constraint(m, DefFpTr[e=tr_edges,t=times],
	X_pu[e]*fp[e,t]== (theta[e_fr[e],t]-theta[e_to[e],t]))
	@variable(m,fpTo[edges,times])
	@constraint(m, DefFpX0[e=tr_edges,t=times],
	fp[e,t] + fpTo[e,t] == 0)

	DefFpDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	DefFpToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	for l=dist_edges for t=times
		if R_pu[l] == 0.0
			get!(DefFpDist,(l,t), @constraint(m, fp[l,t] + fpTo[l,t] == 0))
			@constraint(m, c[(e_fr[l],e_fr[l]),t] ==  c[(e_to[l],e_to[l]),t])
			@constraint(m, c[(e_fr[l],e_fr[l]),t] ==  c[(e_fr[l],e_to[l]),t])
			@constraint(m, s[(e_fr[l],e_to[l]),t] == 0.0)
			@constraint(m, theta[e_fr[l]] == theta[e_to[l]])
		else
			get!(DefFpDist,(l,t),
			@constraint(m,
				fp[l,t] == -Gbus[e_fr[l],e_to[l]]*c[(e_fr[l],e_fr[l]),t]
				+ Gbus[e_fr[l],e_to[l]]*c[(e_fr[l],e_to[l]),t]
				- Bbus[e_fr[l],e_to[l]]*s[(e_fr[l],e_to[l]),t])
			)
			get!(DefFpToDist,(l,t),
			@constraint(m,
				fpTo[l,t]
				+ fp[l,t]
				+ Gbus[e_to[l],e_fr[l]]*(c[(e_fr[l],e_fr[l]),t] - 2*c[(e_fr[l],e_to[l]),t] + c[(e_to[l],e_to[l]),t])
				== 0
			))
		end
	end end
	@variable(m, fq[dist_edges,times])
	@variable(m, fqTo[dist_edges,times])
	DefFqDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	DefFqToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	for l=dist_edges for t=times
		if X_pu[l] == 0.0
			get!(DefFqDist,(l,t), @constraint(m, fq[l,t] + fqTo[l,t] == 0))
		else
			get!(DefFqDist,(l,t),
			@constraint(m,
				fq[l,t] == Bbus[e_fr[l],e_to[l]]*c[(e_fr[l],e_fr[l]),t]
				- Bbus[e_fr[l],e_to[l]]*c[(e_fr[l],e_to[l]),t]
				- Gbus[e_fr[l],e_to[l]]*s[(e_fr[l],e_to[l]),t])
			)
			get!(DefFqToDist,(l,t),
			@constraint(m,
				fq[l,t]
				+ fqTo[l,t]
				- Bbus[e_to[l],e_fr[l]]*(c[(e_fr[l],e_fr[l]),t] - 2*c[(e_fr[l],e_to[l]),t] + c[(e_to[l],e_to[l]),t])
				== 0
			))
		end
	end end

	### CONSTRAINTS ###
	@constraint(m, TrPowerFlowLimitUp[e=intersect(line_limit,tr_edges),t=times],
	fp[e,t] <= rate_a[e])

	@constraint(m, TrPowerFlowLimitDown[e=intersect(line_limit,tr_edges),t=times],
	-fp[e,t] <= rate_a[e])

	@constraint(m, VoltageLimitsDown[i=v_limit,t=times],
	-c[(i,i),t] <= -vsq_min[i])

	@constraint(m, VoltageLimitsUp[i=v_limit,t=times],
	c[(i,i),t] <= vsq_max[i])

	@constraint(m, LinePowerLimit1[e=intersect(line_limit,dist_edges),t=times],
	fp[e,t]^2 + fq[e,t]^2 <= rate_a[e]^2)

	@constraint(m, LinePowerLimit2[e=intersect(line_limit,dist_edges),t=times],
	fpTo[e,t]^2 + fqTo[e,t]^2 <= rate_a[e]^2)

	mod_s_ind=[(e_fr[e],e_to[e]) for e=dist_edges if R_pu[e]!=0.0]
	@NLconstraint(m, Quadconstraint[e=s_index,t=times],
	c[e,t]^2 + s[e,t]^2 == c[(e[1],e[1]),t]*c[(e[2],e[2]),t])

	@NLconstraint(m, Tanconstraint[e=s_index,t=times],
	c[e,t]*sin(theta[e[1],t]-theta[e[2],t])
	+ s[e,t]*cos(theta[e[1],t]-theta[e[2],t]) == 0)

	for i=slack_bus for t=times fix(theta[i,t],0;force=true) end end

	@expression(m, obj,
	+sum(-dualsol.priceP[i,t]*(sum(fpTo[e,t] for e=n_to[i])
		+ sum(fp[e,t] for e=n_fr[i])) for i=nodes for t=times)
	+sum(-dualsol.priceQ[i,t]*(sum(fqTo[e,t] for e=dist_n_to[i])
		+ sum(fq[e,t] for e=dist_n_fr[i])) for i=dist_nodes for t=times)
	-sum(-dualsol.priceP[i,t]*(sum(sol.fpTo[e,t] for e=n_to[i])
		+ sum(sol.fp[e,t] for e=n_fr[i])) for i=nodes for t=times)
	-sum(-dualsol.priceQ[i,t]*(sum(sol.fqTo[e,t] for e=dist_n_to[i])
		+ sum(sol.fq[e,t] for e=dist_n_fr[i])) for i=dist_nodes for t=times)
	)

	#@constraint(m, ObjPos, obj >= 0)
	#ext_edges=[e for e=tr_edges if ((e_fr[e] in tr_nodes && e_to[e] in dist_nodes) || (e_to[e] in tr_nodes && e_fr[e] in dist_nodes))]
	#for e=ext_edges t=times @constraint(fp[e,t]==sol.fp[e,t]) end end
	### OBJECTIVE ###
	@objective(m, Max, obj)

	for t=times
		for l=c_index
			set_start_value(c[l,t],sol.c[l,t])
		end
		for i=nodes
			set_start_value(theta[i,t],sol.theta[i,t])
		end
		for l=s_index
			set_start_value(s[l,t],sol.s[l,t])
		end
		for l=edges
			set_start_value(fp[l,t],sol.fp[l,t])
			set_start_value(fpTo[l,t],sol.fpTo[l,t])
		end
		for l=dist_edges
			set_start_value(fq[l,t],sol.fq[l,t])
			set_start_value(fqTo[l,t],sol.fqTo[l,t])
		end
	end

	#println(m)
	optimize!(m)

	cHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((i,t) => value(c[i,t]) for i=c_index for t=times)
	sHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((l,t) => value(s[l,t]) for l=s_index for t=times)
	thetaHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(theta[i,t]) for i=nodes for t=times)
	fpHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fp[l,t]) for l=edges for t=times)
	fqHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fq[l,t]) for l=dist_edges for t=times)
	fpToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fpTo[l,t]) for l=edges for t=times)
	fqToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fqTo[l,t]) for l=dist_edges for t=times)
	return (Dict{Int64,Float64}(t =>
	+mapreduce(sum,+,-dualsol.priceP[i,t]*(mapreduce(sum,+,fpToHist[e,t] for e=n_to[i];init=0)
		+ mapreduce(sum,+,fpHist[e,t] for e=n_fr[i];init=0)) for i=nodes;init=0)
	+mapreduce(sum,+,-dualsol.priceQ[i,t]*(mapreduce(sum,+,fqToHist[e,t] for e=dist_n_to[i];init=0)
		+ mapreduce(sum,+,fqHist[e,t] for e=dist_n_fr[i];init=0)) for i=dist_nodes;init=0)
	-mapreduce(sum,+,-dualsol.priceP[i,t]*(mapreduce(sum,+,sol.fpTo[e,t] for e=n_to[i];init=0)
		+ mapreduce(sum,+,sol.fp[e,t] for e=n_fr[i];init=0)) for i=nodes;init=0)
	-mapreduce(sum,+,-dualsol.priceQ[i,t]*(mapreduce(sum,+,sol.fqTo[e,t] for e=dist_n_to[i];init=0)
		+ mapreduce(sum,+,sol.fq[e,t] for e=dist_n_fr[i];init=0)) for i=dist_nodes;init=0) for t=times),
	PFsol(false,sol.pg,sol.qg,cHist,sHist,thetaHist,fpHist,fqHist,fpToHist,fqToHist,sol.slackPUp,sol.slackPDown,sol.slackQUp,sol.slackQDown))
end

function operator_subproblem_dec(sol::PFsol,dualsol::Dualsol)::Tuple{Dict{Int64,Float64},PFsol}
	cHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}()
	sHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}()
	thetaHist = Dict{Tuple{Int64,Int64},Float64}()
	fpHist = Dict{Tuple{Int64,Int64},Float64}()
	fqHist = Dict{Tuple{Int64,Int64},Float64}()
	fpToHist = Dict{Tuple{Int64,Int64},Float64}()
	fqToHist = Dict{Tuple{Int64,Int64},Float64}()
	obj_times= Dict{Int64,Float64}()
	tol_var=.1
	for t=times
 	m=Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>"yes",
	#"bound_frac" => 1e-10,
	#"bound_push" => 1e-10,
	# "slack_bound_frac" => 1e-10,
	# "slack_bound_push" => 1e-10,
	#"bound_mult_init_method" => "mu-based",
	# "mu_init" => 1e-5,
	))
	### VARIABLES ###
	@variable(m, 0 <= c[i=c_index,t])
	@variable(m, s[i=s_index,t])
	@variable(m, -pi <= theta[i=nodes,t] <= pi)
	@variable(m, sol.fp[e,t] - tol_var <= fp[e=edges,t] <= sol.fp[e,t] + tol_var )
	@variable(m, sol.fpTo[e,t] - tol_var <= fpTo[e=edges,t] <= sol.fpTo[e,t] + tol_var)
	@variable(m, sol.fq[e,t] - tol_var <= fq[e=dist_edges,t] <= sol.fq[e,t] + tol_var)
	@variable(m, sol.fqTo[e,t] - tol_var <= fqTo[e=dist_edges,t] <= sol.fqTo[e,t] + tol_var)

	### FLOW DEFINITION ###
	@constraint(m, DefFpTr[e=tr_edges],
	#fp[e,t]== Bbus[e_fr[e],e_to[e]]*(theta_tr[e_fr[e],t]-theta_tr[e_to[e],t]))
	X_pu[e]*fp[e,t]== (theta[e_fr[e],t]-theta[e_to[e],t]))
	# @constraint(m, DefFpToTr[e=tr_edges,t=times],
	# #fpTo[e,t]== Bbus[e_to[e],e_fr[e]]*(theta_tr[e_to[e],t]-theta_tr[e_fr[e],t]))
	# X_pu[e]*fpTo[e,t]== (theta_tr[e_to[e],t]-theta_tr[e_fr[e],t]))
	#
	# @constraint(m, DefFpX0[e=[l for l=tr_edges if X_pu[l]==0.0],t=times],
	# fp[e,t] + fpTo[e,t] == 0)
	@constraint(m, DefFpX0[e=tr_edges],
	fp[e,t] + fpTo[e,t] == 0)

	DefFpDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	DefFpToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	for l=dist_edges
		if R_pu[l] == 0.0
			get!(DefFpDist,(l,t), @constraint(m, fp[l,t] + fpTo[l,t] == 0))
			@constraint(m, c[(e_fr[l],e_to[l]),t] ==  c[(e_to[l],e_to[l]),t])
			@constraint(m, c[(e_fr[l],e_fr[l]),t] ==  c[(e_fr[l],e_to[l]),t])
			@constraint(m, s[(e_fr[l],e_to[l]),t] == 0.0)
			@constraint(m, theta[e_fr[l]] == theta[e_to[l]])
		else
			get!(DefFpDist,(l,t),
			@constraint(m,
				fp[l,t] == -Gbus[e_fr[l],e_to[l]]*c[(e_fr[l],e_fr[l]),t]
				+ Gbus[e_fr[l],e_to[l]]*c[(e_fr[l],e_to[l]),t]
				- Bbus[e_fr[l],e_to[l]]*s[(e_fr[l],e_to[l]),t])
			)
			get!(DefFpToDist,(l,t),
			@constraint(m,
				fpTo[l,t] == -Gbus[e_to[l],e_fr[l]]*c[(e_to[l],e_to[l]),t]
				+ Gbus[e_to[l],e_fr[l]]*c[(e_fr[l],e_to[l]),t]
				+ Bbus[e_to[l],e_fr[l]]*s[(e_fr[l],e_to[l]),t]
			))
		end
	end
	DefFqDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	DefFqToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	for l=dist_edges
		if X_pu[l] == 0.0
			get!(DefFqDist,(l,t), @constraint(m, fq[l,t] + fqTo[l,t] == 0))
		else
			get!(DefFqDist,(l,t),
			@constraint(m,
				fq[l,t] == Bbus[e_fr[l],e_to[l]]*c[(e_fr[l],e_fr[l]),t]
				- Bbus[e_fr[l],e_to[l]]*c[(e_fr[l],e_to[l]),t]
				- Gbus[e_fr[l],e_to[l]]*s[(e_fr[l],e_to[l]),t])
			)
			get!(DefFqToDist,(l,t),
			@constraint(m,
				fqTo[l,t] == Bbus[e_to[l],e_fr[l]]*c[(e_to[l],e_to[l]),t]
				- Bbus[e_to[l],e_fr[l]]*c[(e_fr[l],e_to[l]),t]
				+ Gbus[e_to[l],e_fr[l]]*s[(e_fr[l],e_to[l]),t]
			))
		end
	end

	### CONSTRAINTS ###
	@constraint(m, TrPowerFlowLimitUp[e=intersect(line_limit,tr_edges)],
	fp[e,t] <= rate_a[e])

	@constraint(m, TrPowerFlowLimitDown[e=intersect(line_limit,tr_edges)],
	-fp[e,t] <= rate_a[e])

	@constraint(m, VoltageLimitsDown[i=v_limit],
	-c[(i,i),t] <= -vsq_min[i])

	@constraint(m, VoltageLimitsUp[i=v_limit],
	c[(i,i),t] <= vsq_max[i])

	@constraint(m, LinePowerLimit1[e=intersect(line_limit,dist_edges)],
	fp[e,t]^2 + fq[e,t]^2 <= rate_a[e]^2)

	@constraint(m, LinePowerLimit2[e=intersect(line_limit,dist_edges)],
	fpTo[e,t]^2 + fqTo[e,t]^2 <= rate_a[e]^2)

	mod_s_ind=[(e_fr[e],e_to[e]) for e=dist_edges if R_pu[e]!=0.0]
	@NLconstraint(m, Quadconstraint[e=mod_s_ind],
	c[e,t]^2 + s[e,t]^2 == c[(e[1],e[1]),t]*c[(e[2],e[2]),t])

	@NLconstraint(m, Tanconstraint[e=mod_s_ind],
	c[e,t]*sin(theta[e[1],t]-theta[e[2],t])
	+ s[e,t]*cos(theta[e[1],t]-theta[e[2],t]) == 0)

	for i=slack_bus @constraint(m,theta[i,t]==0) end

	@expression(m, obj,
	+sum(-dualsol.priceP[i,t]*(sum(fpTo[e,t] for e=n_to[i])
		+ sum(fp[e,t] for e=n_fr[i])) for i=nodes )
	+sum(-dualsol.priceQ[i,t]*(sum(fqTo[e,t] for e=dist_n_to[i])
		+ sum(fq[e,t] for e=dist_n_fr[i])) for i=dist_nodes )
	-sum(-dualsol.priceP[i,t]*(sum(sol.fpTo[e,t] for e=n_to[i])
		+ sum(sol.fp[e,t] for e=n_fr[i])) for i=nodes )
	-sum(-dualsol.priceQ[i,t]*(sum(sol.fqTo[e,t] for e=dist_n_to[i])
		+ sum(sol.fq[e,t] for e=dist_n_fr[i])) for i=dist_nodes)
	)

	#@constraint(m, ObjPos, obj >= 0)
	#ext_edges=[e for e=tr_edges if ((e_fr[e] in tr_nodes && e_to[e] in dist_nodes) || (e_to[e] in tr_nodes && e_fr[e] in dist_nodes))]
	#for e=ext_edges t=times @constraint(fp[e,t]==sol.fp[e,t]) end end
	### OBJECTIVE ###
	@objective(m, Max, obj)

	for l=c_index
		set_start_value(c[l,t],sol.c[l,t])
	end
	for i=nodes
		set_start_value(theta[i,t],sol.theta[i,t])
	end
	for l=s_index
		set_start_value(s[l,t],sol.s[l,t])
	end
	for l=edges
		set_start_value(fp[l,t],sol.fp[l,t])
		set_start_value(fpTo[l,t],sol.fpTo[l,t])
	end
	for l=dist_edges
		set_start_value(fq[l,t],sol.fq[l,t])
		set_start_value(fqTo[l,t],sol.fqTo[l,t])
	end

	#println(m)
	optimize!(m)

	for i=c_index cHist[i,t]=value(c[i,t]) end
	for i=s_index sHist[i,t]=value(s[i,t]) end
	for i=nodes thetaHist[i,t]=value(theta[i,t]) end
	for l=edges
		fpHist[l,t]=value(fp[l,t])
		fpToHist[l,t]=value(fpTo[l,t])
	end
	for l=dist_edges
		fqHist[l,t]=value(fq[l,t])
		fqToHist[l,t]=value(fqTo[l,t])
	end
	obj_times[t]=(+mapreduce(sum,+,-dualsol.priceP[i,t]*(mapreduce(sum,+,fpToHist[e,t] for e=n_to[i];init=0)
		+ mapreduce(sum,+,fpHist[e,t] for e=n_fr[i];init=0)) for i=nodes;init=0)
	+mapreduce(sum,+,-dualsol.priceQ[i,t]*(mapreduce(sum,+,fqToHist[e,t] for e=dist_n_to[i];init=0)
		+ mapreduce(sum,+,fqHist[e,t] for e=dist_n_fr[i];init=0)) for i=dist_nodes;init=0)
	-mapreduce(sum,+,-dualsol.priceP[i,t]*(mapreduce(sum,+,sol.fpTo[e,t] for e=n_to[i];init=0)
		+ mapreduce(sum,+,sol.fp[e,t] for e=n_fr[i];init=0)) for i=nodes;init=0)
	-mapreduce(sum,+,-dualsol.priceQ[i,t]*(mapreduce(sum,+,sol.fqTo[e,t] for e=dist_n_to[i];init=0)
		+ mapreduce(sum,+,sol.fq[e,t] for e=dist_n_fr[i];init=0)) for i=dist_nodes;init=0))
	end

	return (obj_times,
	PFsol(false,sol.pg,sol.qg,cHist,sHist,thetaHist,fpHist,fqHist,fpToHist,
	fqToHist,sol.slackPUp,sol.slackPDown,sol.slackQUp,sol.slackQDown))
end

function dist_line_subproblem(cur_e::Int64,sol::PFsol,dualsol::Dualsol)::Bool
	m=Model(optimizer_with_attributes(Ipopt.Optimizer,"print_level"=>0))
	### VARIABLES ###
	@variable(m, 0 <= c[[(e_fr[cur_e],e_to[cur_e]),(e_fr[cur_e],e_fr[cur_e]),(e_to[cur_e],e_to[cur_e])],times])
	@variable(m, s[[(e_fr[cur_e],e_to[cur_e])],times])
	@variable(m, -pi <= theta[[e_fr[cur_e],e_to[cur_e]],times] <= pi)

	### FLOW DEFINITION ###
	@variable(m, fp[cur_e,times])
	@variable(m, fpTo[cur_e,times])
	DefFpDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	DefFpToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	for t=times
		if R_pu[cur_e] == 0.0
			get!(DefFpDist,(cur_e,t), @constraint(m, fp[cur_e,t] + fpTo[cur_e,t] == 0))
			@constraint(m, c[(e_fr[cur_e],e_fr[cur_e]),t] ==  c[(e_to[cur_e],e_to[cur_e]),t])
			@constraint(m, c[(e_fr[cur_e],e_fr[cur_e]),t] ==  c[(e_fr[cur_e],e_to[cur_e]),t])
			@constraint(m, s[(e_fr[cur_e],e_to[cur_e]),t] == 0.0)
			@constraint(m, theta[e_fr[cur_e]] == theta[e_to[cur_e]])
		else
			get!(DefFpDist,(cur_e,t),
			@constraint(m,
				fp[cur_e,t] == -Gbus[e_fr[cur_e],e_to[cur_e]]*c[(e_fr[cur_e],e_fr[cur_e]),t]
				+ Gbus[e_fr[cur_e],e_to[cur_e]]*c[(e_fr[cur_e],e_to[cur_e]),t]
				- Bbus[e_fr[cur_e],e_to[cur_e]]*s[(e_fr[cur_e],e_to[cur_e]),t])
			)
			get!(DefFpToDist,(cur_e,t),
			@constraint(m,
				fpTo[cur_e,t] == -Gbus[e_to[cur_e],e_fr[cur_e]]*c[(e_to[cur_e],e_to[cur_e]),t]
				+ Gbus[e_to[cur_e],e_fr[cur_e]]*c[(e_fr[cur_e],e_to[cur_e]),t]
				+ Bbus[e_to[cur_e],e_fr[cur_e]]*s[(e_fr[cur_e],e_to[cur_e]),t]
			))
		end
	end
	@variable(m, fq[cur_e,times])
	@variable(m, fqTo[cur_e,times])
	DefFqDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	DefFqToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	for t=times
		if X_pu[cur_e] == 0.0
			get!(DefFqDist,(cur_e,t), @constraint(m, fq[cur_e,t] + fqTo[cur_e,t] == 0))
		else
			get!(DefFqDist,(cur_e,t),
			@constraint(m,
				fq[cur_e,t] == Bbus[e_fr[cur_e],e_to[cur_e]]*c[(e_fr[cur_e],e_fr[cur_e]),t]
				- Bbus[e_fr[cur_e],e_to[cur_e]]*c[(e_fr[cur_e],e_to[cur_e]),t]
				- Gbus[e_fr[cur_e],e_to[cur_e]]*s[(e_fr[cur_e],e_to[cur_e]),t])
			)
			get!(DefFqToDist,(cur_e,t),
			@constraint(m,
				fqTo[cur_e,t] == Bbus[e_to[cur_e],e_fr[cur_e]]*c[(e_to[cur_e],e_to[cur_e]),t]
				- Bbus[e_to[cur_e],e_fr[cur_e]]*c[(e_fr[cur_e],e_to[cur_e]),t]
				+ Gbus[e_to[cur_e],e_fr[cur_e]]*s[(e_fr[cur_e],e_to[cur_e]),t]
			))
		end
	end

	### CONSTRAINTS ###
	@constraint(m, VoltageLimitsDown[i=[e_fr[cur_e],e_to[cur_e]],t=times],
	-c[(i,i),t] <= -vsq_min[i])

	@constraint(m, VoltageLimitsUp[i=[e_fr[cur_e],e_to[cur_e]],t=times],
	c[(i,i),t] <= vsq_max[i])

	@constraint(m, LinePowerLimit1[e=intersect(line_limit,cur_e),t=times],
	fp[e,t]^2 + fq[e,t]^2 <= rate_a[e]^2)

	@constraint(m, LinePowerLimit2[e=intersect(line_limit,cur_e),t=times],
	fpTo[e,t]^2 + fqTo[e,t]^2 <= rate_a[e]^2)

	mod_s_ind=[(e_fr[e],e_to[e]) for e=cur_e if R_pu[e]!=0.0]
	@NLconstraint(m, Quadconstraint[e=mod_s_ind,t=times],
	c[e,t]^2 + s[e,t]^2 == c[(e[1],e[1]),t]*c[(e[2],e[2]),t])

	@NLconstraint(m, Tanconstraint[e=mod_s_ind,t=times],
	c[e,t]*sin(theta[e[1],t]-theta[e[2],t])
	+ s[e,t]*cos(theta[e[1],t]-theta[e[2],t]) == 0)

	@expression(m, obj,
	+sum(- dualsol.priceP[e_fr[cur_e],t]*(fp[cur_e,t]-sol.fp[cur_e,t])
	- dualsol.priceP[e_to[cur_e],t]*(fpTo[cur_e,t]-sol.fpTo[cur_e,t]) for t=times)
	+sum(- dualsol.priceQ[e_fr[cur_e],t]*(fq[cur_e,t]-sol.fq[cur_e,t])
	- dualsol.priceQ[e_to[cur_e],t]*(fqTo[cur_e,t]-sol.fqTo[cur_e,t]) for t=times)
	)

	#@constraint(m, ObjPos, obj >= 0)
	#ext_edges=[e for e=tr_edges if ((e_fr[e] in tr_nodes && e_to[e] in dist_nodes) || (e_to[e] in tr_nodes && e_fr[e] in dist_nodes))]
	#for e=ext_edges t=times @constraint(fp[e,t]==sol.fp[e,t]) end end
	### OBJECTIVE ###
	@objective(m, Max, obj)

	#println(m)
	optimize!(m)
	return (objective_value(m)<1e-5)
end

function preprocess_loc_op(sol::PFsol,dualsol::Dualsol)::Set{Int64}
	preprocess_lines=Set{Int64}()
	for e=dist_edges
		if dist_line_subproblem(e,sol,dualsol)
			push!(preprocess_lines,e)
		end
	end
	return preprocess_lines
end

function max_diff_sol_inj(sol1::PFsol,sol2::PFsol)
	diff_inj=Dict(
		(i,t) => abs(mapreduce(sum,+,init=0,sol1.pg[sb]-sol2.pg[sb] for sb=seg_bid if sb[1]==i && sb[5]==t))
	for i=nodes for t=times)
	for i=nodes for t=times if diff_inj[i,t] > 1e-6
		diff_inj[i,t]/=sum(p_max[sb]-p_min[sb] for sb=seg_bid if sb[1]==i && sb[5]==t)
	end end end
	return norm(values(diff_inj))/length(diff_inj)
end
