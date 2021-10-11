
function test_feasibility_subnet(dn::String)::Tuple{PFsol,Dualsol,Float64}
	dn_nodes=DNs[dn]
	dn_edges=edges_DN[dn]
	dn_c_index=c_index_DN[dn]
	dn_s_index=s_index_DN[dn]
	dn_seg_bid=seg_bid_DN[dn]
	dn_dist_seg_bid=dist_seg_bid_DN[dn]
	dn_bid=bid_DN[dn]
	dn_ext_nodes=ext_nodes_DN[dn]
	dn_ext_edges=ext_edges_DN[dn]
	dn_RampIndexMinIncr=RampIndexMinIncr_DN[dn]
	dn_RampIndexMaxIncr=RampIndexMaxIncr_DN[dn]
	m=Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
	### VARIABLES ###
	@variable(m, 0 <= c[dn_c_index,times])
	@variable(m, s[dn_s_index,times])
	@variable(m, pg[dn_seg_bid]==0)
	@expression(m, net_pg_bid[bd=dn_bid], sum(pg[sb] for sb=bid_2_seg_bids[bd]))
	@variable(m, qg[dn_seg_bid])
	@expression(m, net_qg_bid[bd=dn_bid], sum(qg[sb] for sb=bid_2_seg_bids[bd]))
	@variable(m, -pi <= theta[union(dn_nodes,dn_ext_nodes),times] <= pi)

	### FLOW DEFINITION ###
	@variable(m, fp[union(dn_edges,dn_ext_edges),times])
	@variable(m, fpTo[union(dn_edges,dn_ext_edges),times])

	DefFpDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	DefFpToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	for l=dn_edges for t=times
	if R_pu[l] == 0.0
		get!(DefFpDist,(l,t), @constraint(m, fp[l,t] + fpTo[l,t] == 0))
		@constraint(m, c[(e_fr[l],e_to[l]),t] ==  c[(e_to[l],e_to[l]),t])
		@constraint(m, c[(e_fr[l],e_fr[l]),t] ==  c[(e_fr[l],e_to[l]),t])
		#fix(s[(e_fr[l],e_to[l]),t],0.0;force=true)
		@constraint(m, s[(e_fr[l],e_to[l]),t]==0.0)
		@constraint(m, theta[e_fr[l],t]==theta[e_to[l],t])
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
			== 0)
		)
	end
	end end

	@variable(m, fq[dn_edges,times])
	@variable(m, fqTo[dn_edges,times])

	DefFqDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	DefFqToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	for l=dn_edges for t=times
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
			== 0)
		)
	end
	end end

	### CONSTRAINTS ###
	@variable(m, 0 <= relax_sa[dn_seg_bid] <= 1)
	@constraint(m, PDown[sb=dn_seg_bid], -pg[sb] <= -relax_sa[sb]*p_min[sb])
	@constraint(m, PUp[sb=dn_seg_bid], pg[sb] <= relax_sa[sb]*p_max[sb])

	@constraint(m, PowerBalanceP[i=dn_nodes,t=times],
	sum(sum(pg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
	+ net_injP[i,t]
	#+ slackPUp[i,t]
	==
	+ sum(fpTo[e,t] for e=n_to[i])
	+ sum(fp[e,t] for e=n_fr[i])
	#+ slackPDown[i,t]
	)

	@constraint(m, PowerBalanceQ[i=dn_nodes,t=times],
	sum(sum(qg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
	+ net_injQ[i,t]
	#+ slackQUp[i,t]
	==
	+ sum(fqTo[e,t] for e=dist_n_to[i])
	+ sum(fq[e,t] for e=dist_n_fr[i])
	#+ slackQDown[i,t]
	)

	@constraint(m, VoltageLimitsDown[i=intersect(dn_nodes,v_limit),t=times],
	-c[(i,i),t] <= -vsq_min[i])

	@constraint(m, VoltageLimitsUp[i=intersect(dn_nodes,v_limit),t=times],
	c[(i,i),t] <= vsq_max[i])

	@constraint(m, LinePowerLimit1[e=intersect(line_limit,dn_edges),t=times],
	fp[e,t]^2 + fq[e,t]^2 <= rate_a[e]^2
	)

	@constraint(m, LinePowerLimit2[e=intersect(line_limit,dn_edges),t=times],
	fpTo[e,t]^2 + fqTo[e,t]^2 <= rate_a[e]^2)

	mod_s_ind=[(e_fr[e],e_to[e]) for e=dn_edges if R_pu[e]!=0.0]
	@NLconstraint(m, Quadconstraint[e=mod_s_ind,t=times],
	c[e,t]^2 + s[e,t]^2 == c[(e[1],e[1]),t]*c[(e[2],e[2]),t])

	@NLconstraint(m, Tanconstraint[e=mod_s_ind,t=times],
	c[e,t]*sin(theta[e[1],t]-theta[e[2],t])
	+ s[e,t]*cos(theta[e[1],t]-theta[e[2],t]) == 0)


	@constraint(m, RampMinIncr[rcI=dn_RampIndexMinIncr],
	+ sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
	- sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
	- get(RCRealPowerIncr,(rcI[1],rcI[2],rcI[4]),0)
	>= 0)

	@constraint(m, RampMaxIncr[rcI=dn_RampIndexMaxIncr],
	+ sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
	+ get(RCRealPowerIncr,(rcI[1],rcI[2],rcI[4]),0)
	- sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
	>= 0)

	@constraint(m, DiscQtoP[BQpd=[qpd for qpd=QPDiscConstrIndex if qpd[1] in dn_nodes]],
	get(QPDiscMax,(BQpd[2],BQpd[3],BQpd[5]),0)^2
	- net_pg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]^2
	- net_qg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]^2 >= 0
	)

	@constraint(m, ConstrainToUpFromLineQtoP[BHpc=[qpd for qpd=QPHPConstrIndex1 if qpd[1] in dn_nodes]],
	net_qg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	-QPHPOffset[(BHpc[5],BHpc[2])]
	-QPHPSlope[(BHpc[5],BHpc[2])]*net_pg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	>= 0)

	@constraint(m, ConstrainToDownFromLineQtoP[BHpc=[qpd for qpd=QPHPConstrIndex0 if qpd[1] in dn_nodes]],
	QPHPOffset[(BHpc[5],BHpc[2])]
	+QPHPSlope[(BHpc[5],BHpc[2])]*net_pg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	-net_qg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])] >= 0)

	### Flow definition for interconnection ###
	@constraint(m, DefFpTr[e=dn_ext_edges,t=times],
	X_pu[e]*fp[e,t] == (theta[e_fr[e],t]-theta[e_to[e],t]))
	@constraint(m, FlowLimUp[e=dn_ext_edges,t=times],
	fp[e,t] <= rate_a[e])
	@constraint(m, FlowLimDown[e=dn_ext_edges,t=times],
	-fp[e,t] <= rate_a[e])
	@constraint(m, DefFpX0[e=dn_ext_edges, t=times],
	fp[e,t] + fpTo[e,t] == 0)

	# ### OBJECTIVE ###
	if !isempty(dn_seg_bid)
		@objective(m, Min, sum(cost[sb][1]*pg[sb] for sb=dn_seg_bid))
	end

	optimize!(m)
	#@assert primal_status(m)==MOI.FEASIBLE_POINT
	print(primal_status(m),"  ")
	## only build a partial solution for this DN ##
	pgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(pg[j]) for j=dn_seg_bid)
	qgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(qg[j]) for j=dn_dist_seg_bid)
	cHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((i,t) => value(c[i,t]) for i=dn_c_index for t=times)
	sHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((l,t) => value(s[l,t]) for l=dn_s_index for t=times)
	thetaHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(theta[i,t]) for i=dn_nodes for t=times)
	fpHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fp[l,t]) for l=dn_edges for t=times)
	fqHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fq[l,t]) for l=dn_edges for t=times)
	fpToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fpTo[l,t]) for l=dn_edges for t=times)
	fqToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fqTo[l,t]) for l=dn_edges for t=times)
	priceP=Dict{Tuple{Int64,Int64},Float64}((i,t) => dual(PowerBalanceP[i,t])
	for i=dn_nodes for t=times)
	priceQ=Dict{Tuple{Int64,Int64},Float64}((i,t) => dual(PowerBalanceQ[i,t]) for i=intersect(dn_nodes,dist_nodes) for t=times)
	slackPUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => 0.0
	for i=dn_nodes for t=times)
	slackPDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => 0.0
	for i=dn_nodes for t=times)
	slackQUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => 0.0
	for i=intersect(dn_nodes,dist_nodes) for t=times)
	slackQDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => 0.0
	for i=intersect(dn_nodes,dist_nodes) for t=times)
	pfsol=PFsol(true,pgHist,qgHist,cHist,sHist,thetaHist,
	fpHist,fqHist,fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
	dual_sol=Dualsol(priceP,priceQ)
	return (pfsol,dual_sol,solve_time(m))
end

function test_feasibility_TN(solver::Symbol=:Mosek)::Tuple{PFsol,Dualsol,Float64}
	m=Model(optimizer_with_attributes(Mosek.Optimizer,
	"MSK_IPAR_PRESOLVE_USE" => 0,
	"MSK_IPAR_LOG" => 0,
	"MSK_IPAR_INTPNT_BASIS" => 0
	# MSK_IPAR_PRESOLVE_MAX_NUM_REDUCTIONS=0,
	# MSK_IPAR_LOG_PRESOLVE=5,
	# MSK_DPAR_MIO_TOL_FEAS=1e-8,
	# MSK_DPAR_MIO_TOL_ABS_RELAX_INT=1e-8,
	# MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-8,
	# MSK_DPAR_INTPNT_CO_TOL_PFEAS=1e-8,
	# MSK_DPAR_ANA_SOL_INFEAS_TOL=1e-8,
	# MSK_DPAR_MIO_TOL_REL_GAP=1e-2,
	))
	if solver==:Gurobi
		m=Model(optimizer_with_attributes(Gurobi.Optimizer,
		#"NumericFocus" => 2,
		"BarHomogeneous" => 1,
		"Presolve" => -1,
		"QCPDual" => 1,
		#Aggregate=0
		# # BarQCPConvTol=1e-6,
		# # FeasibilityTol=1e-6,
		# # IntFeasTol=1e-6,
		#"PreQLinearize" => 0,
		# "MIQCPMethod" => 0
		# #PreMIQCPForm=0
		))
	end
	### CONTINUOUS VARIABLES ###
	@variable(m, pg[sb=tr_seg_bid])
	@expression(m, net_pg_bid[bd=tr_bid], sum(pg[sb] for sb=bid_2_seg_bids[bd]))
	@variable(m, -pi <= theta_tr[union(tr_nodes_ext,slack_bus),times] <= pi)
	@variable(m, 0 <= relax_sa[tr_seg_bid] <= 1)

	### FLOW DEFINITION ###
	@variable(m, fp[tr_edges,times])
	@constraint(m, DefFpTr[e=tr_edges,t=times],
	X_pu[e]*fp[e,t] == (theta_tr[e_fr[e],t]-theta_tr[e_to[e],t]))

	@variable(m,fpTo[tr_edges,times])
	@constraint(m, DefFpX0[e=tr_edges,t=times],
	fp[e,t] + fpTo[e,t] == 0)

	### CONSTRAINTS ###
	@constraint(m, PDownSA[sb=tr_seg_bid], -pg[sb] <= -relax_sa[sb]*p_min[sb])
	@constraint(m, PUpSA[sb=tr_seg_bid], pg[sb] <= relax_sa[sb]*p_max[sb])

	@constraint(m, TrPowerFlowLimitUp[e=intersect(line_limit,tr_edges),t=times],
	fp[e,t] <= rate_a[e])

	@constraint(m, TrPowerFlowLimitDown[e=intersect(line_limit,tr_edges),t=times],
	-fp[e,t] <= rate_a[e])

	for i=slack_bus for t=times @constraint(m, theta_tr[i,t]==0) end end

	@constraint(m, PowerBalanceP[i=tr_nodes,t=times],
	sum(sum(pg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
	+ net_injP[i,t]
	==
	+ sum(fpTo[e,t] for e=n_to[i])
	+ sum(fp[e,t] for e=n_fr[i])
	)

	@constraint(m, RampMinIncr[rcI=tr_RampIndexMinIncr],
	+ sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
	- sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
	- get(RCRealPowerIncr,(rcI[1],rcI[2],rcI[4]),0)
	>= 0)

	@constraint(m, RampMaxIncr[rcI=tr_RampIndexMaxIncr],
	+ sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
	+ get(RCRealPowerIncr,(rcI[1],rcI[2],rcI[4]),0)
	- sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
	>= 0)

	### OBJECTIVE ###
	@objective(m, Min,
	+ sum(cost[sb][1]*pg[sb] for sb=tr_seg_bid)
	)

	optimize!(m)
	# @assert termination_status(m)==MOI.OPTIMAL
	# @assert primal_status(m)==MOI.FEASIBLE_POINT
	println(termination_status(m))
	println(primal_status(m))

	pgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(pg[j]) for j=tr_seg_bid)
	qgHist = 0
	cHist = 0
	sHist = 0
	theta_trHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(theta_tr[i,t]) for i=union(slack_bus,tr_nodes_ext) for t=times)
	fpHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fp[l,t]) for l=tr_edges for t=times)
	fqHist = 0
	fpToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fpTo[l,t]) for l=tr_edges for t=times)
	fqToHist = 0
	slackPUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => 0.0 for i=tr_nodes for t=times)
	slackPDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => 0.0 for i=tr_nodes for t=times)
	slackQUpHist = 0
	slackQDownHist = 0
	dc_sol=PFsol(false,pgHist,qgHist,cHist,sHist,theta_trHist,fpHist,fqHist,
	fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)

	dual_sol=Dualsol(0,0)
	if has_duals(m)
		priceP=Dict{Tuple{Int64,Int64},Float64}((i,t) => dual(PowerBalanceP[i,t])
		for i=tr_nodes for t=times)
		priceQ=0
		dual_sol=Dualsol(priceP,priceQ)
	else
		println(">> Dual solution not computed.")
	end
	return (dc_sol,dual_sol,solve_time(m))
end

function test_feasibility()::Tuple{PFsol,Dualsol,Dict{String,Float64}}
	solve_times=Dict{String,Float64}()
	println("Starting decomposition...")
	(sol_TN,dualsol_TN,solve_times["TN"])=test_feasibility_TN()
	println("> SOLVED Transmission network problem.")
	sol_DN=Dict{String,PFsol}()
	dualsol_DN=Dict{String,Dualsol}()
	println("Starting $(length(root_dist_node)) Distribution networks problems.")
	for (num,dn) in enumerate(keys(DNs))
		(sol_DN[dn],dualsol_DN[dn],solve_times[dn])=test_feasibility_subnet(dn)
		println(">$num. SOLVED Distribution network problem: $dn.")
	end
	pgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> sol_TN.pg[j] for j=tr_seg_bid)
	for dn=keys(DNs) for sb=seg_bid_DN[dn]
			pgHist[sb]=sol_DN[dn].pg[sb]
	end end
	qgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}()
	for dn=keys(DNs) for j=seg_bid_DN[dn] qgHist[j] = sol_DN[dn].qg[j] end end
	cHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}()
	for dn=keys(DNs) for i=c_index_DN[dn] for t=times cHist[i,t] = sol_DN[dn].c[i,t] end end end
	sHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}()
	for dn=keys(DNs) for i=s_index_DN[dn] for t=times sHist[i,t] = sol_DN[dn].s[i,t] end end end
	thetaHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => sol_TN.theta[i,t] for i=tr_nodes for t=times)
	for dn=keys(DNs) for i=DNs[dn] for t=times
			thetaHist[i,t]=sol_DN[dn].theta[i,t]
	end end end
	fpHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => sol_TN.fp[l,t] for l=tr_edges for t=times)
	fpToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => sol_TN.fpTo[l,t] for l=tr_edges for t=times)
	for dn=keys(DNs) for e=edges_DN[dn] for t=times
		fpHist[e,t]=sol_DN[dn].fp[e,t]
		fpToHist[e,t]=sol_DN[dn].fpTo[e,t]
	end end end
	fqHist = Dict{Tuple{Int64,Int64},Float64}()
	fqToHist = Dict{Tuple{Int64,Int64},Float64}()
	for dn=keys(DNs) for e=edges_DN[dn] for t=times
		fqHist[e,t]=sol_DN[dn].fq[e,t]
		fqToHist[e,t]=sol_DN[dn].fqTo[e,t]
	end end end
	pricePHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => dualsol_TN.priceP[i,t] for i=tr_nodes for t=times)
	for dn=keys(DNs) for i=DNs[dn] for t=times
		pricePHist[i,t]=dualsol_DN[dn].priceP[i,t]
	end end end
	priceQHist= Dict{Tuple{Int64,Int64},Float64}()
	for dn=keys(DNs) for i=DNs[dn] for t=times priceQHist[i,t]=dualsol_DN[dn].priceQ[i,t] end end end
	slackPUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => sol_TN.slackPUp[i,t] for i=tr_nodes for t=times)
	slackPDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => sol_TN.slackPDown[i,t] for i=tr_nodes for t=times)
	for dn=keys(DNs) for i=DNs[dn] for t=times
		slackPUpHist[i,t]=sol_DN[dn].slackPUp[i,t]
		slackPDownHist[i,t]=sol_DN[dn].slackPDown[i,t]
	end end end
	slackQUpHist = Dict{Tuple{Int64,Int64},Float64}()
	slackQDownHist = Dict{Tuple{Int64,Int64},Float64}()
	for dn=keys(DNs) for i=DNs[dn] for t=times
		slackQUpHist[i,t] = sol_DN[dn].slackQUp[i,t]
		slackQDownHist[i,t] = sol_DN[dn].slackQDown[i,t]
	end end end

	println("\n\n >>> Informations about the solution <<<")
	MaxSlack=maximum(abs,union(values(slackPUpHist),values(slackQUpHist),values(slackPDownHist), values(slackQDownHist)))
	println("Max Slack: $MaxSlack")
	pfsol=PFsol(MaxSlack>1e-5,pgHist,qgHist,cHist,sHist,thetaHist,
	fpHist,fqHist,fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
	dual_sol=Dualsol(pricePHist,priceQHist)
	ComputeGaps(pfsol)
	println("SW: $(SW(pfsol))")
	println("\nTotal solve time: $(sum(values(solve_times)))")
	println("Max solve time: $(findmax(solve_times))")
	return (pfsol,dual_sol,solve_times)
end
