__precompile__()
push!(LOAD_PATH,".")

using Gurobi, CSV, JuMP, DataFrames, Ipopt, Mosek, MosekTools,LinearAlgebra

ChargeData=true

mutable struct PFsol
	isPFsol
	pg
	qg
	c
	s
	theta
	fp
	fq
	fpTo
	fqTo
	slackPUp
	slackPDown
	slackQUp
	slackQDown
end

mutable struct Binsol
	sa
	qa
	qta
	alpha
	omega
end

mutable struct Dualsol
	priceP
	priceQ
end

#if ChargeData
print("Charging .csv files...")
#TestCase="../data/Small/7nodes_ex/"
#TestCase="../data/Small/3nodes_ex_no_v_constr/"
#TestCase="../data/Small/small_T_and_D/"
### Test case for which we run the One shot model ###
TestCase="../data/Italy/CSV_Italy_693_T00_original/"
#TestCase="../data/Italy/CSV_Italy_652_T00_simplified/"
#TestCase="../data/Italy/CSV_Italy_652_T66_simplified/"
#TestCase="../data/Denmark/CSV_Denmark_401_T88/"
#TestCase="../data/Small/3nodes_ex_thesis/"



include("data_charging.jl")
println("DONE")

include("feasibility_testing.jl")
include("loc_computation.jl")
include("build_set_tables.jl")
include("res_supp_function.jl")
include("res_supp_function_bin.jl")

function solveMIDCSOCP(exp_flow=0, solver::Symbol=:Mosek)::Tuple{PFsol,Binsol}
	m=Model(optimizer_with_attributes(Mosek.Optimizer,
	"MSK_IPAR_PRESOLVE_USE" => 1,
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
	"MSK_DPAR_OPTIMIZER_MAX_TIME" => 3600
	))
	if solver==:Gurobi
		m=Model(optimizer_with_attributes(Gurobi.Optimizer,
		"NumericFocus" => 2,
		"BarHomogeneous" => 1,
		"Presolve" => -1,
		"QCPDual" => 0,
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
	@variable(m, 0 <= c[e=c_index,times])
	@variable(m, s[s_index,times])
	@variable(m, pg[sb=seg_bid])
	@expression(m, net_pg_bid[bd=bid], sum(pg[sb] for sb=bid_2_seg_bids[bd]))
	#@variable(m, net_pg_bid[bd=bid])
	#@constraint(m, [bd=bid], net_pg_bid[bd] == sum(pg[sb] for sb=bid_2_seg_bids[bd]))
	@variable(m, qg[sb=dist_seg_bid])
	@expression(m, net_qg_bid[bd=dist_bid], sum(qg[sb] for sb=bid_2_seg_bids[bd]))
	#@variable(m, net_qg_bid[bd=dist_bid])
	#@constraint(m, [bd=dist_bid], net_qg_bid[bd] == sum(qg[sb] for sb=bid_2_seg_bids[bd]))
	@variable(m, -pi <= theta_tr[union(tr_nodes_ext,slack_bus),times] <= pi)

	### BINARY VARIABLES ###
	@variable(m, sa[seg_bid], Bin)
	@variable(m, qa[bid], Bin)
	@variable(m, qta[qt_bid], Bin)
	@variable(m, alpha[AlphaOmegaIndex], Bin)
	@variable(m, omega[AlphaOmegaIndex], Bin)

	### SLACK VARIABLES ###
	@variable(m, slackPUp[nodes,times] == 0)
	@variable(m, slackPDown[nodes,times] == 0)
	@variable(m, slackQUp[dist_nodes,times] == 0)
	@variable(m, slackQDown[dist_nodes,times] == 0)

	### FLOW DEFINITION ###
	@variable(m, fp[edges,times])
	@variable(m, fpTo[edges,times])
	@constraint(m, DefFpTr[e=tr_edges,t=times],
	X_pu[e]*fp[e,t]== (theta_tr[e_fr[e],t]-theta_tr[e_to[e],t]))
	@constraint(m, DefFpX0[e=tr_edges,t=times],
	fp[e,t] + fpTo[e,t] == 0)

	DefFpDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	DefFpToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	for l=dist_edges for t=times
		if R_pu[l] == 0.0
			get!(DefFpDist,(l,t), @constraint(m, fp[l,t] + fpTo[l,t] == 0))
			@constraint(m, c[(e_fr[l],e_to[l]),t] ==  c[(e_to[l],e_to[l]),t])
			@constraint(m, c[(e_fr[l],e_fr[l]),t] ==  c[(e_fr[l],e_to[l]),t])
			@constraint(m, s[(e_fr[l],e_to[l]),t] == 0)
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
	@constraint(m, PDownSA[sb=seg_bid], -pg[sb] <= -sa[sb]*p_min[sb])
	@constraint(m, PUpSA[sb=seg_bid], pg[sb] <= sa[sb]*p_max[sb])

	@constraint(m, TrPowerFlowLimitUp[e=intersect(line_limit,tr_edges),t=times],
	fp[e,t] <= rate_a[e])

	@constraint(m, TrPowerFlowLimitDown[e=intersect(line_limit,tr_edges),t=times],
	-fp[e,t] <= rate_a[e])

	for i=slack_bus for t=times @constraint(m, theta_tr[i,t] == 0) end end

	@constraint(m, PowerBalanceP[i=nodes,t=times],
	sum(sum(pg[sb] for sb=bid_2_seg_bids[B])
		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
	+ net_injP[i,t]
	+ slackPUp[i,t]
	==
	+ sum(fpTo[e,t] for e=n_to[i])
	+ sum(fp[e,t] for e=n_fr[i])
	+ slackPDown[i,t]
	)

	# ### tentative constraint ###
	# @constraint(m, TotalPBalance[t=times],
	# + sum(pg[sb] for sb=seg_bid if sb[5]==t)
	# + sum(net_injP[i,t] for i=nodes)
	# + sum(slackPUp[i,t] for i=nodes)
	# ==
	# + sum(fpTo[e,t] for e=edges)
	# + sum(fp[e,t] for e=edges)
	# + sum(slackPDown[i,t] for i=nodes)
	# )
	# ###

	@constraint(m, PowerBalanceQ[i=dist_nodes,t=times],
	sum(sum(qg[sb] for sb=bid_2_seg_bids[B])
		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
	+ net_injQ[i,t]
	+ slackQUp[i,t]
	==
	+ sum(fqTo[e,t] for e=dist_n_to[i])
	+ sum(fq[e,t] for e=dist_n_fr[i])
	+ slackQDown[i,t]
	)

	# ### tentative constraint ###
	# @constraint(m, TotalQBalance[t=times],
	# + sum(qg[sb] for sb=dist_seg_bid if sb[5]==t)
	# + sum(net_injQ[i,t] for i=dist_nodes)
	# + sum(slackQUp[i,t] for i=dist_nodes)
	# ==
	# + sum(fqTo[e,t] for e=dist_edges)
	# + sum(fq[e,t] for e=dist_edges)
	# + sum(slackQDown[i,t] for i=dist_nodes)
	# )
	# ###

	@constraint(m, VoltageLimitsDown[i=v_limit,t=times],
	-c[(i,i),t] <= -vsq_min[i])

	@constraint(m, VoltageLimitsUp[i=v_limit,t=times],
	c[(i,i),t] <= vsq_max[i])

	mod_s_ind=[(e_fr[e],e_to[e]) for e=dist_edges if R_pu[e]!=0.0]
	if solver!=:Gurobi
		@constraint(m, LinePowerLimit1[e=intersect(line_limit,dist_edges),t=times],
		[rate_a[e],fp[e,t],fq[e,t]] in SecondOrderCone())

		@constraint(m, LinePowerLimit2[e=intersect(line_limit,dist_edges),t=times],
		[rate_a[e],fpTo[e,t],fqTo[e,t]] in SecondOrderCone())

		@constraint(m, SOCPconstraint[e=mod_s_ind,t=times],
		[0.5*c[(e[1],e[1]),t],c[(e[2],e[2]),t],c[e,t],s[e,t]] in RotatedSecondOrderCone())

		@constraint(m, DiscQtoP[BQpd=QPDiscConstrIndex],
		[QPDiscMax[(BQpd[2],BQpd[3],BQpd[5])],
		net_pg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])],
		net_qg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]]
		in SecondOrderCone())
	else
		@constraint(m, LinePowerLimit1[e=intersect(line_limit,dist_edges),t=times],
		fp[e,t]^2+fq[e,t]^2 <= rate_a[e]^2)

		@constraint(m, LinePowerLimit2[e=intersect(line_limit,dist_edges),t=times],
		fpTo[e,t]^2+fqTo[e,t]^2 <= rate_a[e]^2)

		@constraint(m, SOCPconstraint[e=mod_s_ind,t=times],
		c[e,t]^2+s[e,t]^2 <= c[(e[1],e[1]),t]*c[(e[2],e[2]),t])

		@constraint(m, DiscQtoP[BQpd=QPDiscConstrIndex],
		QPDiscMax[(BQpd[2],BQpd[3],BQpd[5])]^2
		- net_pg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]^2
		- net_qg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]^2 >= 0
		)
	end

	@constraint(m, RampMinIncr[rcI=RampIndexMinIncr],
	+ sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
	- sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
	- RCRealPowerIncr[(rcI[1],rcI[2],rcI[4])]
	>= 0)

	@constraint(m, RampMaxIncr[rcI=RampIndexMaxIncr],
	+ sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
	+ RCRealPowerIncr[(rcI[1],rcI[2],rcI[4])]
	- sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
	>= 0)

	@constraint(m, ConstrainToUpFromLineQtoP[BHpc=QPHPConstrIndex1],
	net_qg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	-QPHPOffset[(BHpc[5],BHpc[2])]
	-QPHPSlope[(BHpc[5],BHpc[2])]*net_pg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	>= 0)

	@constraint(m, ConstrainToDownFromLineQtoP[BHpc=QPHPConstrIndex0],
	QPHPOffset[(BHpc[5],BHpc[2])]
	+QPHPSlope[(BHpc[5],BHpc[2])]*net_pg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	-net_qg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])] >= 0)

	### BINARY CONSTRAINTS ###
	@constraint(m, ExQtBidsConstr[i=ex_qt_bid_id],
	sum(qta[ExQt] for ExQt=ex_qt_bid[i]) <= 1)

	@constraint(m, QtaDefLow[qt=qt_bid],
	qta[qt] <= sum(qa[bd] for bd=bid if bd[2]==qt))

	@constraint(m, QtaDefHigh[bd=bid],
	qa[bd] <= qta[bd[2]])

	@constraint(m, DefQa[sb=seg_bid],
	sa[sb] <= qa[(sb[1],sb[2],sb[3],sb[5])])

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

	### Fix interconnection flow if provided ###
	if exp_flow!=0
		@constraint(m, FixInterFlow[e=ext_edges,t=times], fp[e,t]==exp_flow[e,t])
	end
	### OBJECTIVE ###
	@objective(m, Min,
	+ sum(cost[sb][1]*pg[sb] for sb=seg_bid)
	+ cost_ls*(sum(slackPUp[i,t] + slackPDown[i,t] for i=nodes for t=times)
		 + sum(slackQUp[i,t] + slackQDown[i,t] for i=dist_nodes for t=times))
	)

	optimize!(m)
	println(termination_status(m))
	println(primal_status(m))

	pgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(pg[j]) for j=seg_bid)
	qgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(qg[j]) for j=dist_seg_bid)
	cHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((i,t) => value(c[i,t]) for i=c_index for t=times)
	sHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((l,t) => value(s[l,t]) for l=s_index for t=times)
	theta_trHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(theta_tr[i,t]) for i=union(slack_bus,tr_nodes_ext) for t=times)
	merge!(theta_trHist,Dict{Tuple{Int64,Int64},Float64}((i,t) => 0.0 for i=setdiff(nodes,tr_nodes_ext) for t=times))
	fpHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fp[l,t]) for l=edges for t=times)
	fqHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fq[l,t]) for l=dist_edges for t=times)
	fpToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fpTo[l,t]) for l=edges for t=times)
	fqToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fqTo[l,t]) for l=dist_edges for t=times)
	slackPUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPUp[i,t]) for i=nodes for t=times)
	slackPDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPDown[i,t]) for i=nodes for t=times)
	slackQUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackQUp[i,t]) for i=dist_nodes for t=times)
	slackQDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackQDown[i,t]) for i=dist_nodes for t=times)
	socp_sol=PFsol(false,pgHist,qgHist,cHist,sHist,theta_trHist,fpHist,fqHist,
	fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)

	socp_sol=heurTheta2(socp_sol)

	bin_sol=Binsol(
		Dict{Any,Any}(sb=>value(sa[sb]) for sb=seg_bid),
		Dict{Any,Any}(bd=>value(qa[bd]) for bd=bid),
		Dict{Any,Any}(qt=>value(qta[qt]) for qt=qt_bid),
		Dict{Any,Any}(sb=>value(alpha[sb]) for sb=AlphaOmegaIndex),
		Dict{Any,Any}(sb=>value(omega[sb]) for sb=AlphaOmegaIndex)
	)
	ComputeGaps(socp_sol)
	println("SW: $(SW(socp_sol))")
	MaxSlack=(maximum(abs,union(value.(slackPUp), value.(slackQUp),value.(slackPDown), value.(slackQDown))))
	SumSlack=(sum(value.(slackPUp)) + sum(value.(slackPDown)) + sum(value.(slackQUp)) + sum(value.(slackQDown)))
	println("Max Slack: $MaxSlack")
	println("Sum Slack: $SumSlack")
	return (socp_sol,bin_sol)
end

function solveDCSOCP(solver::Symbol=:Mosek,binsol=0)::Tuple{PFsol,Binsol,Dualsol,Any}
	m=Model(optimizer_with_attributes(Mosek.Optimizer,
	"MSK_IPAR_PRESOLVE_USE" => 1,
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
	"MSK_DPAR_OPTIMIZER_MAX_TIME" => 3600
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
	@variable(m, 0 <= c[e=c_index,times])
	@variable(m, s[s_index,times])
	@variable(m, pg[sb=seg_bid])
	@expression(m, net_pg_bid[bd=bid], sum(pg[sb] for sb=bid_2_seg_bids[bd]))
	@variable(m, qg[sb=dist_seg_bid])
	@expression(m, net_qg_bid[bd=dist_bid], sum(qg[sb] for sb=bid_2_seg_bids[bd]))
	@variable(m, -pi <= theta_tr[union(tr_nodes_ext,slack_bus),times] <= pi)

	### BINARY VARIABLES ###
	@variable(m, 0 <= sa[seg_bid] <= 1)
	@variable(m, 0 <= qa[bid] <= 1)
	@variable(m, 0 <= qta[qt_bid] <= 1)
	@variable(m, 0 <= alpha[AlphaOmegaIndex] <= 1)
	@variable(m, 0 <= omega[AlphaOmegaIndex] <= 1)
	if binsol!=0
		for sb=seg_bid fix(sa[sb],binsol.sa[sb], force=true) end
		for sb=bid fix(qa[sb],binsol.qa[sb], force=true) end
		for sb=qt_bid fix(qta[sb],binsol.qta[sb], force=true) end
		for sb=AlphaOmegaIndex fix(alpha[sb],binsol.alpha[sb], force=true) end
		for sb=AlphaOmegaIndex fix(omega[sb],binsol.omega[sb], force=true) end
	end

	### SLACK VARIABLES ###
	@variable(m, slackPUp[nodes,times] == 0)
	@variable(m, slackPDown[nodes,times] == 0)
	@variable(m, slackQUp[dist_nodes,times] == 0)
	@variable(m, slackQDown[dist_nodes,times] == 0)

	### FLOW DEFINITION ###
	@variable(m, fp[edges,times])
	@variable(m, fpDist[ext_edges,times])
	@variable(m,fpTo[edges,times])
	@constraint(m, DefFpTr[e=tr_edges,t=times],
	X_pu[e]*fp[e,t]== (theta_tr[e_fr[e],t]-theta_tr[e_to[e],t]))
	@constraint(m, DefFpX0[e=tr_edges,t=times],
	fp[e,t] + fpTo[e,t] == 0)
	# We add this constraint to obtain the interface prices #
	@constraint(m, InterfaceConstr[e=ext_edges,t=times], fp[e,t] == fpDist[e,t])

	DefFpDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	DefFpToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
	for l=dist_edges for t=times
		if R_pu[l] == 0.0
			get!(DefFpDist,(l,t), @constraint(m, fp[l,t] + fpTo[l,t] == 0))
			@constraint(m, c[(e_fr[l],e_to[l]),t] ==  c[(e_to[l],e_to[l]),t])
			@constraint(m, c[(e_fr[l],e_fr[l]),t] ==  c[(e_fr[l],e_to[l]),t])
			@constraint(m, s[(e_fr[l],e_to[l]),t] == 0)
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
	@constraint(m, PDownSA[sb=seg_bid], -pg[sb] <= -sa[sb]*p_min[sb])
	@constraint(m, PUpSA[sb=seg_bid], pg[sb] <= sa[sb]*p_max[sb])

	@constraint(m, TrPowerFlowLimitUp[e=intersect(line_limit,tr_edges),t=times],
	fp[e,t] <= rate_a[e])

	@constraint(m, TrPowerFlowLimitDown[e=intersect(line_limit,tr_edges),t=times],
	-fp[e,t] <= rate_a[e])

	for i=slack_bus for t=times fix(theta_tr[i,t],0;force=true) end end
	#for i=slack_bus for t=times @constraint(m,theta_tr[i,t]==0) end end

	PowerBalanceP=Dict{Tuple{Int64,Int64},ConstraintRef}()
	for i=tr_nodes for t=times
		PowerBalanceP[i,t]=@constraint(m,
		sum(sum(pg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
			for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
		+ net_injP[i,t]
		+ slackPUp[i,t]
		==
		+ sum(fpTo[e,t] for e=n_to[i])
		+ sum(fp[e,t] for e=n_fr[i])
		+ slackPDown[i,t]
		)
	end end
	for i=dist_nodes for t=times
		PowerBalanceP[i,t]=@constraint(m,
		sum(sum(pg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
			for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
		+ net_injP[i,t]
		+ slackPUp[i,t]
		==
		+ sum(fpTo[e,t] for e=n_to[i])
		+ sum(fp[e,t] for e=n_fr[i] if e in dist_edges)
		+ sum(fpDist[e,t] for e=n_fr[i] if e in ext_edges)
		+ slackPDown[i,t]
		)
	end end

	@constraint(m, PowerBalanceQ[i=dist_nodes,t=times],
	sum(sum(qg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
	+ net_injQ[i,t]
	+ slackQUp[i,t]
	==
	+ sum(fqTo[e,t] for e=dist_n_to[i])
	+ sum(fq[e,t] for e=dist_n_fr[i])
	+ slackQDown[i,t]
	)

	# ### tentative constraint ###
	# @constraint(m, TotalPBalance[t=times],
	# + sum(pg[sb] for sb=seg_bid if sb[5]==t)
	# + sum(net_injP[i,t] for i=nodes)
	# + sum(slackPUp[i,t] for i=nodes)
	# ==
	# + sum(fpTo[e,t] for e=edges)
	# + sum(fp[e,t] for e=edges)
	# + sum(slackPDown[i,t] for i=nodes)
	# )
	# ###
	#
	# ### tentative constraint ###
	# @constraint(m, TotalQBalance[t=times],
	# + sum(qg[sb] for sb=dist_seg_bid if sb[5]==t)
	# + sum(net_injQ[i,t] for i=dist_nodes)
	# + sum(slackQUp[i,t] for i=dist_nodes)
	# ==
	# + sum(fqTo[e,t] for e=dist_edges)
	# + sum(fq[e,t] for e=dist_edges)
	# + sum(slackQDown[i,t] for i=dist_nodes)
	# )
	# ###


	@constraint(m, VoltageLimitsDown[i=v_limit,t=times],
	-c[(i,i),t] <= -vsq_min[i])

	@constraint(m, VoltageLimitsUp[i=v_limit,t=times],
	c[(i,i),t] <= vsq_max[i])

	mod_s_ind=[(e_fr[e],e_to[e]) for e=dist_edges if R_pu[e]!=0.0]
	if solver==:Gurobi
		@constraint(m, LinePowerLimit1[e=intersect(line_limit,dist_edges),t=times],
		fp[e,t]^2+fq[e,t]^2 <= rate_a[e]^2)

		@constraint(m, LinePowerLimit2[e=intersect(line_limit,dist_edges),t=times],
		fpTo[e,t]^2+fqTo[e,t]^2 <= rate_a[e]^2)

		@constraint(m, SOCPconstraint[e=mod_s_ind,t=times],
		c[e,t]^2+s[e,t]^2 <= c[(e[1],e[1]),t]*c[(e[2],e[2]),t])

		@constraint(m, DiscQtoP[BQpd=QPDiscConstrIndex],
		get(QPDiscMax,(BQpd[2],BQpd[3],BQpd[5]),0)^2
		- net_pg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]^2
		- net_qg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]^2 >= 0
		)
	else
		@constraint(m, LinePowerLimit1[e=intersect(line_limit,dist_edges),t=times],
		[rate_a[e],fp[e,t],fq[e,t]] in SecondOrderCone())

		@constraint(m, LinePowerLimit2[e=intersect(line_limit,dist_edges),t=times],
		[rate_a[e],fpTo[e,t],fqTo[e,t]] in SecondOrderCone())

		@constraint(m, SOCPconstraint[e=mod_s_ind,t=times],
		[0.5*c[(e[1],e[1]),t],c[(e[2],e[2]),t],c[e,t],s[e,t]] in RotatedSecondOrderCone())

		@constraint(m, DiscQtoP[BQpd=QPDiscConstrIndex],
		[get(QPDiscMax,(BQpd[2],BQpd[3],BQpd[5]),0),
		net_pg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])],
		net_qg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]]
		in SecondOrderCone()
		)
	end

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

	@constraint(m, ConstrainToUpFromLineQtoP[BHpc=QPHPConstrIndex1],
	net_qg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	-QPHPOffset[BHpc[5],BHpc[2]]
	-QPHPSlope[BHpc[5],BHpc[2]]*net_pg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	>= 0)

	@constraint(m, ConstrainToDownFromLineQtoP[BHpc=QPHPConstrIndex0],
	QPHPOffset[BHpc[5],BHpc[2]]
	+QPHPSlope[BHpc[5],BHpc[2]]*net_pg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
	-net_qg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])] >= 0)

	### BINARY CONSTRAINTS ###
	if binsol!=0
		@constraint(m, ExQtBidsConstr[i=ex_qt_bid_id],
		sum(qta[ExQt] for ExQt=ex_qt_bid[i]) <= 1)

		@constraint(m, QtaDefLow[qt=qt_bid],
		qta[qt] <= sum(qa[bd] for bd=bid if bd[2]==qt))

		@constraint(m, QtaDefHigh[bd=bid],
		qa[bd] <= qta[bd[2]])

		@constraint(m, DefQa[sb=seg_bid],
		sa[sb] <= qa[(sb[1],sb[2],sb[3],sb[5])])

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
	end

	### OBJECTIVE ###
	@objective(m, Min,
	+ sum(cost[sb][1]*pg[sb] for sb=seg_bid)
	+ cost_ls*(sum(slackPUp[i,t] + slackPDown[i,t] for i=nodes for t=times)
		 + sum(slackQUp[i,t] + slackQDown[i,t] for i=dist_nodes for t=times))
	)

	optimize!(m)
	println(termination_status(m))
	println(primal_status(m))

	pgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(pg[j]) for j=seg_bid)
	qgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(qg[j]) for j=dist_seg_bid)
	cHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((i,t) => value(c[i,t]) for i=c_index for t=times)
	sHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((l,t) => value(s[l,t]) for l=s_index for t=times)
	theta_trHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(theta_tr[i,t]) for i=union(slack_bus,tr_nodes_ext) for t=times)
	merge!(theta_trHist,Dict{Tuple{Int64,Int64},Float64}((i,t) => 0.0 for i=setdiff(nodes,tr_nodes_ext) for t=times))
	fpHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fp[l,t]) for l=edges for t=times)
	fqHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fq[l,t]) for l=dist_edges for t=times)
	fpToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fpTo[l,t]) for l=edges for t=times)
	fqToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fqTo[l,t]) for l=dist_edges for t=times)
	slackPUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPUp[i,t]) for i=nodes for t=times)
	slackPDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPDown[i,t]) for i=nodes for t=times)
	slackQUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackQUp[i,t]) for i=dist_nodes for t=times)
	slackQDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackQDown[i,t]) for i=dist_nodes for t=times)
	socp_sol=PFsol(false,pgHist,qgHist,cHist,sHist,theta_trHist,fpHist,fqHist,
	fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)

	socp_sol=heurTheta2(socp_sol)

	bin_sol=Binsol(
		Dict{Any,Any}(sb=>value(sa[sb]) for sb=seg_bid),
		Dict{Any,Any}(bd=>value(qa[bd]) for bd=bid),
		Dict{Any,Any}(qt=>value(qta[qt]) for qt=qt_bid),
		Dict{Any,Any}(sb=>value(alpha[sb]) for sb=AlphaOmegaIndex),
		Dict{Any,Any}(sb=>value(omega[sb]) for sb=AlphaOmegaIndex)
	)

	dual_sol=Dualsol(0,0)
	prices_interface=Dict{Tuple{Int64,Int64},Float64}()
	if has_duals(m) && solver != :Gurobi
		priceP=Dict{Tuple{Int64,Int64},Float64}((i,t) => dual(PowerBalanceP[i,t]) for i=nodes for t=times)
		priceQ=Dict{Tuple{Int64,Int64},Float64}((i,t) => dual(PowerBalanceQ[i,t]) for i=dist_nodes for t=times)
		prices_interface=Dict{Tuple{Int64,Int64},Float64}((e,t) => dual(InterfaceConstr[e,t]) for e=ext_edges for t=times)
		dual_sol=Dualsol(priceP,priceQ)
	else
		println(">> Dual solution not computed.")
	end

	ComputeGaps(socp_sol)
	println("SW: $(SW(socp_sol))")
	MaxSlack=(maximum(abs,union(value.(slackPUp), value.(slackQUp),value.(slackPDown), value.(slackQDown))))
	SumSlack=(sum(value.(slackPUp)) + sum(value.(slackPDown)) + sum(value.(slackQUp)) + sum(value.(slackQDown)))
	println("Max Slack: $MaxSlack")
	println("Sum Slack: $SumSlack")
	return (socp_sol,bin_sol,dual_sol,prices_interface)
end

function solveDCAC_bin(initsol,binsol,warmstart::String="yes")::Tuple{PFsol,Dualsol}
 	m=Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>warmstart,
	"bound_frac" => 1e-30,
	"bound_push" => 1e-30,
	"slack_bound_frac" => 1e-30,
	"slack_bound_push" => 1e-30,
	"bound_mult_init_method" => "mu-based",
	"mu_init" => 0.01,
	"max_cpu_time" => 3.6e+3
	))
	### VARIABLES ###
	@variable(m, 0 <= c[c_index,times])
	@variable(m, s[s_index,times])
	@variable(m, pg[sb=seg_bid])
	@expression(m, net_pg_bid[bd=bid], sum(pg[sb] for sb=bid_2_seg_bids[bd]))
	@variable(m, injP[nodes,times])
	@variable(m, qg[sb=dist_seg_bid])
	@expression(m, net_qg_bid[bd=dist_bid], sum(qg[sb] for sb=bid_2_seg_bids[bd]))
	@variable(m, injQ[dist_nodes,times])
	@variable(m, -pi <= theta[nodes,times] <= pi)
	@variable(m, slackPUp[nodes,times] >= 0)
	@variable(m, slackPDown[nodes,times] >= 0)
	@variable(m, slackQUp[dist_nodes,times] >= 0)
	@variable(m, slackQDown[dist_nodes,times] >= 0)

	### FLOW DEFINITION ###
	@variable(m, fp[edges,times])
	@variable(m,fpTo[edges,times])
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
			#fix(s[(e_fr[l],e_to[l]),t],0.0;force=true)
			@constraint(m,s[(e_fr[l],e_to[l]),t] == 0.0)
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
	@constraint(m, PDown[sb=seg_bid], -pg[sb] <= -binsol.sa[sb]*p_min[sb])
	@constraint(m, PUp[sb=seg_bid], pg[sb] <= binsol.sa[sb]*p_max[sb])

	@constraint(m, TrPowerFlowLimitUp[e=intersect(line_limit,tr_edges),t=times],
	fp[e,t] <= rate_a[e])

	@constraint(m, TrPowerFlowLimitDown[e=intersect(line_limit,tr_edges),t=times],
	-fp[e,t] <= rate_a[e])

	@constraint(m, PowerBalanceP[i=nodes,t=times],
	sum(sum(pg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
	+ net_injP[i,t]
	+ slackPUp[i,t]
	==
	+ sum(fpTo[e,t] for e=n_to[i])
	+ sum(fp[e,t] for e=n_fr[i])
	+ slackPDown[i,t]
	)

	@constraint(m, PowerBalanceQ[i=dist_nodes,t=times],
	sum(sum(qg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
	+ net_injQ[i,t]
	+ slackQUp[i,t]
	==
	+ sum(fqTo[e,t] for e=dist_n_to[i])
	+ sum(fq[e,t] for e=dist_n_fr[i])
	+ slackQDown[i,t]
	)

	@constraint(m, VoltageLimitsDown[i=v_limit,t=times],
	-c[(i,i),t] <= -vsq_min[i])

	@constraint(m, VoltageLimitsUp[i=v_limit,t=times],
	c[(i,i),t] <= vsq_max[i])

	@constraint(m, LinePowerLimit1[e=intersect(line_limit,dist_edges),t=times],
	fp[e,t]^2 + fq[e,t]^2 <= rate_a[e]^2)

	@constraint(m, LinePowerLimit2[e=intersect(line_limit,dist_edges),t=times],
	fpTo[e,t]^2 + fqTo[e,t]^2 <= rate_a[e]^2)

	mod_s_ind=[(e_fr[e],e_to[e]) for e=dist_edges if R_pu[e]!=0.0]
	@NLconstraint(m, Quadconstraint[e=mod_s_ind,t=times],
	c[e,t]^2 + s[e,t]^2 == c[(e[1],e[1]),t]*c[(e[2],e[2]),t])

	#@NLconstraint(m, Tanconstraint[e=mod_s_ind,t=times],
	#c[e,t]*sin(theta[e[1],t]-theta[e[2],t])
	#+ s[e,t]*cos(theta[e[1],t]-theta[e[2],t]) == 0)

	@constraint(m, FixTheta[i=slack_bus,t=times],theta[i,t]==0)

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
	get(QPDiscMax,(BQpd[2],BQpd[3],BQpd[5]),0)^2
	- net_pg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]^2
	- net_qg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]^2 >= 0
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

	### OBJECTIVE ###
	@objective(m, Min,(
	+ sum(cost[sb][1]*pg[sb] for sb=seg_bid)
	+ cost_ls*(sum(slackPUp[i,t]+ slackPDown[i,t] for i=nodes for t=times)
		 + sum(slackQUp[i,t] + slackQDown[i,t] for i=dist_nodes for t=times))
	))

	if warmstart=="yes"
		println(" >> WARMSTARTING <<")
		for sb=seg_bid set_start_value(pg[sb],initsol.pg[sb]) end
		for sb=dist_seg_bid set_start_value(qg[sb],initsol.qg[sb]) end
		for t=times
			for l=c_index set_start_value(c[l,t],initsol.c[l,t]) end
			for i=nodes
				set_start_value(theta[i,t],initsol.theta[i,t])
				set_start_value(slackPUp[i,t],0.0)
				set_start_value(slackPDown[i,t],0.0)
			end
			for i=dist_nodes
				set_start_value(slackQUp[i,t],0.0)
				set_start_value(slackQDown[i,t],0.0)
			end
			for l=s_index set_start_value(s[l,t],initsol.s[l,t]) end
			for e=edges
				set_start_value(fp[e,t],initsol.fp[e,t])
				set_start_value(fpTo[e,t],initsol.fpTo[e,t])
			end
			for e=dist_edges
				set_start_value(fq[e,t],initsol.fq[e,t])
				set_start_value(fqTo[e,t],initsol.fqTo[e,t])
			end
		end
	end

	optimize!(m)

	pgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(pg[j]) for j=seg_bid)
	qgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(qg[j]) for j=dist_seg_bid)
	cHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((i,t) => value(c[i,t]) for i=c_index for t=times)
	sHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((l,t) => value(s[l,t]) for l=s_index for t=times)
	thetaHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(theta[i,t]) for i=nodes for t=times)
	fpHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fp[l,t]) for l=edges for t=times)
	fqHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fq[l,t]) for l=dist_edges for t=times)
	fpToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fpTo[l,t]) for l=edges for t=times)
	fqToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fqTo[l,t]) for l=dist_edges for t=times)
	priceP=Dict{Tuple{Int64,Int64},Float64}((i,t) => dual(PowerBalanceP[i,t]) for i=nodes for t=times)
	priceQ=Dict{Tuple{Int64,Int64},Float64}((i,t) => dual(PowerBalanceQ[i,t]) for i=dist_nodes for t=times)
	slackPUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPUp[i,t]) for i=nodes for t=times)
	slackPDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPDown[i,t]) for i=nodes for t=times)
	slackQUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackQUp[i,t]) for i=dist_nodes for t=times)
	slackQDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackQDown[i,t]) for i=dist_nodes for t=times)
	MaxSlack=(maximum(abs,union(value.(slackPUp), value.(slackQUp),value.(slackPDown), value.(slackQDown))))
	println("Max Slack: $MaxSlack")
	pfsol=PFsol(MaxSlack>1e-5,pgHist,qgHist,cHist,sHist,thetaHist,
	fpHist,fqHist,fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
	dual_sol=Dualsol(priceP,priceQ)
	ComputeGaps(pfsol)
	println("SW: $(SW(pfsol))")
	return (pfsol,dual_sol)
end

function solveDC_TN(sol::PFsol,
					binsol,
					prices_int::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}(),
					solver::Symbol=:Mosek)
	m=Model(optimizer_with_attributes(Mosek.Optimizer,
	"MSK_IPAR_PRESOLVE_USE" => 0,
	"MSK_IPAR_LOG" => 0,
	"MSK_IPAR_INTPNT_BASIS" => 0,
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
	@variable(m, slackPUp[tr_nodes,times] == 0)
	@variable(m, slackPDown[tr_nodes,times] == 0)

	### FLOW DEFINITION ###
	@variable(m, fp[tr_edges,times])
	@constraint(m, DefFpTr[e=tr_edges,t=times],
	X_pu[e]*fp[e,t] == (theta_tr[e_fr[e],t]-theta_tr[e_to[e],t]))

	@variable(m,fpTo[tr_edges,times])
	@constraint(m, DefFpX0[e=tr_edges,t=times],
	fp[e,t] + fpTo[e,t] == 0)

	### CONSTRAINTS ###
	if binsol!=0
		@constraint(m, PDownSA[sb=tr_seg_bid], -pg[sb] <= -binsol.sa[sb]*p_min[sb])
		@constraint(m, PUpSA[sb=tr_seg_bid], pg[sb] <= binsol.sa[sb]*p_max[sb])
	else
		@variable(m, 0 <= relax_sa[tr_seg_bid] <= 1)
	    @constraint(m, PDown[sb=tr_seg_bid], -pg[sb] <= -relax_sa[sb]*p_min[sb])
	    @constraint(m, PUp[sb=tr_seg_bid], pg[sb] <= relax_sa[sb]*p_max[sb])
	end

	@constraint(m, TrPowerFlowLimitUp[e=intersect(line_limit,tr_edges),t=times],
	fp[e,t] <= rate_a[e])

	@constraint(m, TrPowerFlowLimitDown[e=intersect(line_limit,tr_edges),t=times],
	-fp[e,t] <= rate_a[e])

	for i=slack_bus for t=times @constraint(m, theta_tr[i,t]==0) end end

	@constraint(m, PowerBalanceP[i=tr_nodes,t=times],
	sum(sum(pg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
	+ net_injP[i,t]
	+ slackPUp[i,t]
	==
	+ sum(fpTo[e,t] for e=n_to[i])
	+ sum(fp[e,t] for e=n_fr[i])
	+ slackPDown[i,t]
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

	### for the results of the thesis and the transaction paper, this constraint
	# was placed here which is probably wrong when computing the prices
	# is corrected for now. Provide better results for the small Italian test case
	@constraint(m,thetaFix[i=root_dist_node,t=times],theta_tr[i,t]-sol.theta[i,t] == 0)

	### OBJECTIVE ###
	@objective(m, Min,
	+ sum(cost[sb][1]*pg[sb] for sb=tr_seg_bid)
	+ cost_ls*sum(slackPUp[i,t] + slackPDown[i,t] for i=tr_nodes for t=times)
	)
	if isempty(prices_int)
		@constraint(m,fpFix[e=ext_edges,t=times], fp[e,t]-sol.fp[e,t] == 0)
		### for the results of the thesis and the transaction paper, this constraint
		# was placed ABOVE which is probably wrong when computing the prices.
		#@constraint(m,thetaFix[i=root_dist_node,t=times],theta_tr[i,t]-sol.theta[i,t] == 0)
	else
		@objective(m, Min,
		+ sum(cost[sb][1]*pg[sb] for sb=tr_seg_bid)
		- sum(prices_int[e,t]*fp[e,t] for e=ext_edges for t=times)
		+ cost_ls*sum(slackPUp[i,t] + slackPDown[i,t] for i=tr_nodes for t=times)
		)
	end

	#println(m)
	optimize!(m)
	@assert termination_status(m)==MOI.OPTIMAL
	@assert primal_status(m)==MOI.FEASIBLE_POINT

	pgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(pg[j]) for j=tr_seg_bid)
	qgHist = 0
	cHist = 0
	sHist = 0
	theta_trHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(theta_tr[i,t]) for i=union(slack_bus,tr_nodes_ext) for t=times)
	fpHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fp[l,t]) for l=tr_edges for t=times)
	fqHist = 0
	fpToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fpTo[l,t]) for l=tr_edges for t=times)
	fqToHist = 0
	slackPUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPUp[i,t]) for i=tr_nodes for t=times)
	slackPDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPDown[i,t]) for i=tr_nodes for t=times)
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
	new_prices_int=copy(prices_int)
	if isempty(prices_int) for e=ext_edges for t=times
			new_prices_int[e,t]=dual(fpFix[e,t])
	end end end
	return (dc_sol,dual_sol,new_prices_int,solve_time(m))
end

function solveDC_pert_interface(exp_flow,disturb=0.1,solver::Symbol=:Mosek)
	m=Model(optimizer_with_attributes(Mosek.Optimizer,
	"MSK_IPAR_PRESOLVE_USE" => 0,
	"MSK_IPAR_LOG" => 0,
	"MSK_IPAR_INTPNT_BASIS" => 0,
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
	@variable(m, slackPUp[tr_nodes,times] >= 0)
	@variable(m, slackPDown[tr_nodes,times] >= 0)

	### BINARY VARIABLES ###
	@variable(m, sa[tr_seg_bid], Bin)
	@variable(m, qa[tr_bid], Bin)
	@variable(m, qta[tr_qt_bid], Bin)
	@variable(m, alpha[sb=[aoi for aoi=AlphaOmegaIndex if aoi[1] in tr_nodes]], Bin)
	@variable(m, omega[sb=[aoi for aoi=AlphaOmegaIndex if aoi[1] in tr_nodes]], Bin)

	### FLOW DEFINITION ###
	@variable(m, fp[tr_edges,times])
	@constraint(m, DefFpTr[e=tr_edges,t=times],
	X_pu[e]*fp[e,t] == (theta_tr[e_fr[e],t]-theta_tr[e_to[e],t]))

	@variable(m,fpTo[tr_edges,times])
	@constraint(m, DefFpX0[e=tr_edges,t=times],
	fp[e,t] + fpTo[e,t] == 0)

	### CONSTRAINTS ###
	@constraint(m, PDownSA[sb=tr_seg_bid], -pg[sb] <= -sa[sb]*p_min[sb])
	@constraint(m, PUpSA[sb=tr_seg_bid], pg[sb] <= sa[sb]*p_max[sb])

	@constraint(m, TrPowerFlowLimitUp[e=intersect(line_limit,tr_edges),t=times],
	fp[e,t] <= rate_a[e])

	@constraint(m, TrPowerFlowLimitDown[e=intersect(line_limit,tr_edges),t=times],
	-fp[e,t] <= rate_a[e])

	for i=slack_bus for t=times @constraint(m, theta_tr[i,t]==0) end end

	@constraint(m, PowerBalanceP[i=tr_nodes,t=times],
	sum(sum(pg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
	+ net_injP[i,t]
	+ slackPUp[i,t]
	==
	+ sum(fpTo[e,t] for e=n_to[i])
	+ sum(fp[e,t] for e=n_fr[i])
	+ slackPDown[i,t]
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

	@constraint(m,fpFix[e=ext_edges,t=times], fp[e,t] == (1-sign(exp_flow[e,t])*disturb)*exp_flow[e,t])

	@constraint(m, ExQtBidsConstr[i=tr_ex_qt_bid_id],
	sum(qta[ExQt] for ExQt=ex_qt_bid[i]) <= 1)

	@constraint(m, QtaDefLow[qt=tr_qt_bid],
	qta[qt] <= sum(qa[bd] for bd=bid if bd[2]==qt))

	@constraint(m, QtaDefHigh[bd=tr_bid],
	qa[bd] <= qta[bd[2]])

	@constraint(m, DefQa[sb=tr_seg_bid],
	sa[sb] <= qa[(sb[1],sb[2],sb[3],sb[5])])

	@constraint(m, AlphaOmegaDef[B=[aol for aol=AlphaOmegaLate if aol[1] in tr_nodes]],
	qa[B]
	- sum(qa[bd] for bd=bid if bd[1]==B[1] && bd[1]==B[2] && bd[4]==B[4]-1)
	- alpha[B]
	+ omega[B]
	== 0)

	@constraint(m, MutexAlphaOmega[B=[aoi for aoi=AlphaOmegaIndex if aoi[1] in tr_nodes]],
	alpha[B] + omega[B] <= 1)

	@constraint(m, AlphaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in tr_nodes]],
	qa[B] - alpha[B] == 0)

	@constraint(m, OmegaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in tr_nodes]],
	omega[B] == 0)

	@constraint(m, MinDurationConstr[mdp=tr_MinDurIndex],
	+ sum(qa[bd] for bd=bid if bd[2]==mdp[1] && bd[3]==mdp[3] && bd[4]==mdp[5])
		- sum(alpha[bd] for bd=bid if bd[2]==mdp[1] && bd[3]==mdp[2] && bd[4]==mdp[4])
	>= 0)

	@constraint(m, NoNewAct[B=[nnai for nnai=NoNewActIndex if nnai[1] in tr_nodes]],
	alpha[B] == 0)

	### OBJECTIVE ###
	@objective(m, Min,
	+ sum(cost[sb][1]*pg[sb] for sb=tr_seg_bid)
	+ cost_ls*sum(slackPUp[i,t] + slackPDown[i,t] for i=tr_nodes for t=times)
	)

	optimize!(m)
	@assert termination_status(m)==MOI.OPTIMAL
	@assert primal_status(m)==MOI.FEASIBLE_POINT

	pgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(pg[j]) for j=tr_seg_bid)
	qgHist = 0
	cHist = 0
	sHist = 0
	theta_trHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(theta_tr[i,t]) for i=union(slack_bus,tr_nodes_ext) for t=times)
	fpHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fp[l,t]) for l=tr_edges for t=times)
	fqHist = 0
	fpToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fpTo[l,t]) for l=tr_edges for t=times)
	fqToHist = 0
	slackPUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPUp[i,t]) for i=tr_nodes for t=times)
	slackPDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPDown[i,t]) for i=tr_nodes for t=times)
	slackQUpHist = 0
	slackQDownHist = 0
	dc_sol=PFsol(false,pgHist,qgHist,cHist,sHist,theta_trHist,fpHist,fqHist,
	fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
	return dc_sol
end

function solveDC_no_coordination(solver::Symbol=:Mosek)
	m=Model(optimizer_with_attributes(Mosek.Optimizer,
	"MSK_IPAR_PRESOLVE_USE" => 0,
	"MSK_IPAR_LOG" => 0,
	"MSK_IPAR_INTPNT_BASIS" => 0,
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
	@variable(m, slackPUp[tr_nodes,times] >= 0)
	@variable(m, slackPDown[tr_nodes,times] >= 0)

	### BINARY VARIABLES ###
	@variable(m, sa[tr_seg_bid], Bin)
	@variable(m, qa[tr_bid], Bin)
	@variable(m, qta[tr_qt_bid], Bin)
	@variable(m, alpha[sb=[aoi for aoi=AlphaOmegaIndex if aoi[1] in tr_nodes]], Bin)
	@variable(m, omega[sb=[aoi for aoi=AlphaOmegaIndex if aoi[1] in tr_nodes]], Bin)

	### FLOW DEFINITION ###
	@variable(m, fp[tr_edges,times])
	@constraint(m, DefFpTr[e=tr_edges,t=times],
	X_pu[e]*fp[e,t] == (theta_tr[e_fr[e],t]-theta_tr[e_to[e],t]))

	@variable(m,fpTo[tr_edges,times])
	@constraint(m, DefFpX0[e=tr_edges,t=times],
	fp[e,t] + fpTo[e,t] == 0)

	### CONSTRAINTS ###
	@constraint(m, PDownSA[sb=tr_seg_bid], -pg[sb] <= -sa[sb]*p_min[sb])
	@constraint(m, PUpSA[sb=tr_seg_bid], pg[sb] <= sa[sb]*p_max[sb])

	@constraint(m, TrPowerFlowLimitUp[e=intersect(line_limit,tr_edges),t=times],
	fp[e,t] <= rate_a[e])

	@constraint(m, TrPowerFlowLimitDown[e=intersect(line_limit,tr_edges),t=times],
	-fp[e,t] <= rate_a[e])

	for i=slack_bus for t=times @constraint(m, theta_tr[i,t]==0) end end

	@constraint(m, PowerBalanceP[i=tr_nodes,t=times],
	sum(sum(pg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
	+ net_injP[i,t]
	+ slackPUp[i,t]
	==
	+ sum(fpTo[e,t] for e=n_to[i])
	+ sum(fp[e,t] for e=n_fr[i])
	+ slackPDown[i,t]
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

	@constraint(m,fpFix[e=ext_edges,t=times], fp[e,t] == 0)#(1-sign(exp_flow[e,t])*disturb)*exp_flow[e,t])

	@constraint(m, ExQtBidsConstr[i=tr_ex_qt_bid_id],
	sum(qta[ExQt] for ExQt=ex_qt_bid[i]) <= 1)

	@constraint(m, QtaDefLow[qt=tr_qt_bid],
	qta[qt] <= sum(qa[bd] for bd=bid if bd[2]==qt))

	@constraint(m, QtaDefHigh[bd=tr_bid],
	qa[bd] <= qta[bd[2]])

	@constraint(m, DefQa[sb=tr_seg_bid],
	sa[sb] <= qa[(sb[1],sb[2],sb[3],sb[5])])

	@constraint(m, AlphaOmegaDef[B=[aol for aol=AlphaOmegaLate if aol[1] in tr_nodes]],
	qa[B]
	- sum(qa[bd] for bd=bid if bd[1]==B[1] && bd[1]==B[2] && bd[4]==B[4]-1)
	- alpha[B]
	+ omega[B]
	== 0)

	@constraint(m, MutexAlphaOmega[B=[aoi for aoi=AlphaOmegaIndex if aoi[1] in tr_nodes]],
	alpha[B] + omega[B] <= 1)

	@constraint(m, AlphaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in tr_nodes]],
	qa[B] - alpha[B] == 0)

	@constraint(m, OmegaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in tr_nodes]],
	omega[B] == 0)

	@constraint(m, MinDurationConstr[mdp=tr_MinDurIndex],
	+ sum(qa[bd] for bd=bid if bd[2]==mdp[1] && bd[3]==mdp[3] && bd[4]==mdp[5])
		- sum(alpha[bd] for bd=bid if bd[2]==mdp[1] && bd[3]==mdp[2] && bd[4]==mdp[4])
	>= 0)

	@constraint(m, NoNewAct[B=[nnai for nnai=NoNewActIndex if nnai[1] in tr_nodes]],
	alpha[B] == 0)

	### OBJECTIVE ###
	@objective(m, Min,
	+ sum(cost[sb][1]*pg[sb] for sb=tr_seg_bid)
	+ cost_ls*sum(slackPUp[i,t] + slackPDown[i,t] for i=tr_nodes for t=times)
	)

	optimize!(m)
	@assert termination_status(m)==MOI.OPTIMAL
	@assert primal_status(m)==MOI.FEASIBLE_POINT

	pgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(pg[j]) for j=tr_seg_bid)
	qgHist = 0
	cHist = 0
	sHist = 0
	theta_trHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(theta_tr[i,t]) for i=union(slack_bus,tr_nodes_ext) for t=times)
	fpHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fp[l,t]) for l=tr_edges for t=times)
	fqHist = 0
	fpToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fpTo[l,t]) for l=tr_edges for t=times)
	fqToHist = 0
	slackPUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPUp[i,t]) for i=tr_nodes for t=times)
	slackPDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPDown[i,t]) for i=tr_nodes for t=times)
	slackQUpHist = 0
	slackQDownHist = 0
	dc_sol=PFsol(false,pgHist,qgHist,cHist,sHist,theta_trHist,fpHist,fqHist,
	fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
	return dc_sol
end

function solveAC_subnet(dn::String,
						sol::PFsol,
						binsol=0,
						prices_int::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}(),
						warmstart::String="no")::Tuple{PFsol,Dualsol,Float64}
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
	m=Model(optimizer_with_attributes(Ipopt.Optimizer,
	"warm_start_init_point"=>warmstart,
	# "bound_frac" => 1e-5,
	# "bound_push" => 1e-5,
	# "slack_bound_frac" => 1e-5,
	# "slack_bound_push" => 1e-5,
	#"bound_mult_init_method" => "mu-based",
	"mu_init" => 1e-4,
	"print_level" => 0,
	#"tol" => 1e-7,
	# "dual_inf_tol" => 1e-7,
	# "constr_viol_tol" => 1e-7,
	# "compl_inf_tol" => 1e-7,
	))
	if !isempty(prices_int)
	m=Model(optimizer_with_attributes(Ipopt.Optimizer,
	"warm_start_init_point"=>warmstart,
	#"bound_frac" => 1e-5,
	#"bound_push" => 1e-5,
	#"slack_bound_frac" => 1e-5,
	#"slack_bound_push" => 1e-5,
	#"bound_mult_init_method" => "mu-based",
	#"mu_init" => 1e-4,
	"print_level" => 0
	))
	end
	### VARIABLES ###
	@variable(m, 0 <= c[dn_c_index,times])
	@variable(m, s[dn_s_index,times])
	@variable(m, pg[dn_seg_bid])
	@expression(m, net_pg_bid[bd=dn_bid], sum(pg[sb] for sb=bid_2_seg_bids[bd]))
	@variable(m, qg[dn_seg_bid])
	@expression(m, net_qg_bid[bd=dn_bid], sum(qg[sb] for sb=bid_2_seg_bids[bd]))
	@variable(m, -pi <= theta[union(dn_nodes,dn_ext_nodes),times] <= pi)

	# @variable(m, slackPUp[dn_nodes,times] >= 0)
	# @variable(m, slackPDown[dn_nodes,times] >= 0)
	# @variable(m, slackQUp[dn_nodes,times] >= 0)
	# @variable(m, slackQDown[dn_nodes,times] >= 0)

	# bound_prod=Dict{Tuple{Int64,Int64},Float64}((i,t) => +net_injP[i,t]
	# 	for i=dn_nodes for t=times)
	# if binsol!=0
	# 	# for i=dn_nodes for t=times
	# 	# 	bound_prod[i,t]+=(mapreduce(sum,+,init=0,mapreduce(sum,+,init=0,p_min[sb]
	# 	# 		for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()) if binsol.sa[sb] == 1 && p_min[sb]==p_max[sb])
	# 	# 		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}())))
	# 	# end end
	# 	for sb=dn_seg_bid
	# 		if binsol.sa[sb]>0.5 && abs(p_min[sb]-p_max[sb]) < 1e-6
	# 			bound_prod[sb[1],sb[5]]+=p_min[sb]
	# 		end
	# 	end
	# # end
	# @variable(m, abs(net_injP[i,t]) >= slackPUp[i=dn_nodes,t=times] >= 0)
	# @variable(m, abs(net_injP[i,t]) >= slackPDown[i=dn_nodes,t=times] >= 0)
	# @variable(m, abs(net_injQ[i,t]) >= slackQUp[i=dn_nodes,t=times] >= 0)
	# @variable(m, abs(net_injQ[i,t]) >= slackQDown[i=dn_nodes,t=times] >= 0)
	#
	@variable(m, slackPUp[i=dn_nodes,t=times] >= 0)
	@variable(m, slackPDown[i=dn_nodes,t=times] >= 0)
	@variable(m, slackQUp[i=dn_nodes,t=times] >= 0)
	@variable(m, slackQDown[i=dn_nodes,t=times] >= 0)
	# for i=dn_nodes for t=times
	# 	if abs(bound_prod[i,t]) > 1e-6
	# 		if sign(bound_prod[i,t]) > 0
	# 			@constraint(m,slackPUp[i,t]==0)
	# 		else
	# 			@constraint(m,slackPDown[i,t]==0)
	# 		end
	# 	else
	# 		@constraint(m,slackPUp[i,t]==0)
	# 		@constraint(m,slackPDown[i,t]==0)
	# 	end
	# 	if abs(net_injQ[i,t]) > 1e-6
	# 		if sign(net_injQ[i,t]) > 0
	# 			@constraint(m,slackQUp[i,t]==0)
	# 		else
	# 			@constraint(m,slackQDown[i,t]==0)
	# 		end
	# 	else
	# 		@constraint(m,slackQUp[i,t]==0)
	# 		@constraint(m,slackQDown[i,t]==0)
	# 	end
	# end end

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
	if binsol!=0
		@constraint(m, PDown[sb=dn_seg_bid], -pg[sb] <= -binsol.sa[sb]*p_min[sb])
		@constraint(m, PUp[sb=dn_seg_bid], pg[sb] <= binsol.sa[sb]*p_max[sb])
	else
		@variable(m, 0 <= relax_sa[dn_seg_bid] <= 1)
		@constraint(m, PDown[sb=dn_seg_bid], -pg[sb] <= -relax_sa[sb]*p_min[sb])
		@constraint(m, PUp[sb=dn_seg_bid], pg[sb] <= relax_sa[sb]*p_max[sb])
	end

	@constraint(m, PowerBalanceP[i=dn_nodes,t=times],
	sum(sum(pg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
	+ net_injP[i,t]
	+ slackPUp[i,t]
	==
	+ sum(fpTo[e,t] for e=n_to[i])
	+ sum(fp[e,t] for e=n_fr[i])
	+ slackPDown[i,t]
	)

	@constraint(m, PowerBalanceQ[i=dn_nodes,t=times],
	sum(sum(qg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
		for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
	+ net_injQ[i,t]
	+ slackQUp[i,t]
	==
	+ sum(fqTo[e,t] for e=dist_n_to[i])
	+ sum(fq[e,t] for e=dist_n_fr[i])
	+ slackQDown[i,t]
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
	if !isempty(prices_int)
	@constraint(m, FlowLimUp[e=dn_ext_edges,t=times],
	fp[e,t] <= rate_a[e])
	@constraint(m, FlowLimDown[e=dn_ext_edges,t=times],
	-fp[e,t] <= rate_a[e])
	end
	@constraint(m, DefFpX0[e=dn_ext_edges, t=times],
	fp[e,t] + fpTo[e,t] == 0)

	# ### OBJECTIVE ###
	@objective(m, Min,
	+ sum(cost[sb][1]*pg[sb] for sb=dn_seg_bid)
	+ cost_ls*sum(slackPUp[i,t] + slackPDown[i,t]
			+ slackQUp[i,t] + slackQDown[i,t] for i=dn_nodes for t=times)
	)

	if isempty(prices_int)
		### Fixed transmission variables ###
		@constraint(m,fixTheta[i=dn_ext_nodes,t=times], theta[i,t]-sol.theta[i,t]==0)
		@constraint(m,fpFix[e=dn_ext_edges,t=times], fp[e,t]-sol.fp[e,t] == 0)
	else
		@objective(m, Min,
		+ sum(cost[sb][1]*pg[sb] for sb=dn_seg_bid)
		+ cost_ls*sum(slackPUp[i,t] + slackPDown[i,t]
				+ slackQUp[i,t] + slackQDown[i,t] for i=dn_nodes for t=times)
		+ sum(prices_int[e,t]*fp[e,t] for e=dn_ext_edges for t=times)
		)
	end

	if warmstart=="yes"
		#println(" >> WARMSTARTING <<")
		for sb=dn_seg_bid
			set_start_value(pg[sb],sol.pg[sb])
			set_start_value(qg[sb],sol.qg[sb])
		end
		for t=times
			for l=dn_c_index set_start_value(c[l,t],sol.c[l,t]) end
			for l=dn_s_index set_start_value(s[l,t],sol.s[l,t]) end
			for i=dn_nodes
				set_start_value(theta[i,t],sol.theta[i,t])
				set_start_value(slackPUp[i,t],0.0)
				set_start_value(slackPDown[i,t],0.0)
				set_start_value(slackQUp[i,t],0.0)
				set_start_value(slackQDown[i,t],0.0)
			end
			for l=dn_edges
				set_start_value(fp[l,t],sol.fp[l,t])
				set_start_value(fpTo[l,t],sol.fpTo[l,t])
				set_start_value(fq[l,t],sol.fq[l,t])
				set_start_value(fqTo[l,t],sol.fqTo[l,t])
			end
			for l=dn_ext_edges
				set_start_value(fp[l,t],sol.fp[l,t])
				set_start_value(fpTo[l,t],sol.fpTo[l,t])
			end
			for i=dn_ext_nodes set_start_value(theta[i,t],sol.theta[i,t]) end
		end
	end

	optimize!(m)
	#@assert primal_status(m)==MOI.FEASIBLE_POINT
	if !(primal_status(m)==MOI.FEASIBLE_POINT || primal_status(m)==MOI.NEARLY_FEASIBLE_POINT)
		set_optimizer_attribute(m, "warm_start_init_point", "no")
		optimize!(m)
	end
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
	slackPUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPUp[i,t])
	for i=dn_nodes for t=times)
	slackPDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPDown[i,t])
	for i=dn_nodes for t=times)
	slackQUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackQUp[i,t])
	for i=intersect(dn_nodes,dist_nodes) for t=times)
	slackQDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackQDown[i,t])
	for i=intersect(dn_nodes,dist_nodes) for t=times)
	MaxSlack=(maximum(abs,union(value.(slackPUp), value.(slackQUp),value.(slackPDown), value.(slackQDown))))
	#println("Max Slack: $MaxSlack")
	pfsol=PFsol(MaxSlack>1e-5,pgHist,qgHist,cHist,sHist,thetaHist,
	fpHist,fqHist,fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
	dual_sol=Dualsol(priceP,priceQ)
	return (pfsol,dual_sol,solve_time(m))
end

function solve_DN_relax(dn::String,sol_TN,binsol=0,
						prices_int::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}())
        dn_nodes=DNs[dn]
        dn_edges=edges_DN[dn]
        dn_c_index=c_index_DN[dn]
        dn_s_index=s_index_DN[dn]
        dn_seg_bid=seg_bid_DN[dn]
        dn_dist_seg_bid=dist_seg_bid_DN[dn]
        dn_bid=bid_DN[dn]
        dn_qt_bid=qt_bid_DN[dn]
        dn_ex_qt_bid_id=ex_qt_bid_id_DN[dn]
        dn_ext_nodes=ext_nodes_DN[dn]
        dn_ext_edges=ext_edges_DN[dn]
        dn_RampIndexMinIncr=RampIndexMinIncr_DN[dn]
        dn_RampIndexMaxIncr=RampIndexMaxIncr_DN[dn]
        dn_MinDurIndex=MinDurIndex_DN[dn]
        m=Model(optimizer_with_attributes(Mosek.Optimizer,
       "MSK_IPAR_PRESOLVE_USE" => 1,
       "MSK_IPAR_LOG" => 0,
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
        @variable(m, 0 <= c[dn_c_index,times])
        @variable(m, s[dn_s_index,times])
        @variable(m, pg[dn_seg_bid])
        @expression(m, net_pg_bid[bd=dn_bid], sum(pg[sb] for sb=bid_2_seg_bids[bd]))
        @variable(m, qg[dn_seg_bid])
        @expression(m, net_qg_bid[bd=dn_bid], sum(qg[sb] for sb=bid_2_seg_bids[bd]))
        @variable(m, -pi <= theta[union(dn_nodes,dn_ext_nodes),times] <= pi)

        @variable(m, slackPUp[i=dn_nodes,t=times] >= 0)
        @variable(m, slackPDown[i=dn_nodes,t=times] >= 0)
        @variable(m, slackQUp[i=dn_nodes,t=times] >= 0)
        @variable(m, slackQDown[i=dn_nodes,t=times] >= 0)

        ### BINARY VARIABLES ###
		if binsol==0 #&& isempty(prices_int)
		   	@variable(m, sa[dn_seg_bid], Bin)
 		  	@variable(m, qa[dn_bid], Bin)
 		  	@variable(m, qta[dn_qt_bid], Bin)
 		  	@variable(m, alpha[[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]], Bin)
 		  	@variable(m, omega[[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]], Bin)
	   end
	   #if !isempty(prices_int)
		#   @variable(m, 0 <= relax_sa[dn_seg_bid] <= 1)
	  	#end

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
		if binsol==0 #&& isempty(prices_int)
        	@constraint(m, PDown[sb=dn_seg_bid], -pg[sb] <= -sa[sb]*p_min[sb])
        	@constraint(m, PUp[sb=dn_seg_bid], pg[sb] <= sa[sb]*p_max[sb])
		# elseif !isempty(prices_int)
		# 	@constraint(m, PDown[sb=dn_seg_bid], -pg[sb] <= -relax_sa[sb]*p_min[sb])
        # 	@constraint(m, PUp[sb=dn_seg_bid], pg[sb] <= relax_sa[sb]*p_max[sb])
		else
			@constraint(m, PDown[sb=dn_seg_bid], -pg[sb] <= -binsol.sa[sb]*p_min[sb])
			@constraint(m, PUp[sb=dn_seg_bid], pg[sb] <= binsol.sa[sb]*p_max[sb])
		end

        @constraint(m, PowerBalanceP[i=dn_nodes,t=times],
        sum(sum(pg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
                for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
        + net_injP[i,t]
        + slackPUp[i,t]
        ==
        + sum(fpTo[e,t] for e=n_to[i])
        + sum(fp[e,t] for e=n_fr[i])
        + slackPDown[i,t]
        )

        @constraint(m, PowerBalanceQ[i=dn_nodes,t=times],
        sum(sum(qg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
                for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
        + net_injQ[i,t]
        + slackQUp[i,t]
        ==
        + sum(fqTo[e,t] for e=dist_n_to[i])
        + sum(fq[e,t] for e=dist_n_fr[i])
        + slackQDown[i,t]
        )

        @constraint(m, VoltageLimitsDown[i=intersect(dn_nodes,v_limit),t=times],
        -c[(i,i),t] <= -vsq_min[i])

        @constraint(m, VoltageLimitsUp[i=intersect(dn_nodes,v_limit),t=times],
        c[(i,i),t] <= vsq_max[i])

        @constraint(m, LinePowerLimit1[e=intersect(line_limit,dn_edges),t=times],
        [rate_a[e],fp[e,t],fq[e,t]] in SecondOrderCone())

        @constraint(m, LinePowerLimit2[e=intersect(line_limit,dn_edges),t=times],
        [rate_a[e],fpTo[e,t],fqTo[e,t]] in SecondOrderCone())

        mod_s_ind=[(e_fr[e],e_to[e]) for e=dn_edges if R_pu[e]!=0.0]
       @constraint(m, SOCPconstraint[e=mod_s_ind,t=times],
       [0.5*c[(e[1],e[1]),t],c[(e[2],e[2]),t],c[e,t],s[e,t]] in RotatedSecondOrderCone())

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
       [get(QPDiscMax,(BQpd[2],BQpd[3],BQpd[5]),0),
       net_pg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])],
       net_qg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]]
       in SecondOrderCone()
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
        @constraint(m, DefFpX0[e=dn_ext_edges, t=times],
        fp[e,t] + fpTo[e,t] == 0)

        ### BINARY CONSTRAINTS ###
		if binsol==0 #&& isempty(prices_int)
       @constraint(m, ExQtBidsConstr[i=dn_ex_qt_bid_id],
       sum(qta[ExQt] for ExQt=ex_qt_bid[i]) <= 1)

       @constraint(m, QtaDefLow[qt=dn_qt_bid],
       qta[qt] <= sum(qa[bd] for bd=bid if bd[2]==qt))

       @constraint(m, QtaDefHigh[bd=dn_bid],
       qa[bd] <= qta[bd[2]])

       @constraint(m, DefQa[sb=dn_seg_bid],
       sa[sb] <= qa[(sb[1],sb[2],sb[3],sb[5])])

       @constraint(m, AlphaOmegaDef[B=[aol for aol=AlphaOmegaLate if aol[1] in dn_nodes]],
       qa[B]
       - sum(qa[bd] for bd=bid if bd[1]==B[1] && bd[1]==B[2] && bd[4]==B[4]-1)
       - alpha[B]
       + omega[B]
       == 0)

       @constraint(m, MutexAlphaOmega[B=[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]],
       alpha[B] + omega[B] <= 1)

       @constraint(m, AlphaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in dn_nodes]],
       qa[B] - alpha[B] == 0)

       @constraint(m, OmegaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in dn_nodes]],
       omega[B] == 0)

       @constraint(m, MinDurationConstr[mdp=dn_MinDurIndex],
       + sum(qa[bd] for bd=bid if bd[2]==mdp[1] && bd[3]==mdp[3] && bd[4]==mdp[5])
       - sum(alpha[bd] for bd=bid if bd[2]==mdp[1] && bd[3]==mdp[2] && bd[4]==mdp[4])
       >= 0)

       @constraint(m, NoNewAct[B=[nnai for nnai=NoNewActIndex if nnai[1] in dn_nodes]],
       alpha[B] == 0)
   		end

		@constraint(m,fixTheta[i=dn_ext_nodes,t=times], theta[i,t]-sol_TN.theta[i,t]==0)

        # ### OBJECTIVE ###
		if isempty(prices_int)
        @objective(m, Min,
        + sum(cost[sb][1]*pg[sb] for sb=dn_seg_bid)
        + cost_ls*sum(slackPUp[i,t] + slackPDown[i,t]
                        + slackQUp[i,t] + slackQDown[i,t] for i=dn_nodes for t=times)
        )
        ### Flow at the interconnection is fixed to a given value ###
       @constraint(m,fpFix[e=dn_ext_edges,t=times], fp[e,t]-sol_TN.fp[e,t] == 0)
   		else
		@objective(m, Min,
	     + sum(cost[sb][1]*pg[sb] for sb=dn_seg_bid)
	     + cost_ls*sum(slackPUp[i,t] + slackPDown[i,t]
	                        + slackQUp[i,t] + slackQDown[i,t] for i=dn_nodes for t=times)
		+ sum(prices_int[e,t]*fp[e,t] for e=dn_ext_edges for t=times)
	     )
		 @constraint(m, FlowLimUp[e=dn_ext_edges,t=times],
	 	fp[e,t] <= rate_a[e])
	 	@constraint(m, FlowLimDown[e=dn_ext_edges,t=times],
	 	-fp[e,t] <= rate_a[e])
	 	end

        optimize!(m)
		#println(primal_status(m))
		#println(termination_status(m))
        #@assert primal_status(m)==MOI.FEASIBLE_POINT
		if primal_status(m)!=MOI.FEASIBLE_POINT
			socp_sol=PFsol(false,0,0,0,0,0,0,0,0,0,0,0,0,0)
			return (socp_sol,0,binsol,solve_time(m))
		end
	## only build a partial solution for this DN ##
	pgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(pg[j]) for j=dn_seg_bid)
	qgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(qg[j]) for j=dn_dist_seg_bid)
	cHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((i,t) => value(c[i,t]) for i=dn_c_index for t=times)
	sHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((l,t) => value(s[l,t]) for l=dn_s_index for t=times)
    theta_trHist=Dict{Tuple{Int64,Int64},Float64}((i,t) => value(theta[i,t]) for i=union(dn_nodes,dn_ext_nodes) for t=times)
	fpHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fp[l,t]) for l=dn_edges for t=times)
	fqHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fq[l,t]) for l=dn_edges for t=times)
	fpToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fpTo[l,t]) for l=dn_edges for t=times)
	fqToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fqTo[l,t]) for l=dn_edges for t=times)
	slackPUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPUp[i,t])
	for i=dn_nodes for t=times)
	slackPDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPDown[i,t])
	for i=dn_nodes for t=times)
	slackQUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackQUp[i,t])
	for i=intersect(dn_nodes,dist_nodes) for t=times)
	slackQDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackQDown[i,t])
	for i=intersect(dn_nodes,dist_nodes) for t=times)
	MaxSlack=(maximum(abs,union(value.(slackPUp), value.(slackQUp),value.(slackPDown), value.(slackQDown))))
	#println("Max Slack: $MaxSlack")
        socp_sol=PFsol(false,pgHist,qgHist,cHist,sHist,theta_trHist,fpHist,fqHist,
	fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
		# if 10060114 in dn_ext_edges
		# println(socp_sol.theta[e_fr[10060114],1])
		# println(socp_sol.theta[e_to[10060114],1])
		# end
        socp_sol=heurTheta2(socp_sol,dn_nodes)
		# if 10060114 in dn_ext_edges
		# println(socp_sol.theta[e_fr[10060114],1])
		# println(socp_sol.theta[e_to[10060114],1])
		# end
        gaps=ComputeGaps(socp_sol,dn_s_index)
        if abs(gaps["SOCPGap"]) < 1e-6 && abs(gaps["TanGap"]) < 1e-6
                socp_sol.isPFsol=true
        end
		if binsol==0
        dn_bin_sol=Binsol(
		Dict{Any,Any}(sb=>value(sa[sb]) for sb=dn_seg_bid),
		Dict{Any,Any}(bd=>value(qa[bd]) for bd=dn_bid),
		Dict{Any,Any}(qt=>value(qta[qt]) for qt=dn_qt_bid),
		Dict{Any,Any}(sb=>value(alpha[sb]) for sb=AlphaOmegaIndex if sb[1] in dn_nodes),
		Dict{Any,Any}(sb=>value(omega[sb]) for sb=AlphaOmegaIndex if sb[1] in dn_nodes)
			)
			return (socp_sol,0,dn_bin_sol,solve_time(m))
		end
		priceP=Dict{Tuple{Int64,Int64},Float64}((i,t) => dual(PowerBalanceP[i,t])
		for i=dn_nodes for t=times)
		priceQ=Dict{Tuple{Int64,Int64},Float64}((i,t) => dual(PowerBalanceQ[i,t])
		for i=intersect(dn_nodes,dist_nodes) for t=times)
		dual_sol=Dualsol(priceP,priceQ)
		return (socp_sol,dual_sol,binsol,solve_time(m))
end

function solveDCAC_bin_decomposed(sol::PFsol,
								  binsol,
								  prices_int::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}(),
								  warmstart::String="yes")::Tuple{PFsol,Dualsol,Dict{String,Float64}}
	solve_times=Dict{String,Float64}()
	println("Starting decomposition...")
	(sol_TN,dualsol_TN,new_prices_int,solve_times["TN"])=solveDC_TN(sol,binsol,prices_int)
	#(sol_TN,dualsol_TN,new_prices_int,solve_times["TN"])=solveDC_TN(sol,0,prices_int)
	println("> SOLVED Transmission network problem.")
	sol_DN=Dict{String,Any}()
	dualsol_DN=Dict{String,Any}()
	println("Starting $(length(root_dist_node)) Distribution networks problems.")
	for (num,dn) in enumerate(keys(DNs))
		(sol_DN[dn],dualsol_DN[dn],tmp_bin_sol,solve_times[dn])=solve_DN_relax(dn,sol,binsol,prices_int)
		if isempty(prices_int) && sol_DN[dn].isPFsol==false
			(sol_DN[dn],dualsol_DN[dn],tmp_st)=solveAC_subnet(dn,sol,binsol,prices_int,warmstart)
			solve_times[dn]+=tmp_st
		end
		#(sol_DN[dn],dualsol_DN[dn],solve_times[dn])=solveAC_subnet(dn,sol,binsol,prices_int,warmstart)
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

function solve_pert_interface(relax_sol,disturb=0.1)
	println("Starting decomposition...")
	sol_TN=solveDC_no_coordination(relax_sol,disturb)
	println("> SOLVED Transmission network problem.")
	sol_DN=Dict{String,Any}()
	dualsol_DN=Dict{String,Any}()
	println("Starting $(length(root_dist_node)) Distribution networks problems.")
	for (num,dn) in enumerate(keys(DNs))
		(sol_DN[dn],dualsol_DN[dn],tmp_bin_sol,solve_times)=solve_DN_relax(dn,sol_TN)
		if sol_DN[dn].isPFsol==false
			(sol_DN[dn],dualsol_DN[dn],tmp_st)=solveAC_subnet(dn,sol_TN,tmp_bin_sol)
		end
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
	ComputeGaps(pfsol)
	println("SW: $(SW(pfsol))")
	return (pfsol)
end

function solve_no_coordination()
	println("Starting decomposition...")
	sol_TN=solveDC_no_coordination()
	println("> SOLVED Transmission network problem.")
	sol_DN=Dict{String,Any}()
	dualsol_DN=Dict{String,Any}()
	println("Starting $(length(root_dist_node)) Distribution networks problems.")
	for (num,dn) in enumerate(keys(DNs))
		(sol_DN[dn],dualsol_DN[dn],tmp_bin_sol,solve_times)=solve_DN_relax(dn,sol_TN)
		if sol_DN[dn].isPFsol==false
			(sol_DN[dn],dualsol_DN[dn],tmp_st)=solveAC_subnet(dn,sol_TN,tmp_bin_sol)
		end
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
	ComputeGaps(pfsol)
	println("SW: $(SW(pfsol))")
	return (pfsol)
end

function compute_dec_sol(solver::Symbol=:Mosek)::Tuple{PFsol,Dualsol,Dict{String,Float64},Dict{String,Float64}}
	println("\n##############################\n### Step 1: Solve MIDCSOCP ###\n##############################\n")
	(socp_sol,bin_sol)=solveMIDCSOCP()
	println("\n##########################\n### Step 2: Solve SOCP ###\n##########################\n")
	(cont_socp_sol,cont_bin_sol,cont_dual_sol,prices_int)=solveDCSOCP(solver,bin_sol)
	println("\n################################\n### Step 3: Decompose Primal ###\n################################\n")
	(dec_sol,dec_dual_sol_bis,dec_solve_times_pr)=solveDCAC_bin_decomposed(cont_socp_sol,bin_sol)
	# println(dec_sol.theta[e_fr[10060114],1])
	# println(dec_sol.theta[e_to[10060114],1])
	println("\n##############################\n### Step 4: Decompose Dual ###\n##############################\n")
	(dec_sol_bis,dec_dual_sol,dec_solve_times_du)=solveDCAC_bin_decomposed(cont_socp_sol,bin_sol,prices_int)

	total_st_pr=sum(sum(values(dec_solve_times_pr[i])) for i=keys(dec_solve_times_pr))
    max_par_st_pr=maximum(values(dec_solve_times_pr))

	total_st_du=sum(sum(values(dec_solve_times_du[i])) for i=keys(dec_solve_times_du))
    max_par_st_du=maximum(values(dec_solve_times_du))

    println("\n>>> Total primal solve time: $total_st_pr ($max_par_st_pr)\n")
	println("\n>>> Total dual solve time: $total_st_du ($max_par_st_du)\n")
	return (dec_sol,dec_dual_sol,dec_solve_times_pr,dec_solve_times_du)
end

function heurTheta2(SOCPsol::PFsol,set_nodes=dist_nodes)::PFsol
	set_s_index=[s_i for s_i=s_index if s_i[1] in set_nodes && s_i[2] in set_nodes]
	m=Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 0))
	@variable(m,pi >= theta[set_nodes,times] >= -pi)
	for n=root_dist_node if n in set_nodes for t=times fix(theta[n,t],SOCPsol.theta[n,t]; force=true) end end end
	@constraint(m, DefS[l=set_s_index,t=times],
	asin(-SOCPsol.s[l,t]/sqrt(SOCPsol.c[(l[1],l[1]),t]*SOCPsol.c[(l[2],l[2]),t])) == (theta[l[1],t]-theta[l[2],t])
	)
	optimize!(m)
	@assert primal_status(m) == MOI.FEASIBLE_POINT
	newTheta=copy(SOCPsol.theta)
	for i=set_nodes
		for t=times
			newTheta[i,t]=value(theta[i,t])
		end
	end
	return PFsol(false,SOCPsol.pg,SOCPsol.qg,SOCPsol.c,SOCPsol.s,newTheta,
	SOCPsol.fp,SOCPsol.fq,SOCPsol.fpTo,SOCPsol.fqTo,SOCPsol.slackPUp,SOCPsol.slackPDown,SOCPsol.slackQUp,SOCPsol.slackQDown)
end

function ComputeGaps(sol::PFsol,s_ind=s_index)::Dict{String,Float64}
	SOCPGap = mapreduce(maximum,max,abs(sol.c[l,t]^2 + sol.s[l,t]^2
	- sol.c[(l[1],l[1]),t]*sol.c[(l[2],l[2]),t]) for l=s_ind for t=times;init=0)
	TanGap = mapreduce(maximum,max,abs(sol.c[l,t]*sin(sol.theta[l[1],t]-sol.theta[l[2],t])
	+ sol.s[l,t]*cos(sol.theta[l[1],t]-sol.theta[l[2],t]))
	for l=s_ind for t=times; init=0)
	println("SOCPGap: $SOCPGap")
	println("TanGap: $TanGap")
	return Dict("SOCPGap" => SOCPGap, "TanGap" => TanGap)
end

function ComputeGaps_ext(sol::PFsol)
	SOCPGap = findmax(Dict((l,t) => abs(sol.c[l,t]^2 + sol.s[l,t]^2
               - sol.c[(l[1],l[1]),t]*sol.c[(l[2],l[2]),t]) for l=s_index for t=times))
	TanGap = findmax(Dict((l,t) => abs(sol.c[l,t]*sin(sol.theta[l[1],t]-sol.theta[l[2],t])
	+ sol.s[l,t]*cos(sol.theta[l[1],t]-sol.theta[l[2],t])) for l=s_index for t=times))
	println("SOCPGap: $SOCPGap")
	println("TanGap: $TanGap")
	return Dict("SOCPGap" => SOCPGap, "TanGap" => TanGap)
end

function SW(sol::PFsol, set_node=nodes)::Float64
	return (+mapreduce(sum,+,init=0,cost[sb][1]*sol.pg[sb] for sb=seg_bid if sb[1] in set_node)
	+cost_ls*mapreduce(sum,+,init=0,sol.slackPUp[i,t]+sol.slackPDown[i,t] for i=set_node for t=times)
	+cost_ls*mapreduce(sum,+,init=0,sol.slackQUp[i,t]+sol.slackQDown[i,t] for i=set_node for t=times if i in dist_nodes)
	)
end

function SlackInfo(sol::PFsol)::Dict{String,Float64}
	MaxSlack=maximum(union(values(sol.slackPUp),values(sol.slackPDown),values(sol.slackQUp),values(sol.slackQDown)))
	SumSlack=sum(values(sol.slackPUp))+sum(values(sol.slackPDown))+sum(values(sol.slackQUp))+sum(values(sol.slackQDown))
	println("MaxSlack: ", MaxSlack)
	println("SumSlack: ", SumSlack)
	return Dict("MaxSlack"=>MaxSlack, "SumSlack" => SumSlack)
end

function in_DN(node::Int64)::String
	for dn=keys(DNs)
		if node in DNs[dn] return dn end
	end
	return "?"
end

function correct_sol(sol::PFsol,dual_sol::Dualsol,bin_sol::Binsol)
	max_slack=findmax(ComputeGaps_ext(sol))
	corr_st=Dict{String,Float64}()
	while max_slack[1][1] > 1e-6
		dn=in_DN(max_slack[1][2][1][1])
		(sol_DN,dualsol_DN,corr_st[dn])=solveAC_subnet(dn,sol,bin_sol)
		for sb=seg_bid_DN[dn]
			sol.pg[sb]=sol_DN.pg[sb]
			sol.qg[sb]=sol_DN.qg[sb]
		end
		for i=c_index_DN[dn] for t=times sol.c[i,t] = sol_DN.c[i,t] end end
		for i=s_index_DN[dn] for t=times sol.s[i,t] = sol_DN.s[i,t] end end
		for i=DNs[dn] for t=times
			 sol.theta[i,t]=sol_DN.theta[i,t]
			 sol.slackPUp[i,t]=sol_DN.slackPUp[i,t]
	 		 sol.slackPDown[i,t]=sol_DN.slackPDown[i,t]
			 sol.slackQUp[i,t]=sol_DN.slackQUp[i,t]
	 		 sol.slackQDown[i,t]=sol_DN.slackQDown[i,t]
			 #dual_sol.priceP[i,t]=dualsol_DN.priceP[i,t]
			 #dual_sol.priceQ[i,t]=dualsol_DN.priceQ[i,t]
		end end
		for e=edges_DN[dn] for t=times
			sol.fp[e,t]=sol_DN.fp[e,t]
			sol.fpTo[e,t]=sol_DN.fpTo[e,t]
		end end
		for e=edges_DN[dn] for t=times
			sol.fq[e,t]=sol_DN.fq[e,t]
			sol.fqTo[e,t]=sol_DN.fqTo[e,t]
		end end
		println(">Corrected Distribution network problem: $dn.")
		max_slack=findmax(ComputeGaps_ext(sol))
	end
	return (sol,dual_sol,corr_st)
end

function exec_benchmark(dir=0)
	(socp_sol,bin_sol)=(socp_sol,bin_sol)=solveMIDCSOCP()
	(origpfsol,origdualsol)=solveDCAC_bin(socp_sol,bin_sol)
	additional_info(origpfsol,origdualsol)
	if dir!=0 set_ads_tables_to_CSV(dir,origpfsol,origdualsol) end
	return (origpfsol,origdualsol)
end

function new_exec_benchmark(dir=0)
	(socp_sol,bin_sol)=solveMIDCSOCP()
	(cont_socp_sol,cont_bin_sol,cont_dual_sol,socp_prices_int)=solveDCSOCP(:Mosek,bin_sol)
	(origpfsol,origdualsol,corr_st)=correct_sol(cont_socp_sol,cont_dual_sol,bin_sol)
	println("	>>> Total correction solve time: $(sum(values(corr_st))) ($(maximum(values(corr_st))))")
	additional_info(origpfsol,origdualsol)
	if dir!=0 set_ads_tables_to_CSV(dir,origpfsol,origdualsol) end
	return (origpfsol,origdualsol)
end

function exec_relax(dir=0)
	(socp_sol,bin_sol)=(socp_sol,bin_sol)=solveMIDCSOCP()
	(cont_socp_sol,cont_bin_sol,cont_dual_sol,socp_prices_int)=solveDCSOCP(:Mosek,bin_sol)
	additional_info(cont_socp_sol,cont_dual_sol)
	if dir!=0 set_ads_tables_to_CSV(dir,cont_socp_sol,cont_dual_sol) end
	return (cont_socp_sol,cont_dual_sol)
end

function exec_relax_dec(dir=0)
	(dec_sol,dec_dual_sol,ST_primal,ST_dual)=compute_dec_sol()
	additional_info(dec_sol,dec_dual_sol)
	if dir!=0 set_ads_tables_to_CSV(dir,dec_sol,dec_dual_sol) end
	return (dec_sol,dec_dual_sol)
end

function exec_rsf(nb_points::Int64, dir=0, exp_flow::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}())
	(rsf_sol, rsf_binsol,rsf_dual,prices_int,rsf_st)=compute_decentralized_bin(nb_points,exp_flow);
	additional_info(rsf_sol,rsf_dual)
	additional_time_info(rsf_st)
	if dir!=0 set_ads_tables_to_CSV(dir,rsf_sol,rsf_dual,prices_int) end
	return (rsf_sol,rsf_dual)
end

function exec_enhance_rsf(nb_points::Int64, dir=0, exp_flow::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}())
	(rsf_sol, rsf_binsol,rsf_dual, prices_int,rsf_st)=enhance_compute_decentralized_bin(nb_points,exp_flow);
	additional_info(rsf_sol,rsf_dual)
	additional_time_info(rsf_st)
	if dir!=0 set_ads_tables_to_CSV(dir,rsf_sol,rsf_dual,prices_int) end
	return (rsf_sol,rsf_dual)
end

function exec_no_last_step(nb_points::Int64, dir=0, exp_flow::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}())
	(rsf_sol, rsf_binsol,rsf_dual, prices_int,rsf_st)=enhance_compute_no_rsf_last_step(nb_points,exp_flow);
	additional_info(rsf_sol,rsf_dual)
	additional_time_info(rsf_st)
	if dir!=0 set_ads_tables_to_CSV(dir,rsf_sol,rsf_dual,prices_int) end
	return (rsf_sol,rsf_dual)
end

function exec_enhance_rsf_oneill(nb_points::Int64, dir=0, exp_flow::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}())
	(rsf_sol, rsf_binsol,rsf_dual, prices_int,rsf_st)=enhance_compute_decentralized_oneill(nb_points,exp_flow);
	additional_info(rsf_sol,rsf_dual)
	additional_time_info(rsf_st)
	if dir!=0 set_ads_tables_to_CSV(dir,rsf_sol,rsf_dual,prices_int) end
	return (rsf_sol,rsf_dual)
end

function exec_enhance_rsf_exact(nb_points::Int64, dir=0, exp_flow::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}())
	(rsf_sol, rsf_binsol,rsf_dual, prices_int,rsf_st)=enhance_compute_decentralized_exact_bin(nb_points,exp_flow)
	additional_info(rsf_sol,rsf_dual)
	additional_time_info(rsf_st)
	if dir!=0 set_ads_tables_to_CSV(dir,rsf_sol,rsf_dual,prices_int) end
	return (rsf_sol,rsf_dual)
end

function n_proc_st(nb_proc::Int64,solve_times)
	total_st=0.0
	cur_max=0.0
	for i=keys(solve_times)
		if ceil(Int64,length(solve_times[i])/nb_proc) <= 1
			cur_max=findmax(solve_times[i])
			total_st+=cur_max[1]
		else
			copy_st=copy(solve_times[i])
			for n=1:ceil(Int64,length(solve_times[i])/nb_proc)
				cur_max=findmax(copy_st)
				total_st+=cur_max[1]
				delete!(copy_st,cur_max[2])
			end
		end
	end
	return total_st
end

function ConsumerPayment(sol::PFsol, dualsol::Dualsol)::Tuple{Float64,Float64}
	return (tr_ConsumerPayment(sol,dualsol)+dist_ConsumerPayment(sol,dualsol)[1],
			dist_ConsumerPayment(sol,dualsol)[2])
end

function tr_ConsumerPayment(sol::PFsol, dualsol::Dualsol)::Float64
	return -(mapreduce(sum,+,sol.pg[sb]*dualsol.priceP[sb[1],sb[5]] for sb=tr_seg_bid if sol.pg[sb]<=0;init=0.0)
	+ mapreduce(sum,+,dualsol.priceP[i,t]*get(net_injP,(i,t),0) for i=tr_nodes for t=times if get(net_injP,(i,t),0)<=0;init=0.0))
end

function dist_ConsumerPayment(sol::PFsol, dualsol::Dualsol)::Tuple{Float64,Float64}
	return (-mapreduce(sum,+,sol.pg[sb]*dualsol.priceP[sb[1],sb[5]] for sb=dist_seg_bid if sol.pg[sb]<=0;init=0.0)
	- mapreduce(sum,+,dualsol.priceP[i,t]*get(net_injP,(i,t),0) for i=dist_nodes for t=times if get(net_injP,(i,t),0)<=0;init=0.0),
	- mapreduce(sum,+,sol.qg[sb]*dualsol.priceQ[sb[1],sb[5]] for sb=dist_seg_bid if sol.qg[sb]<=0;init=0.0)
	- mapreduce(sum,+,dualsol.priceQ[i,t]*get(net_injQ,(i,t),0) for i=dist_nodes for t=times if get(net_injQ,(i,t),0)<=0;init=0.0))
end

function ProducerRevenue(sol::PFsol, dualsol::Dualsol)::Tuple{Float64,Float64}
	return (tr_ProducerRevenue(sol,dualsol)+dist_ProducerRevenue(sol,dualsol)[1],
			dist_ProducerRevenue(sol,dualsol)[2])
end

function tr_ProducerRevenue(sol::PFsol, dualsol::Dualsol)::Float64
	return (mapreduce(sum,+,sol.pg[sb]*dualsol.priceP[sb[1],sb[5]] for sb=tr_seg_bid if sol.pg[sb]>=0;init=0.0)
	+ mapreduce(sum,+,dualsol.priceP[i,t]*get(net_injP,(i,t),0) for i=tr_nodes for t=times if get(net_injP,(i,t),0)>=0;init=0.0))
end

function dist_ProducerRevenue(sol::PFsol, dualsol::Dualsol)::Tuple{Float64,Float64}
	return (+ mapreduce(sum,+,sol.pg[sb]*dualsol.priceP[sb[1],sb[5]] for sb=dist_seg_bid if sol.pg[sb]>=0;init=0.0)
	+ mapreduce(sum,+,dualsol.priceP[i,t]*get(net_injP,(i,t),0) for i=dist_nodes for t=times if get(net_injP,(i,t),0)>=0;init=0.0),
	+ mapreduce(sum,+,sol.qg[sb]*dualsol.priceQ[sb[1],sb[5]] for sb=dist_seg_bid if sol.qg[sb]>=0;init=0.0)
	+ mapreduce(sum,+,dualsol.priceQ[i,t]*get(net_injQ,(i,t),0) for i=dist_nodes for t=times if get(net_injQ,(i,t),0)>=0;init=0.0))
end

function Merchandising_Surplus(sol::PFsol,dualsol::Dualsol)
	return (sum(ConsumerPayment(sol,dualsol))-sum(ProducerRevenue(sol,dualsol)))
end

function SellerRevenue(sol::PFsol,dualsol::Dualsol)
	return (sum(max(0,sol.pg[sb]*dualsol.priceP[sb[1],sb[5]]) for sb=seg_bid)
	+ sum(max(0,dualsol.priceP[i,t]*net_injP[i,t]) for i=nodes for t=times)
	+ sum(max(0,sol.qg[sb]*dualsol.priceQ[sb[1],sb[5]]) for sb=dist_seg_bid)
	+ sum(max(0,dualsol.priceQ[i,t]*net_injQ[i,t]) for i=dist_nodes for t=times))
end

function additional_info(sol::PFsol,dualsol::Dualsol)
	(loco,operatorsol)=operator_subproblem_relax(sol,dualsol)
	locg=locg_decentralized(sol,dualsol)
	norm_sol=sum((p_min[sb]+p_max[sb])*cost[sb][1]/2 for sb=seg_bid)
	println("\n#################### Additional info ####################\n")
	println(">>> Slack Info: ")
	SlackInfo(sol)
	println(">>> Violations: ")
	ComputeGaps(sol)
	println(">>> LOC = $(max(0.0,sum(values(loco)))+max(0.0,sum(values(locg))))")
	println(">>> PLP = $(SellerRevenue(sol,dualsol))")
	println(">>> SW = $(SW(sol))")
	println(">>> normSW = $norm_sol")
	println(">>> Merchandising Surplus = $(Merchandising_Surplus(sol,dualsol))")
end

function additional_time_info(solve_times)
	println("\n#################### Additional info ####################\n")
	println(">>> Sequential Time: $(sum(sum(values(solve_times[i])) for i=keys(solve_times)))")
	println(">>> Parallel opt Time: $(sum(maximum(values(solve_times[i])) for i=keys(solve_times)))")
	println(">>> Parrallel16 Time: $(n_proc_st(16,solve_times))")
	println(">>> Parrallel32 Time: $(n_proc_st(32,solve_times))")
	println(">>> Parrallel64 Time: $(n_proc_st(64,solve_times))")
end
