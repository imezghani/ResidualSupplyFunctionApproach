function compute_res_fun_bin(dn::String,
                         exp_flow::Dict{Tuple{Int64,Int64},Float64},
                         m::Model=Model())
        dn_ext_edges=ext_edges_DN[dn]
        if isempty(m.obj_dict)
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
         if slack_bus in union(dn_nodes,dn_ext_nodes)
                 for i=slack_bus for t=times fix(theta[i,t],0;force=true) end end
         end
         #for i=dn_ext_nodes for t=times fix(theta[i,t],0;force=true) end end
         #@constraint(m, fixTheta0[i=dn_ext_nodes,t=times],theta[i,t]==0)

         # @variable(m, abs(net_injP[i,t]) >= slackPUp[i=dn_nodes,t=times] >= 0)
         # @variable(m, abs(net_injP[i,t]) >= slackPDown[i=dn_nodes,t=times] >= 0)
         # @variable(m, abs(net_injQ[i,t]) >= slackQUp[i=dn_nodes,t=times] >= 0)
         # @variable(m, abs(net_injQ[i,t]) >= slackQDown[i=dn_nodes,t=times] >= 0)
         @variable(m, slackPUp[i=dn_nodes,t=times] >= 0)
         @variable(m, slackPDown[i=dn_nodes,t=times] >= 0)
         @variable(m, slackQUp[i=dn_nodes,t=times] >= 0)
         @variable(m, slackQDown[i=dn_nodes,t=times] >= 0)
         #@variable(m, slackFUp[dn_ext_edges,times] >= 0)
         #@variable(m, slackFDown[dn_ext_edges,times] >= 0)

         ### BINARY VARIABLES ###
 	@variable(m, 0 <= sa[dn_seg_bid] <= 1)
 	@variable(m, 0 <= qa[dn_bid] <= 1)
 	@variable(m, 0 <= qta[dn_qt_bid] <= 1)
 	@variable(m, 0 <= alpha[[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]] <= 1)
 	@variable(m, 0 <= omega[[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]] <= 1)

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
         @constraint(m, PDown[sb=dn_seg_bid], -pg[sb] <= -sa[sb]*p_min[sb])
         @constraint(m, PUp[sb=dn_seg_bid], pg[sb] <= sa[sb]*p_max[sb])


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
         + sum(net_pg_bid[bd] for bd=dn_bid if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
         - sum(net_pg_bid[bd] for bd=dn_bid if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
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
         X_pu[e]*fp[e,t] #+ slackFUp[e,t] - slackFDown[e,t]
         == (theta[e_fr[e],t]-theta[e_to[e],t]))
         @constraint(m, DefFpX0[e=dn_ext_edges, t=times],
         fp[e,t] + fpTo[e,t] == 0)

         ### BINARY CONSTRAINTS ###
 	@constraint(m, ExQtBidsConstr[i=dn_ex_qt_bid_id],
 	sum(qta[ExQt] for ExQt=ex_qt_bid[i]) <= 1)

        ## this constraint is not really necessary,
        ## I also noticed that it might lead to weird numerical behavior
        ## in terms of feasiblity and dual value returned
 	# @constraint(m, QtaDefLow[qt=dn_qt_bid],
 	# qta[qt] <= sum(qa[bd] for bd=dn_bid if bd[2]==qt))

 	@constraint(m, QtaDefHigh[bd=dn_bid],
 	qa[bd] <= qta[bd[2]])

 	@constraint(m, DefQa[sb=dn_seg_bid],
 	sa[sb] <= qa[(sb[1],sb[2],sb[3],sb[5])])

 	@constraint(m, AlphaOmegaDef[B=[aol for aol=AlphaOmegaLate if aol[1] in dn_nodes]],
 	qa[B]
 	- sum(qa[bd] for bd=dn_bid if bd[1]==B[1] && bd[1]==B[2] && bd[4]==B[4]-1)
 	- alpha[B]
 	+ omega[B]
 	== 0)

 	@constraint(m, MutexAlphaOmega[B=[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]],
 	alpha[B] + omega[B] <= 1)

 	@constraint(m, AlphaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in dn_nodes]],
 	qa[B] - alpha[B] <= 0)

 	@constraint(m, OmegaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in dn_nodes]],
 	omega[B] == 0)

 	@constraint(m, MinDurationConstr[mdp=dn_MinDurIndex],
 	+ sum(qa[bd] for bd=dn_bid if bd[2]==mdp[1] && bd[3]==mdp[3] && bd[4]==mdp[5])
  	- sum(alpha[bd] for bd=dn_bid if bd[2]==mdp[1] && bd[3]==mdp[2] && bd[4]==mdp[4])
 	>= 0)
        #println(MinDurationConstr)

 	@constraint(m, NoNewAct[B=[nnai for nnai=NoNewActIndex if nnai[1] in dn_nodes]],
 	alpha[B] == 0)

         # ### OBJECTIVE ###
         @objective(m, Min,
         + sum(cost[sb][1]*pg[sb] for sb=dn_seg_bid)
         + cost_ls*(sum(slackPUp[i,t] + slackPDown[i,t]
                         + slackQUp[i,t] + slackQDown[i,t] for i=dn_nodes for t=times)
                #+sum(slackFUp[e,t] + slackFDown[e,t] for e=dn_ext_edges for t=times)
                        )
         )

         ### Flow at the interconnection is fixed to a given value ###
        @constraint(m,fpFix[e=dn_ext_edges,t=times], fp[e,t]-exp_flow[e,t] == 0)
        else
                for e=dn_ext_edges for t=times
                        set_normalized_rhs(m.obj_dict[:fpFix][e,t],exp_flow[e,t])
                end end
        end
        #@constraint(m, fpFixUp[e=dn_ext_edges,t=times], sign(exp_flow[e,t])*fp[e,t] >= (1-sign(exp_flow[e,t])*0.001)*exp_flow[e,t])
        #@constraint(m, fpFixDown[e=dn_ext_edges,t=times], sign(exp_flow[e,t])*fp[e,t] <= (1+sign(exp_flow[e,t])*0.001)*exp_flow[e,t])

         optimize!(m)
         #println(primal_status(m))
         # println(dual_status(m))
         # println(raw_status(m))
         #println(termination_status(m),"\n")
         # if termination_status(m)!=MOI.OPTIMAL#primal_status(m)!=MOI.FEASIBLE_POINT
         #         set_optimizer_attribute(m,"MSK_IPAR_PRESOLVE_USE",0)
         #         optimize!(m)
         #         set_optimizer_attribute(m,"MSK_IPAR_PRESOLVE_USE",1)
         # end
         # println(primal_status(m))
         # println(termination_status(m),"\n")
         qty_price_interface=Dict{Tuple{Int64,Int64},Tuple{Float64,Float64}}()

        #qty_price_interface=Dict((e,t) => (termination_status(m),dual(m.obj_dict[:fpFix][e,t]),exp_flow[e,t]) for e=dn_ext_edges for t=times)
        #if primal_status(m)!=MOI.FEASIBLE_POINT
        if termination_status(m)!=MOI.OPTIMAL
                for e=dn_ext_edges for t=times
                        #qty_price_interface[e,t]=(MOI.SLOW_PROGRESS,10*cost_ls,exp_flow[e,t])
                        qty_price_interface[e,t]=(10*cost_ls,exp_flow[e,t])
                end end
        else
                for e=dn_ext_edges for t=times
                        qty_price_interface[e,t]=(dual(m.obj_dict[:fpFix][e,t]),exp_flow[e,t])
                end end
        end
         return (qty_price_interface,m,solve_time(m))
end

function compute_res_fun_bin_oneill(dn::String,
                         exp_flow::Dict{Tuple{Int64,Int64},Float64},
                         m::Model=Model(),
                         relax_m::Model=Model())
        dn_ext_edges=ext_edges_DN[dn]
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
        dn_RampIndexMinIncr=RampIndexMinIncr_DN[dn]
        dn_RampIndexMaxIncr=RampIndexMaxIncr_DN[dn]
        dn_MinDurIndex=MinDurIndex_DN[dn]
        if isempty(m.obj_dict)

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
         if slack_bus in union(dn_nodes,dn_ext_nodes)
                 for i=slack_bus for t=times fix(theta[i,t],0;force=true) end end
         end
         #for i=dn_ext_nodes for t=times fix(theta[i,t],0;force=true) end end
         #@constraint(m, fixTheta0[i=dn_ext_nodes,t=times],theta[i,t]==0)

         # @variable(m, abs(net_injP[i,t]) >= slackPUp[i=dn_nodes,t=times] >= 0)
         # @variable(m, abs(net_injP[i,t]) >= slackPDown[i=dn_nodes,t=times] >= 0)
         # @variable(m, abs(net_injQ[i,t]) >= slackQUp[i=dn_nodes,t=times] >= 0)
         # @variable(m, abs(net_injQ[i,t]) >= slackQDown[i=dn_nodes,t=times] >= 0)
         @variable(m, slackPUp[i=dn_nodes,t=times] >= 0)
         @variable(m, slackPDown[i=dn_nodes,t=times] >= 0)
         @variable(m, slackQUp[i=dn_nodes,t=times] >= 0)
         @variable(m, slackQDown[i=dn_nodes,t=times] >= 0)
         #@variable(m, slackFUp[dn_ext_edges,times] >= 0)
         #@variable(m, slackFDown[dn_ext_edges,times] >= 0)

         ### BINARY VARIABLES ###
 	@variable(m, sa[dn_seg_bid], Bin)
 	@variable(m, qa[dn_bid], Bin)
 	@variable(m, qta[dn_qt_bid], Bin)
 	@variable(m, alpha[[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]], Bin)
 	@variable(m, omega[[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]], Bin)

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
         @constraint(m, PDown[sb=dn_seg_bid], -pg[sb] <= -sa[sb]*p_min[sb])
         @constraint(m, PUp[sb=dn_seg_bid], pg[sb] <= sa[sb]*p_max[sb])


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
         [1rate_a[e],1fp[e,t],1fq[e,t]] in SecondOrderCone())

         @constraint(m, LinePowerLimit2[e=intersect(line_limit,dn_edges),t=times],
         [1rate_a[e],1fpTo[e,t],1fqTo[e,t]] in SecondOrderCone())

         mod_s_ind=[(e_fr[e],e_to[e]) for e=dn_edges if R_pu[e]!=0.0]
        @constraint(m, SOCPconstraint[e=mod_s_ind,t=times],
        [0.5*c[(e[1],e[1]),t],1c[(e[2],e[2]),t],1c[e,t],1s[e,t]] in RotatedSecondOrderCone())

         @constraint(m, RampMinIncr[rcI=dn_RampIndexMinIncr],
         + sum(net_pg_bid[bd] for bd=dn_bid if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
         - sum(net_pg_bid[bd] for bd=dn_bid if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
         - get(RCRealPowerIncr,(rcI[1],rcI[2],rcI[4]),0)
         >= 0)

         @constraint(m, RampMaxIncr[rcI=dn_RampIndexMaxIncr],
         + sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
         + get(RCRealPowerIncr,(rcI[1],rcI[2],rcI[4]),0)
         - sum(net_pg_bid[bd] for bd=bid if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
         >= 0)

         @constraint(m, DiscQtoP[BQpd=[qpd for qpd=QPDiscConstrIndex if qpd[1] in dn_nodes]],
        [1get(QPDiscMax,(BQpd[2],BQpd[3],BQpd[5]),0),
        1net_pg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])],
        1net_qg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]]
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
         X_pu[e]*fp[e,t] #+ slackFUp[e,t] - slackFDown[e,t]
         == (theta[e_fr[e],t]-theta[e_to[e],t]))
         @constraint(m, DefFpX0[e=dn_ext_edges, t=times],
         fp[e,t] + fpTo[e,t] == 0)

         ### BINARY CONSTRAINTS ###
 	@constraint(m, ExQtBidsConstr[i=dn_ex_qt_bid_id],
 	sum(qta[ExQt] for ExQt=ex_qt_bid[i]) <= 1)

        ## this constraint is not really necessary,
        ## I also noticed that it might lead to weird numerical behavior
        ## in terms of feasiblity and dual value returned
 	# @constraint(m, QtaDefLow[qt=dn_qt_bid],
 	# qta[qt] <= sum(qa[bd] for bd=dn_bid if bd[2]==qt))

 	@constraint(m, QtaDefHigh[bd=dn_bid],
 	qa[bd] <= qta[bd[2]])

 	@constraint(m, DefQa[sb=dn_seg_bid],
 	sa[sb] <= qa[(sb[1],sb[2],sb[3],sb[5])])

 	@constraint(m, AlphaOmegaDef[B=[aol for aol=AlphaOmegaLate if aol[1] in dn_nodes]],
 	qa[B]
 	- sum(qa[bd] for bd=dn_bid if bd[1]==B[1] && bd[1]==B[2] && bd[4]==B[4]-1)
 	- alpha[B]
 	+ omega[B]
 	== 0)

 	@constraint(m, MutexAlphaOmega[B=[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]],
 	alpha[B] + omega[B] <= 1)

 	@constraint(m, AlphaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in dn_nodes]],
 	qa[B] - alpha[B] <= 0)

 	@constraint(m, OmegaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in dn_nodes]],
 	omega[B] == 0)

 	@constraint(m, MinDurationConstr[mdp=dn_MinDurIndex],
 	+ sum(qa[bd] for bd=dn_bid if bd[2]==mdp[1] && bd[3]==mdp[3] && bd[4]==mdp[5])
  	- sum(alpha[bd] for bd=dn_bid if bd[2]==mdp[1] && bd[3]==mdp[2] && bd[4]==mdp[4])
 	>= 0)
        #println(MinDurationConstr)

 	@constraint(m, NoNewAct[B=[nnai for nnai=NoNewActIndex if nnai[1] in dn_nodes]],
 	alpha[B] == 0)

         # ### OBJECTIVE ###
         @objective(m, Min,
         + sum(cost[sb][1]*pg[sb] for sb=dn_seg_bid)
         + cost_ls*(sum(slackPUp[i,t] + slackPDown[i,t]
                         + slackQUp[i,t] + slackQDown[i,t] for i=dn_nodes for t=times)
                #+sum(slackFUp[e,t] + slackFDown[e,t] for e=dn_ext_edges for t=times)
                        )
         )

         ### Flow at the interconnection is fixed to a given value ###
        @constraint(m,fpFix[e=dn_ext_edges,t=times], fp[e,t]-exp_flow[e,t] == 0)
        relax_m=copy(m)
        unset_binary.(relax_m.obj_dict[:sa])
        unset_binary.(relax_m.obj_dict[:qa])
        unset_binary.(relax_m.obj_dict[:qta])
        unset_binary.(relax_m.obj_dict[:alpha])
        unset_binary.(relax_m.obj_dict[:omega])
        set_optimizer(relax_m,Mosek.Optimizer)
        set_optimizer_attributes(relax_m,"MSK_IPAR_PRESOLVE_USE" => 1,
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
        "MSK_DPAR_INTPNT_CO_TOL_NEAR_REL" => 1)
        else
                for e=dn_ext_edges for t=times
                        set_normalized_rhs(m.obj_dict[:fpFix][e,t],exp_flow[e,t])
                        set_normalized_rhs(relax_m.obj_dict[:fpFix][e,t],exp_flow[e,t])
                end end
        end
         optimize!(m)

         qty_price_interface=Dict{Tuple{Int64,Int64},Tuple{Float64,Float64}}()
         if termination_status(m)!=MOI.OPTIMAL
                 for e=dn_ext_edges for t=times
                         qty_price_interface[e,t]=(10*cost_ls,exp_flow[e,t])
                 end end
                 return (qty_price_interface,m,solve_time(m))
         end


         FixSa=@constraint(relax_m, [sb=dn_seg_bid], relax_m.obj_dict[:sa][sb]==round(value(m.obj_dict[:sa][sb])))
         FixQa=@constraint(relax_m, [bd=dn_bid], relax_m.obj_dict[:qa][bd]==round(value(m.obj_dict[:qa][bd])))
         FixQta=@constraint(relax_m, [bd=dn_qt_bid], relax_m.obj_dict[:qta][bd]==round(value(m.obj_dict[:qta][bd])))
         FixAlpha=@constraint(relax_m, [sb=[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]], relax_m.obj_dict[:alpha][sb]==round(value(m.obj_dict[:alpha][sb])))
         FixOmega=@constraint(relax_m, [sb=[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]], relax_m.obj_dict[:omega][sb]==round(value(m.obj_dict[:omega][sb])))
        optimize!(relax_m)
        delete.(relax_m,FixSa)
        delete.(relax_m,FixQa)
        delete.(relax_m,FixQta)
        delete.(relax_m,FixAlpha)
        delete.(relax_m,FixOmega)
        if termination_status(relax_m)!=MOI.OPTIMAL || primal_status(relax_m)!=MOI.FEASIBLE_POINT
                for e=dn_ext_edges for t=times
                        qty_price_interface[e,t]=(10*cost_ls,exp_flow[e,t])
                end end
        else
                for e=dn_ext_edges for t=times
                        qty_price_interface[e,t]=(dual(relax_m.obj_dict[:fpFix][e,t]),exp_flow[e,t])
                end end
        end
         return (qty_price_interface,m,relax_m,solve_time(m))
end

function compute_res_fun_bin_t(dn::String,
                         exp_flow::Dict{Tuple{Int64,Int64},Float64},
                         timing::Int64,
                         m::Model=Model())
        dn_ext_edges=ext_edges_DN[dn]
        if isempty(m.obj_dict)
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
         if slack_bus in union(dn_nodes,dn_ext_nodes)
                 for i=slack_bus for t=times fix(theta[i,t],0;force=true) end end
         end
         #for i=dn_ext_nodes for t=times fix(theta[i,t],0;force=true) end end
         #@constraint(m, fixTheta0[i=dn_ext_nodes,t=times],theta[i,t]==0)

         # @variable(m, abs(net_injP[i,t]) >= slackPUp[i=dn_nodes,t=times] >= 0)
         # @variable(m, abs(net_injP[i,t]) >= slackPDown[i=dn_nodes,t=times] >= 0)
         # @variable(m, abs(net_injQ[i,t]) >= slackQUp[i=dn_nodes,t=times] >= 0)
         # @variable(m, abs(net_injQ[i,t]) >= slackQDown[i=dn_nodes,t=times] >= 0)
         @variable(m, slackPUp[i=dn_nodes,t=times] >= 0)
         @variable(m, slackPDown[i=dn_nodes,t=times] >= 0)
         @variable(m, slackQUp[i=dn_nodes,t=times] >= 0)
         @variable(m, slackQDown[i=dn_nodes,t=times] >= 0)
         #@variable(m, slackFUp[dn_ext_edges,times] >= 0)
         #@variable(m, slackFDown[dn_ext_edges,times] >= 0)

         ### BINARY VARIABLES ###
 	@variable(m, 0 <= sa[dn_seg_bid] <= 1)
 	@variable(m, 0 <= qa[dn_bid] <= 1)
 	@variable(m, 0 <= qta[dn_qt_bid] <= 1)
 	@variable(m, 0 <= alpha[[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]] <= 1)
 	@variable(m, 0 <= omega[[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]] <= 1)

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
         @constraint(m, PDown[sb=dn_seg_bid], -pg[sb] <= -sa[sb]*p_min[sb])
         @constraint(m, PUp[sb=dn_seg_bid], pg[sb] <= sa[sb]*p_max[sb])


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
         + sum(net_pg_bid[bd] for bd=dn_bid if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
         - sum(net_pg_bid[bd] for bd=dn_bid if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
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
         X_pu[e]*fp[e,t] #+ slackFUp[e,t] - slackFDown[e,t]
         == (theta[e_fr[e],t]-theta[e_to[e],t]))
         @constraint(m, DefFpX0[e=dn_ext_edges, t=times],
         fp[e,t] + fpTo[e,t] == 0)

         ### BINARY CONSTRAINTS ###
 	@constraint(m, ExQtBidsConstr[i=dn_ex_qt_bid_id],
 	sum(qta[ExQt] for ExQt=ex_qt_bid[i]) <= 1)

        ## this constraint is not really necessary,
        ## I also noticed that it might lead to weird numerical behavior
        ## in terms of feasiblity and dual value returned
 	# @constraint(m, QtaDefLow[qt=dn_qt_bid],
 	# qta[qt] <= sum(qa[bd] for bd=dn_bid if bd[2]==qt))

 	@constraint(m, QtaDefHigh[bd=dn_bid],
 	qa[bd] <= qta[bd[2]])

 	@constraint(m, DefQa[sb=dn_seg_bid],
 	sa[sb] <= qa[(sb[1],sb[2],sb[3],sb[5])])

 	@constraint(m, AlphaOmegaDef[B=[aol for aol=AlphaOmegaLate if aol[1] in dn_nodes]],
 	qa[B]
 	- sum(qa[bd] for bd=dn_bid if bd[1]==B[1] && bd[1]==B[2] && bd[4]==B[4]-1)
 	- alpha[B]
 	+ omega[B]
 	== 0)

 	@constraint(m, MutexAlphaOmega[B=[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]],
 	alpha[B] + omega[B] <= 1)

 	@constraint(m, AlphaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in dn_nodes]],
 	qa[B] - alpha[B] <= 0)

 	@constraint(m, OmegaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in dn_nodes]],
 	omega[B] == 0)

 	@constraint(m, MinDurationConstr[mdp=dn_MinDurIndex],
 	+ sum(qa[bd] for bd=dn_bid if bd[2]==mdp[1] && bd[3]==mdp[3] && bd[4]==mdp[5])
  	- sum(alpha[bd] for bd=dn_bid if bd[2]==mdp[1] && bd[3]==mdp[2] && bd[4]==mdp[4])
 	>= 0)
        #println(MinDurationConstr)

 	@constraint(m, NoNewAct[B=[nnai for nnai=NoNewActIndex if nnai[1] in dn_nodes]],
 	alpha[B] == 0)

         # ### OBJECTIVE ###
         @objective(m, Min,
         + sum(cost[sb][1]*pg[sb] for sb=dn_seg_bid)
         + cost_ls*(sum(slackPUp[i,t] + slackPDown[i,t]
                         + slackQUp[i,t] + slackQDown[i,t] for i=dn_nodes for t=times)
                #+sum(slackFUp[e,t] + slackFDown[e,t] for e=dn_ext_edges for t=times)
                        )
         )

         ### Flow at the interconnection is fixed to a given value ###
         @constraint(m, TrPowerFlowLimitUp[e=intersect(line_limit,dn_ext_edges),t=times],
 	fp[e,t] <= rate_a[e])

 	@constraint(m, TrPowerFlowLimitDown[e=intersect(line_limit,dn_ext_edges),t=times],
 	-fp[e,t] <= rate_a[e])
        @constraint(m,fpFix[e=dn_ext_edges,t=timing], fp[e,t]-exp_flow[e,t] == 0)
        else
                for e=dn_ext_edges for t=timing
                        set_normalized_rhs(m.obj_dict[:fpFix][e,t],exp_flow[e,t])
                end end
        end
        #@constraint(m, fpFixUp[e=dn_ext_edges,t=times], sign(exp_flow[e,t])*fp[e,t] >= (1-sign(exp_flow[e,t])*0.001)*exp_flow[e,t])
        #@constraint(m, fpFixDown[e=dn_ext_edges,t=times], sign(exp_flow[e,t])*fp[e,t] <= (1+sign(exp_flow[e,t])*0.001)*exp_flow[e,t])

         optimize!(m)
         #println(primal_status(m))
         # println(dual_status(m))
         # println(raw_status(m))
         #println(termination_status(m),"\n")
         if termination_status(m)!=MOI.OPTIMAL#primal_status(m)!=MOI.FEASIBLE_POINT
                 set_optimizer_attribute(m,"MSK_IPAR_PRESOLVE_USE",0)
                 optimize!(m)
                 set_optimizer_attribute(m,"MSK_IPAR_PRESOLVE_USE",1)
         end
         # println(primal_status(m))
         # println(termination_status(m),"\n")
         qty_price_interface=Dict{Tuple{Int64,Int64},Tuple{Float64,Float64}}()

        #qty_price_interface=Dict((e,t) => (termination_status(m),dual(m.obj_dict[:fpFix][e,t]),exp_flow[e,t]) for e=dn_ext_edges for t=times)
        #if primal_status(m)!=MOI.FEASIBLE_POINT
        if termination_status(m)!=MOI.OPTIMAL
                for e=dn_ext_edges for t=timing
                        #qty_price_interface[e,t]=(MOI.SLOW_PROGRESS,10*cost_ls,exp_flow[e,t])
                        qty_price_interface[e,t]=(10*cost_ls,exp_flow[e,t])
                end end
        else
                for e=dn_ext_edges for t=timing
                        qty_price_interface[e,t]=(dual(m.obj_dict[:fpFix][e,t]),exp_flow[e,t])
                end end
        end
         return (qty_price_interface,m,solve_time(m))
end

function compute_res_fun_bin_approx_t(dn::String,
                         exp_flow::Dict{Tuple{Int64,Int64},Float64},
                         timing::Int64,
                         m::Model=Model())
        dn_ext_edges=ext_edges_DN[dn]
        if isempty(m.obj_dict)
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
        dn_t_seg_bid=[sb for sb=dn_seg_bid if sb[5]==timing]
        dn_t_bid=[sb for sb=dn_bid if sb[4]==timing]
        dn_t_dist_seg_bid=[sb for sb=dn_dist_seg_bid if sb[5]==timing]
         @variable(m, 0 <= c[dn_c_index,timing])
         @variable(m, s[dn_s_index,timing])
         @variable(m, pg[dn_t_seg_bid])
         @expression(m, net_pg_bid[bd=dn_t_bid], sum(pg[sb] for sb=bid_2_seg_bids[bd]))
         @variable(m, qg[dn_seg_bid])
         @expression(m, net_qg_bid[bd=dn_t_bid], sum(qg[sb] for sb=bid_2_seg_bids[bd]))
         @variable(m, -pi <= theta[union(dn_nodes,dn_ext_nodes),timing] <= pi)
         if slack_bus in union(dn_nodes,dn_ext_nodes)
                 for i=slack_bus for t=timing fix(theta[i,t],0;force=true) end end
         end
         @variable(m, slackPUp[i=dn_nodes,t=timing] >= 0)
         @variable(m, slackPDown[i=dn_nodes,t=timing] >= 0)
         @variable(m, slackQUp[i=dn_nodes,t=timing] >= 0)
         @variable(m, slackQDown[i=dn_nodes,t=timing] >= 0)

         ### BINARY VARIABLES ###
 	@variable(m, 0 <= sa[dn_t_seg_bid] <= 1)

         ### FLOW DEFINITION ###
         @variable(m, fp[union(dn_edges,dn_ext_edges),timing])
         @variable(m, fpTo[union(dn_edges,dn_ext_edges),timing])

         DefFpDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
         DefFpToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
         for l=dn_edges for t=timing
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

         @variable(m, fq[dn_edges,timing])
         @variable(m, fqTo[dn_edges,timing])

         DefFqDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
         DefFqToDist=Dict{Tuple{Int64,Int64},ConstraintRef}()
         for l=dn_edges for t=timing
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
         @constraint(m, PDown[sb=dn_t_seg_bid], -pg[sb] <= -sa[sb]*p_min[sb])
         @constraint(m, PUp[sb=dn_t_seg_bid], pg[sb] <= sa[sb]*p_max[sb])

         @constraint(m, PowerBalanceP[i=dn_nodes,t=timing],
         sum(sum(pg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
                 for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
         + net_injP[i,t]
         + slackPUp[i,t]
         ==
         + sum(fpTo[e,t] for e=n_to[i])
         + sum(fp[e,t] for e=n_fr[i])
         + slackPDown[i,t]
         )

         @constraint(m, PowerBalanceQ[i=dn_nodes,t=timing],
         sum(sum(qg[sb] for sb=get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
                 for B=get(nt_2_bids,(i,t),Set{Tuple{Int64,Int64}}()))
         + net_injQ[i,t]
         + slackQUp[i,t]
         ==
         + sum(fqTo[e,t] for e=dist_n_to[i])
         + sum(fq[e,t] for e=dist_n_fr[i])
         + slackQDown[i,t]
         )

         @constraint(m, VoltageLimitsDown[i=intersect(dn_nodes,v_limit),t=timing],
         -c[(i,i),t] <= -vsq_min[i])

         @constraint(m, VoltageLimitsUp[i=intersect(dn_nodes,v_limit),t=timing],
         c[(i,i),t] <= vsq_max[i])

         @constraint(m, LinePowerLimit1[e=intersect(line_limit,dn_edges),t=timing],
         [rate_a[e],fp[e,t],fq[e,t]] in SecondOrderCone())

         @constraint(m, LinePowerLimit2[e=intersect(line_limit,dn_edges),t=timing],
         [rate_a[e],fpTo[e,t],fqTo[e,t]] in SecondOrderCone())

         mod_s_ind=[(e_fr[e],e_to[e]) for e=dn_edges if R_pu[e]!=0.0]
        @constraint(m, SOCPconstraint[e=mod_s_ind,t=timing],
        [0.5*c[(e[1],e[1]),t],c[(e[2],e[2]),t],c[e,t],s[e,t]] in RotatedSecondOrderCone())

         @constraint(m, DiscQtoP[BQpd=[qpd for qpd=QPDiscConstrIndex if qpd[1] in dn_nodes && qpd[4]==timing]],
        [get(QPDiscMax,(BQpd[2],BQpd[3],BQpd[5]),0),
        net_pg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])],
        net_qg_bid[(BQpd[1],BQpd[2],BQpd[3],BQpd[4])]]
        in SecondOrderCone()
        )

         @constraint(m, ConstrainToUpFromLineQtoP[BHpc=[qpd for qpd=QPHPConstrIndex1 if qpd[1] in dn_nodes && qpd[4]==timing]],
         net_qg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
         -QPHPOffset[(BHpc[5],BHpc[2])]
         -QPHPSlope[(BHpc[5],BHpc[2])]*net_pg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
         >= 0)

         @constraint(m, ConstrainToDownFromLineQtoP[BHpc=[qpd for qpd=QPHPConstrIndex0 if qpd[1] in dn_nodes && qpd[4]==timing]],
         QPHPOffset[(BHpc[5],BHpc[2])]
         +QPHPSlope[(BHpc[5],BHpc[2])]*net_pg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])]
         -net_qg_bid[(BHpc[1],BHpc[2],BHpc[3],BHpc[4])] >= 0)

         ### Flow definition for interconnection ###
         @constraint(m, DefFpTr[e=dn_ext_edges,t=timing],
         X_pu[e]*fp[e,t] #+ slackFUp[e,t] - slackFDown[e,t]
         == (theta[e_fr[e],t]-theta[e_to[e],t]))
         @constraint(m, DefFpX0[e=dn_ext_edges, t=timing],
         fp[e,t] + fpTo[e,t] == 0)

         # ### OBJECTIVE ###
         @objective(m, Min,
         + sum(cost[sb][1]*pg[sb] for sb=dn_t_seg_bid)
         + cost_ls*(sum(slackPUp[i,t] + slackPDown[i,t]
                         + slackQUp[i,t] + slackQDown[i,t] for i=dn_nodes for t=timing)
                #+sum(slackFUp[e,t] + slackFDown[e,t] for e=dn_ext_edges for t=times)
                        )
         )

         ### Flow at the interconnection is fixed to a given value ###
         @constraint(m, TrPowerFlowLimitUp[e=intersect(line_limit,dn_ext_edges),t=timing],
 	fp[e,t] <= rate_a[e])

 	@constraint(m, TrPowerFlowLimitDown[e=intersect(line_limit,dn_ext_edges),t=timing],
 	-fp[e,t] <= rate_a[e])
        @constraint(m,fpFix[e=dn_ext_edges,t=timing], fp[e,t]-exp_flow[e,t] == 0)
        else
                for e=dn_ext_edges for t=timing
                        set_normalized_rhs(m.obj_dict[:fpFix][e,t],exp_flow[e,t])
                end end
        end
         optimize!(m)

         if termination_status(m)!=MOI.OPTIMAL#primal_status(m)!=MOI.FEASIBLE_POINT
                 set_optimizer_attribute(m,"MSK_IPAR_PRESOLVE_USE",0)
                 optimize!(m)
                 set_optimizer_attribute(m,"MSK_IPAR_PRESOLVE_USE",1)
         end

         qty_price_interface=Dict{Tuple{Int64,Int64},Tuple{Float64,Float64}}()

        if termination_status(m)!=MOI.OPTIMAL
                for e=dn_ext_edges for t=timing
                        qty_price_interface[e,t]=(10*cost_ls,exp_flow[e,t])
                end end
        else
                for e=dn_ext_edges for t=timing
                        qty_price_interface[e,t]=(dual(m.obj_dict[:fpFix][e,t]),exp_flow[e,t])
                end end
        end
         return (qty_price_interface,m,solve_time(m))
end

function compute_res_fun_bin_exact_t(dn::String,
                         exp_flow::Dict{Tuple{Int64,Int64},Float64},
                         timing::Int64,
                         rsf,nb_seg::Int64,
                         m::Model=Model())
        dn_ext_edges=ext_edges_DN[dn]
        if isempty(m.obj_dict)
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
         if slack_bus in union(dn_nodes,dn_ext_nodes)
                 for i=slack_bus for t=times fix(theta[i,t],0;force=true) end end
         end
         #for i=dn_ext_nodes for t=times fix(theta[i,t],0;force=true) end end
         #@constraint(m, fixTheta0[i=dn_ext_nodes,t=times],theta[i,t]==0)

         # @variable(m, abs(net_injP[i,t]) >= slackPUp[i=dn_nodes,t=times] >= 0)
         # @variable(m, abs(net_injP[i,t]) >= slackPDown[i=dn_nodes,t=times] >= 0)
         # @variable(m, abs(net_injQ[i,t]) >= slackQUp[i=dn_nodes,t=times] >= 0)
         # @variable(m, abs(net_injQ[i,t]) >= slackQDown[i=dn_nodes,t=times] >= 0)
         @variable(m, slackPUp[i=dn_nodes,t=times] >= 0)
         @variable(m, slackPDown[i=dn_nodes,t=times] >= 0)
         @variable(m, slackQUp[i=dn_nodes,t=times] >= 0)
         @variable(m, slackQDown[i=dn_nodes,t=times] >= 0)
         #@variable(m, slackFUp[dn_ext_edges,times] >= 0)
         #@variable(m, slackFDown[dn_ext_edges,times] >= 0)

         ### BINARY VARIABLES ###
 	@variable(m, 0 <= sa[dn_seg_bid] <= 1)
 	@variable(m, 0 <= qa[dn_bid] <= 1)
 	@variable(m, 0 <= qta[dn_qt_bid] <= 1)
 	@variable(m, 0 <= alpha[[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]] <= 1)
 	@variable(m, 0 <= omega[[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]] <= 1)

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
         @constraint(m, PDown[sb=dn_seg_bid], -pg[sb] <= -sa[sb]*p_min[sb])
         @constraint(m, PUp[sb=dn_seg_bid], pg[sb] <= sa[sb]*p_max[sb])


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
         + sum(net_pg_bid[bd] for bd=dn_bid if bd[2]==rcI[1] && bd[4]==rcI[3]+1)
         - sum(net_pg_bid[bd] for bd=dn_bid if bd[2]==rcI[1] && bd[3]==rcI[2] && bd[4]==rcI[3])
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
         X_pu[e]*fp[e,t] #+ slackFUp[e,t] - slackFDown[e,t]
         == (theta[e_fr[e],t]-theta[e_to[e],t]))
         @constraint(m, DefFpX0[e=dn_ext_edges, t=times],
         fp[e,t] + fpTo[e,t] == 0)

         ### BINARY CONSTRAINTS ###
 	@constraint(m, ExQtBidsConstr[i=dn_ex_qt_bid_id],
 	sum(qta[ExQt] for ExQt=ex_qt_bid[i]) <= 1)

        ## this constraint is not really necessary,
        ## I also noticed that it might lead to weird numerical behavior
        ## in terms of feasiblity and dual value returned
 	# @constraint(m, QtaDefLow[qt=dn_qt_bid],
 	# qta[qt] <= sum(qa[bd] for bd=dn_bid if bd[2]==qt))

 	@constraint(m, QtaDefHigh[bd=dn_bid],
 	qa[bd] <= qta[bd[2]])

 	@constraint(m, DefQa[sb=dn_seg_bid],
 	sa[sb] <= qa[(sb[1],sb[2],sb[3],sb[5])])

 	@constraint(m, AlphaOmegaDef[B=[aol for aol=AlphaOmegaLate if aol[1] in dn_nodes]],
 	qa[B]
 	- sum(qa[bd] for bd=dn_bid if bd[1]==B[1] && bd[1]==B[2] && bd[4]==B[4]-1)
 	- alpha[B]
 	+ omega[B]
 	== 0)

 	@constraint(m, MutexAlphaOmega[B=[aoi for aoi=AlphaOmegaIndex if aoi[1] in dn_nodes]],
 	alpha[B] + omega[B] <= 1)

 	@constraint(m, AlphaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in dn_nodes]],
 	qa[B] - alpha[B] <= 0)

 	@constraint(m, OmegaInit[B=[aoe for aoe=AlphaOmegaEarly if aoe[1] in dn_nodes]],
 	omega[B] == 0)

 	@constraint(m, MinDurationConstr[mdp=dn_MinDurIndex],
 	+ sum(qa[bd] for bd=dn_bid if bd[2]==mdp[1] && bd[3]==mdp[3] && bd[4]==mdp[5])
  	- sum(alpha[bd] for bd=dn_bid if bd[2]==mdp[1] && bd[3]==mdp[2] && bd[4]==mdp[4])
 	>= 0)
        #println(MinDurationConstr)

 	@constraint(m, NoNewAct[B=[nnai for nnai=NoNewActIndex if nnai[1] in dn_nodes]],
 	alpha[B] == 0)

        ### Flows defined as a function of segments ###
        @variable(m, fp_n[1:nb_seg,dn_ext_edges,setdiff(times,timing)])
        @variable(m, 0 <= x_n[1:nb_seg,dn_ext_edges,setdiff(times,timing)] <= 1)

        @constraint(m, DefFp_n[n=1:nb_seg,e=dn_ext_edges,t=setdiff(times,timing)],
        fp_n[n,e,t] == x_n[n,e,t]*(rsf[e,t][n][2][2]-rsf[e,t][n][2][1]))

        @constraint(m, def_fp_n[e=dn_ext_edges,t=setdiff(times,timing)],
        fp[e,t] == sum(fp_n[n,e,t] for n=1:nb_seg) + rsf[e,t][1][2][1])

         # ### OBJECTIVE ###
         @objective(m, Min,
         + sum(cost[sb][1]*pg[sb] for sb=dn_seg_bid)
         + cost_ls*(sum(slackPUp[i,t] + slackPDown[i,t]
                         + slackQUp[i,t] + slackQDown[i,t] for i=dn_nodes for t=times)
                #+sum(slackFUp[e,t] + slackFDown[e,t] for e=dn_ext_edges for t=times)
                        )
        - sum(fp_n[n,e,t]*rsf[e,t][n][1] for n=1:nb_seg for e=dn_ext_edges for t=setdiff(times,timing))
         )

         ### Flow at the interconnection is fixed to a given value ###
         @constraint(m, TrPowerFlowLimitUp[e=intersect(line_limit,dn_ext_edges),t=times],
 	fp[e,t] <= rate_a[e])

 	@constraint(m, TrPowerFlowLimitDown[e=intersect(line_limit,dn_ext_edges),t=times],
 	-fp[e,t] <= rate_a[e])
        @constraint(m,fpFix[e=dn_ext_edges,t=timing], fp[e,t]-exp_flow[e,t] == 0)
        else
                for e=dn_ext_edges for t=timing
                        set_normalized_rhs(m.obj_dict[:fpFix][e,t],exp_flow[e,t])
                end end
        end
        #@constraint(m, fpFixUp[e=dn_ext_edges,t=times], sign(exp_flow[e,t])*fp[e,t] >= (1-sign(exp_flow[e,t])*0.001)*exp_flow[e,t])
        #@constraint(m, fpFixDown[e=dn_ext_edges,t=times], sign(exp_flow[e,t])*fp[e,t] <= (1+sign(exp_flow[e,t])*0.001)*exp_flow[e,t])

         optimize!(m)
         #println(primal_status(m))
         # println(dual_status(m))
         # println(raw_status(m))
         #println(termination_status(m),"\n")
         if termination_status(m)!=MOI.OPTIMAL#primal_status(m)!=MOI.FEASIBLE_POINT
                 set_optimizer_attribute(m,"MSK_IPAR_PRESOLVE_USE",0)
                 optimize!(m)
                 set_optimizer_attribute(m,"MSK_IPAR_PRESOLVE_USE",1)
         end
         # println(primal_status(m))
         # println(termination_status(m),"\n")
         qty_price_interface=Dict{Tuple{Int64,Int64},Tuple{Float64,Float64}}()

        #qty_price_interface=Dict((e,t) => (termination_status(m),dual(m.obj_dict[:fpFix][e,t]),exp_flow[e,t]) for e=dn_ext_edges for t=times)
        #if primal_status(m)!=MOI.FEASIBLE_POINT
        if termination_status(m)!=MOI.OPTIMAL
                for e=dn_ext_edges for t=timing
                        #qty_price_interface[e,t]=(MOI.SLOW_PROGRESS,10*cost_ls,exp_flow[e,t])
                        qty_price_interface[e,t]=(10*cost_ls,exp_flow[e,t])
                end end
        else
                for e=dn_ext_edges for t=timing
                        qty_price_interface[e,t]=(dual(m.obj_dict[:fpFix][e,t]),exp_flow[e,t])
                end end
        end
         return (qty_price_interface,m,solve_time(m))
end

function is_isolated(dn::String)::Bool
        if length(DNs[dn])==1 #&& length(seg_bid_DN[dn])==0
                return true
        end
        return false
end

function compute_n_points_bin(nb_point::Int64,exp_flow::Dict{Tuple{Int64,Int64},Float64})
        solve_times=Dict{Tuple{String,Int64},Float64}()
        println("Starting decomposition...")
	price_qty_DN=Dict{Tuple{Int64,Int64},Any}((e,t)=> [] for e=ext_edges for t=times)
	println("Starting $(length(root_dist_node)) Distribution networks problems.")
	for (num,dn) in enumerate(keys(DNs))
                #if dn=="DN-22650039"
                st_array=[]
                inter_e=ext_edges_DN[dn]
                inter_e=first(inter_e)
                copy_exp_flow=copy(exp_flow)
                maxFlow=0.0
                model_dn=Model()
                if !(inter_e in line_limit) maxFlow=maximum(abs,[exp_flow[inter_e,t] for t=times]) end
                if is_isolated(dn)
                        for t=times copy_exp_flow[inter_e,t]=net_injP[first(DNs[dn]),t] end
                        (tmp_price_qty,model_dn,tmp_st)=compute_res_fun_bin(dn,copy_exp_flow,model_dn)
                        for n=0:(nb_point+1)
                                solve_times[dn,n]=0.0
                                for i=keys(tmp_price_qty)
                                        push!(price_qty_DN[i],tmp_price_qty[i])
                                end
                        end
                        solve_times[dn,1]=tmp_st
                else
                        for n=0:(nb_point+1)
                                if inter_e in line_limit
                                        qty=-rate_a[inter_e]+n*(2*rate_a[inter_e]/(nb_point+1))
                                        for t=times copy_exp_flow[inter_e,t]=qty end
                                else
                                        #frac=-1.5+3*n/(nb_point+1)
                                        frac=-2+4*n/(nb_point+1)
                                        for t=times copy_exp_flow[inter_e,t]=frac*maxFlow end
                                end
                                (tmp_price_qty,model_dn,solve_times[dn,n])=compute_res_fun_bin(dn,copy_exp_flow,model_dn)
                                for i=keys(tmp_price_qty)
                                        push!(price_qty_DN[i],tmp_price_qty[i])
                                end
                        end
                end
		println(">$num. SOLVED Distribution network problem: $dn.")
                #end
	end
        price_qty=transform_price_qty(nb_point+2,price_qty_DN)
	return (price_qty,solve_times)
        #return (price_qty_DN,solve_times)
end

function compute_n_points_bin_oneill(nb_point::Int64,exp_flow::Dict{Tuple{Int64,Int64},Float64})
        solve_times=Dict{Tuple{String,Int64},Float64}()
        println("Starting decomposition...")
	price_qty_DN=Dict{Tuple{Int64,Int64},Any}((e,t)=> [] for e=ext_edges for t=times)
	println("Starting $(length(root_dist_node)) Distribution networks problems.")
	for (num,dn) in enumerate(keys(DNs))
                #if dn=="DN-22650039"
                st_array=[]
                inter_e=ext_edges_DN[dn]
                inter_e=first(inter_e)
                copy_exp_flow=copy(exp_flow)
                maxFlow=0.0
                model_dn=Model()
                relax_model_dn=Model()
                if !(inter_e in line_limit) maxFlow=maximum(abs,[exp_flow[inter_e,t] for t=times]) end
                if is_isolated(dn)
                        for t=times copy_exp_flow[inter_e,t]=net_injP[first(DNs[dn]),t] end
                        (tmp_price_qty,model_dn,relax_model_dn,tmp_st)=compute_res_fun_bin_oneill(dn,copy_exp_flow,model_dn,relax_model_dn)
                        for n=0:(nb_point+1)
                                solve_times[dn,n]=0.0
                                for i=keys(tmp_price_qty)
                                        push!(price_qty_DN[i],tmp_price_qty[i])
                                end
                        end
                        solve_times[dn,1]=tmp_st
                else
                        for n=0:(nb_point+1)
                                if inter_e in line_limit
                                        qty=-rate_a[inter_e]+n*(2*rate_a[inter_e]/(nb_point+1))
                                        for t=times copy_exp_flow[inter_e,t]=qty end
                                else
                                        #frac=-1.5+3*n/(nb_point+1)
                                        frac=-2+4*n/(nb_point+1)
                                        for t=times copy_exp_flow[inter_e,t]=frac*maxFlow end
                                end
                                (tmp_price_qty,model_dn,relax_model_dn,solve_times[dn,n])=compute_res_fun_bin_oneill(dn,copy_exp_flow,model_dn,relax_model_dn)
                                for i=keys(tmp_price_qty)
                                        push!(price_qty_DN[i],tmp_price_qty[i])
                                end
                        end
                end
		println(">$num. SOLVED Distribution network problem: $dn.")
                #end
	end
        price_qty=transform_price_qty(nb_point+2,price_qty_DN)
	return (price_qty,solve_times)
        #return (price_qty_DN,solve_times)
end

function compute_n_points_bin_t(nb_point::Int64,exp_flow::Dict{Tuple{Int64,Int64},Float64})
        solve_times=Dict((dn,n)=>0.0 for dn=keys(DNs) for n=0:(nb_point+1))
        println("Starting decomposition...")
	price_qty_DN=Dict{Tuple{Int64,Int64},Any}((e,t)=> [] for e=ext_edges for t=times)
	println("Starting $(length(root_dist_node)) Distribution networks problems.")
	for (num,dn) in enumerate(keys(DNs))
                #if dn=="DN-22650039"
                st_array=[]
                inter_e=ext_edges_DN[dn]
                inter_e=first(inter_e)
                copy_exp_flow=copy(exp_flow)
                maxFlow=0.0
                if !(inter_e in line_limit) maxFlow=maximum(abs,[exp_flow[inter_e,t] for t=times]) end
                if is_isolated(dn)
                        for t=times copy_exp_flow[inter_e,t]=net_injP[first(DNs[dn]),t] end
                        (tmp_price_qty,model_dn,tmp_st)=compute_res_fun_bin(dn,copy_exp_flow,model_dn)
                        for n=0:(nb_point+1)
                                solve_times[dn,n]=0.0
                                for i=keys(tmp_price_qty)
                                        push!(price_qty_DN[i],tmp_price_qty[i])
                                end
                        end
                        solve_times[dn,1]=tmp_st
                else
                        model_dn=Dict(t=>Model() for t=times)
                        for n=0:(nb_point+1)
                                if inter_e in line_limit
                                        qty=-rate_a[inter_e]+n*(2*rate_a[inter_e]/(nb_point+1))
                                        for t=times copy_exp_flow[inter_e,t]=qty end
                                else
                                        #frac=-1.5+3*n/(nb_point+1)
                                        frac=-2+4*n/(nb_point+1)
                                        for t=times copy_exp_flow[inter_e,t]=frac*maxFlow end
                                end
                                tmp_price_qty=Dict{Int64,Any}()
                                for t=times
                                (tmp_price_qty[t],model_dn[t],tmp_st)=compute_res_fun_bin_t(dn,copy_exp_flow,t,model_dn[t])
                                solve_times[dn,n]+=tmp_st
                                end
                                for t=times for i=keys(tmp_price_qty[t])
                                        push!(price_qty_DN[i],tmp_price_qty[t][i])
                                end end
                        end
                end
		println(">$num. SOLVED Distribution network problem: $dn.")
                #end
	end
        price_qty=transform_price_qty(nb_point+2,price_qty_DN)
	return (price_qty,solve_times)
        #return (price_qty_DN,solve_times)
end

function compute_n_points_bin_approx_t(nb_point::Int64,exp_flow::Dict{Tuple{Int64,Int64},Float64})
        solve_times=Dict((dn,n)=>0.0 for dn=keys(DNs) for n=0:(nb_point+1))
        println("Starting decomposition...")
	price_qty_DN=Dict{Tuple{Int64,Int64},Any}((e,t)=> [] for e=ext_edges for t=times)
	println("Starting $(length(root_dist_node)) Distribution networks problems.")
	for (num,dn) in enumerate(keys(DNs))
                #if dn=="DN-22650039"
                st_array=[]
                inter_e=ext_edges_DN[dn]
                inter_e=first(inter_e)
                copy_exp_flow=copy(exp_flow)
                maxFlow=0.0
                if !(inter_e in line_limit) maxFlow=maximum(abs,[exp_flow[inter_e,t] for t=times]) end
                model_dn=Dict(t=>Model() for t=times)
                for n=0:(nb_point+1)
                        if inter_e in line_limit
                                qty=-rate_a[inter_e]+n*(2*rate_a[inter_e]/(nb_point+1))
                                for t=times copy_exp_flow[inter_e,t]=qty end
                        else
                                #frac=-1.5+3*n/(nb_point+1)
                                frac=-2+4*n/(nb_point+1)
                                for t=times copy_exp_flow[inter_e,t]=frac*maxFlow end
                        end
                        tmp_price_qty=Dict{Int64,Any}()
                        for t=times
                        (tmp_price_qty[t],model_dn[t],tmp_st)=compute_res_fun_bin_approx_t(dn,copy_exp_flow,t,model_dn[t])
                        solve_times[dn,n]+=tmp_st
                        end
                        for t=times for i=keys(tmp_price_qty[t])
                                push!(price_qty_DN[i],tmp_price_qty[t][i])
                        end end
                end
		println(">$num. SOLVED Distribution network problem: $dn.")
                #end
	end
        price_qty=transform_price_qty(nb_point+2,price_qty_DN)
	return (price_qty,solve_times)
        #return (price_qty_DN,solve_times)
end

function compute_n_points_bin_exact_t(nb_point::Int64,rsf,exp_flow::Dict{Tuple{Int64,Int64},Float64})
        solve_times=Dict((dn,n)=>0.0 for dn=keys(DNs) for n=0:(nb_point+1))
        println("Starting decomposition...")
	price_qty_DN=Dict{Tuple{Int64,Int64},Any}((e,t)=> [] for e=ext_edges for t=times)
	println("Starting $(length(root_dist_node)) Distribution networks problems.")
	for (num,dn) in enumerate(keys(DNs))
                #if dn=="DN-22650039"
                st_array=[]
                inter_e=ext_edges_DN[dn]
                inter_e=first(inter_e)
                copy_exp_flow=copy(exp_flow)
                maxFlow=0.0
                if !(inter_e in line_limit) maxFlow=maximum(abs,[exp_flow[inter_e,t] for t=times]) end
                model_dn=Dict(t=>Model() for t=times)
                for n=0:(nb_point+1)
                        if inter_e in line_limit
                                qty=-rate_a[inter_e]+n*(2*rate_a[inter_e]/(nb_point+1))
                                for t=times copy_exp_flow[inter_e,t]=qty end
                        else
                                #frac=-1.5+3*n/(nb_point+1)
                                frac=-2+4*n/(nb_point+1)
                                for t=times copy_exp_flow[inter_e,t]=frac*maxFlow end
                        end
                        tmp_price_qty=Dict{Int64,Any}()
                        for t=times
                        (tmp_price_qty[t],model_dn[t],tmp_st)=compute_res_fun_bin_exact_t(dn,copy_exp_flow,t,rsf,nb_point+1,model_dn[t])
                        solve_times[dn,n]+=tmp_st
                        end
                        for t=times for i=keys(tmp_price_qty[t])
                                push!(price_qty_DN[i],tmp_price_qty[t][i])
                        end end
                end
		println(">$num. SOLVED Distribution network problem: $dn.")
                #end
	end
        price_qty=transform_price_qty(nb_point+2,price_qty_DN)
	return (price_qty,solve_times)
        #return (price_qty_DN,solve_times)
end

function add_points_price_qty(nb_point::Int64, price_qty)
        new_price_qty=copy(price_qty)
        solve_times=Dict{Tuple{String,Int64},Float64}()
        for (num,dn) in enumerate(keys(DNs))
                e=ext_edges_DN[dn]
                e=first(e)
                up_bound=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                low_bound=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                insert_n=Dict{Tuple{Int64,Int64},Int64}((e,t) => 3 for t=times)
                model_dn=Model()
                if is_isolated(dn)
                        for t=times for n=0:nb_point
                                push!(new_price_qty[e,t],price_qty[e,t][1])
                        end end
                else
                for t=times
                        for n=2:length(price_qty[e,t])
                                if abs(price_qty[e,t][n][1]) > 1e3 && abs(price_qty[e,t][n-1][1]) < 1e3
                                        insert_n[e,t]=n
                                        continue
                                end
                        end
                        up_bound[e,t]=price_qty[e,t][insert_n[e,t]][2][2]
                        low_bound[e,t]=price_qty[e,t][insert_n[e,t]-2][2][1]
                        new_price_qty[e,t][insert_n[e,t]]=(price_qty[e,t][insert_n[e,t]][1],(up_bound[e,t] - (up_bound[e,t]-low_bound[e,t])/(nb_point+3),up_bound[e,t]))
                        deleteat!(new_price_qty[e,t],insert_n[e,t]-1)
                        deleteat!(new_price_qty[e,t],insert_n[e,t]-2)
                end
                qty_hi=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                qty_lo=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                for n=2:(nb_point+3)
                        for t=times
                                qty_hi[e,t]=up_bound[e,t] - (n-1)*(up_bound[e,t]-low_bound[e,t])/(nb_point+3)
                                qty_lo[e,t]=up_bound[e,t] - n*(up_bound[e,t]-low_bound[e,t])/(nb_point+3)
                        end
                        (tmp_price_qty,model_dn,solve_times[dn,n])=compute_res_fun_bin(dn,qty_hi,model_dn)
                        for t=times
                                insert!(new_price_qty[e,t],insert_n[e,t]-2,(tmp_price_qty[e,t][1],(qty_lo[e,t],qty_hi[e,t])))
                        end
                end
                end
                println(">$num. SOLVED Distribution network problem: $dn.")
        end
        return (new_price_qty,solve_times)
end

function add_points_price_qty_oneill(nb_point::Int64, price_qty)
        new_price_qty=copy(price_qty)
        solve_times=Dict{Tuple{String,Int64},Float64}()
        for (num,dn) in enumerate(keys(DNs))
                e=ext_edges_DN[dn]
                e=first(e)
                up_bound=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                low_bound=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                insert_n=Dict{Tuple{Int64,Int64},Int64}((e,t) => 3 for t=times)
                model_dn=Model()
                relax_model_dn=Model()
                if is_isolated(dn)
                        for t=times for n=0:nb_point
                                push!(new_price_qty[e,t],price_qty[e,t][1])
                        end end
                else
                for t=times
                        for n=2:length(price_qty[e,t])
                                if abs(price_qty[e,t][n][1]) > 1e3 && abs(price_qty[e,t][n-1][1]) < 1e3
                                        insert_n[e,t]=n
                                        continue
                                end
                        end
                        up_bound[e,t]=price_qty[e,t][insert_n[e,t]][2][2]
                        low_bound[e,t]=price_qty[e,t][insert_n[e,t]-2][2][1]
                        new_price_qty[e,t][insert_n[e,t]]=(price_qty[e,t][insert_n[e,t]][1],(up_bound[e,t] - (up_bound[e,t]-low_bound[e,t])/(nb_point+3),up_bound[e,t]))
                        deleteat!(new_price_qty[e,t],insert_n[e,t]-1)
                        deleteat!(new_price_qty[e,t],insert_n[e,t]-2)
                end
                qty_hi=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                qty_lo=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                for n=2:(nb_point+3)
                        for t=times
                                qty_hi[e,t]=up_bound[e,t] - (n-1)*(up_bound[e,t]-low_bound[e,t])/(nb_point+3)
                                qty_lo[e,t]=up_bound[e,t] - n*(up_bound[e,t]-low_bound[e,t])/(nb_point+3)
                        end
                        (tmp_price_qty,model_dn,relax_model_dn,solve_times[dn,n])=compute_res_fun_bin_oneill(dn,qty_hi,model_dn,relax_model_dn)
                        for t=times
                                insert!(new_price_qty[e,t],insert_n[e,t]-2,(tmp_price_qty[e,t][1],(qty_lo[e,t],qty_hi[e,t])))
                        end
                end
                end
                println(">$num. SOLVED Distribution network problem: $dn.")
        end
        return (new_price_qty,solve_times)
end

function add_points_price_qty_t(nb_point::Int64, price_qty)
        new_price_qty=copy(price_qty)
        solve_times=Dict((dn,n)=>0.0 for dn=keys(DNs) for n=2:(nb_point+3))
        for (num,dn) in enumerate(keys(DNs))
                e=ext_edges_DN[dn]
                e=first(e)
                up_bound=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                low_bound=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                insert_n=Dict{Tuple{Int64,Int64},Int64}((e,t) => 3 for t=times)
                if is_isolated(dn)
                        for t=times for n=0:nb_point
                                push!(new_price_qty[e,t],price_qty[e,t][1])
                        end end
                else
                for t=times
                        for n=2:length(price_qty[e,t])
                                if abs(price_qty[e,t][n][1]) > 1e3 && abs(price_qty[e,t][n-1][1]) < 1e3
                                        insert_n[e,t]=n
                                        continue
                                end
                        end
                        up_bound[e,t]=price_qty[e,t][insert_n[e,t]][2][2]
                        low_bound[e,t]=price_qty[e,t][insert_n[e,t]-2][2][1]
                        new_price_qty[e,t][insert_n[e,t]]=(price_qty[e,t][insert_n[e,t]][1],(up_bound[e,t] - (up_bound[e,t]-low_bound[e,t])/(nb_point+3),up_bound[e,t]))
                        deleteat!(new_price_qty[e,t],insert_n[e,t]-1)
                        deleteat!(new_price_qty[e,t],insert_n[e,t]-2)
                end
                qty_hi=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                qty_lo=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                model_dn=Dict(t=>Model() for t=times)
                for n=2:(nb_point+3)
                        for t=times
                                qty_hi[e,t]=up_bound[e,t] - (n-1)*(up_bound[e,t]-low_bound[e,t])/(nb_point+3)
                                qty_lo[e,t]=up_bound[e,t] - n*(up_bound[e,t]-low_bound[e,t])/(nb_point+3)
                        end
                        tmp_price_qty=Dict{Int64,Any}()
                        for t=times
                                (tmp_price_qty[t],model_dn[t],tmp_st)=compute_res_fun_bin_t(dn,qty_hi,t,model_dn[t])
                                insert!(new_price_qty[e,t],insert_n[e,t]-2,(tmp_price_qty[t][e,t][1],(qty_lo[e,t],qty_hi[e,t])))
                                solve_times[dn,n]+=tmp_st
                        end
                        # for t=times
                        #         insert!(new_price_qty[e,t],insert_n[e,t]-2,(tmp_price_qty[e,t][1],(qty_lo[e,t],qty_hi[e,t])))
                        # end
                end
                end
                println(">$num. SOLVED Distribution network problem: $dn.")
        end
        return (new_price_qty,solve_times)
end

function add_points_price_qty_approx_t(nb_point::Int64, price_qty)
        new_price_qty=copy(price_qty)
        solve_times=Dict((dn,n)=>0.0 for dn=keys(DNs) for n=2:(nb_point+3))
        for (num,dn) in enumerate(keys(DNs))
                e=ext_edges_DN[dn]
                e=first(e)
                up_bound=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                low_bound=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                insert_n=Dict{Tuple{Int64,Int64},Int64}((e,t) => 3 for t=times)
                for t=times
                        for n=2:length(price_qty[e,t])
                                if abs(price_qty[e,t][n][1]) > 1e3 && abs(price_qty[e,t][n-1][1]) < 1e3
                                        insert_n[e,t]=n
                                        continue
                                end
                        end
                        up_bound[e,t]=price_qty[e,t][insert_n[e,t]][2][2]
                        low_bound[e,t]=price_qty[e,t][insert_n[e,t]-2][2][1]
                        new_price_qty[e,t][insert_n[e,t]]=(price_qty[e,t][insert_n[e,t]][1],(up_bound[e,t] - (up_bound[e,t]-low_bound[e,t])/(nb_point+3),up_bound[e,t]))
                        deleteat!(new_price_qty[e,t],insert_n[e,t]-1)
                        deleteat!(new_price_qty[e,t],insert_n[e,t]-2)
                end
                qty_hi=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                qty_lo=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                model_dn=Dict(t=>Model() for t=times)
                for n=2:(nb_point+3)
                        for t=times
                                qty_hi[e,t]=up_bound[e,t] - (n-1)*(up_bound[e,t]-low_bound[e,t])/(nb_point+3)
                                qty_lo[e,t]=up_bound[e,t] - n*(up_bound[e,t]-low_bound[e,t])/(nb_point+3)
                        end
                        tmp_price_qty=Dict{Int64,Any}()
                        for t=times
                                (tmp_price_qty[t],model_dn[t],tmp_st)=compute_res_fun_bin_approx_t(dn,qty_hi,t,model_dn[t])
                                insert!(new_price_qty[e,t],insert_n[e,t]-2,(tmp_price_qty[t][e,t][1],(qty_lo[e,t],qty_hi[e,t])))
                                solve_times[dn,n]+=tmp_st
                        end
                        # for t=times
                        #         insert!(new_price_qty[e,t],insert_n[e,t]-2,(tmp_price_qty[e,t][1],(qty_lo[e,t],qty_hi[e,t])))
                        # end
                end
                println(">$num. SOLVED Distribution network problem: $dn.")
        end
        return (new_price_qty,solve_times)
end

function add_points_price_qty_exact_t(nb_point::Int64, rsf, price_qty)
        new_price_qty=copy(price_qty)
        solve_times=Dict((dn,n)=>0.0 for dn=keys(DNs) for n=2:(nb_point+3))
        for (num,dn) in enumerate(keys(DNs))
                e=ext_edges_DN[dn]
                e=first(e)
                up_bound=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                low_bound=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                insert_n=Dict{Tuple{Int64,Int64},Int64}((e,t) => 3 for t=times)
                for t=times
                        for n=2:length(price_qty[e,t])
                                if abs(price_qty[e,t][n][1]) > 1e3 && abs(price_qty[e,t][n-1][1]) < 1e3
                                        insert_n[e,t]=n
                                        continue
                                end
                        end
                        up_bound[e,t]=price_qty[e,t][insert_n[e,t]][2][2]
                        low_bound[e,t]=price_qty[e,t][insert_n[e,t]-2][2][1]
                        new_price_qty[e,t][insert_n[e,t]]=(price_qty[e,t][insert_n[e,t]][1],(up_bound[e,t] - (up_bound[e,t]-low_bound[e,t])/(nb_point+3),up_bound[e,t]))
                        deleteat!(new_price_qty[e,t],insert_n[e,t]-1)
                        deleteat!(new_price_qty[e,t],insert_n[e,t]-2)
                end
                qty_hi=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                qty_lo=Dict{Tuple{Int64,Int64},Float64}((e,t)=>0.0 for t=times)
                model_dn=Dict(t=>Model() for t=times)
                for n=2:(nb_point+3)
                        for t=times
                                qty_hi[e,t]=up_bound[e,t] - (n-1)*(up_bound[e,t]-low_bound[e,t])/(nb_point+3)
                                qty_lo[e,t]=up_bound[e,t] - n*(up_bound[e,t]-low_bound[e,t])/(nb_point+3)
                        end
                        tmp_price_qty=Dict{Int64,Any}()
                        for t=times
                                (tmp_price_qty[t],model_dn[t],tmp_st)=compute_res_fun_bin_exact_t(dn,qty_hi,t,rsf,nb_point+1,model_dn[t])
                                insert!(new_price_qty[e,t],insert_n[e,t]-2,(tmp_price_qty[t][e,t][1],(qty_lo[e,t],qty_hi[e,t])))
                                solve_times[dn,n]+=tmp_st
                        end
                        # for t=times
                        #         insert!(new_price_qty[e,t],insert_n[e,t]-2,(tmp_price_qty[e,t][1],(qty_lo[e,t],qty_hi[e,t])))
                        # end
                end
                println(">$num. SOLVED Distribution network problem: $dn.")
        end
        return (new_price_qty,solve_times)
end

function correct_price_qty(price_qty)
        # for e=ext_edges for t=times
        #         n_mod=0
        #         for n=1:(length(price_qty[e,t])-1)
        #                 if price_qty[e,t][n][1] > 1e3 && price_qty[e,t][n+1][1] < 5e2
        #                         n_mod=n
        #                         continue
        #                 end
        #         end
        #         if n_mod!=0
        #                 for n=1:n_mod
        #                         price_qty[e,t][n]=(-price_qty[e,t][n][1],price_qty[e,t][n][2])
        #                 end
        #         end
        # end end
        for e=ext_edges for t=times
                for n=2:(length(price_qty[e,t])-1)
                        if price_qty[e,t][n][1] > 1e3 && price_qty[e,t][n+1][1] < 5e2
                                n_bis=n
                                while price_qty[e,t][n_bis][1] > 1e3 && n_bis > 1
                                        price_qty[e,t][n_bis]=(price_qty[e,t][n+1][1],price_qty[e,t][n_bis][2])
                                        n_bis-=1
                                end
                        end
                end
        end end
        return price_qty
end

function solve_TN_res_sup_fun_bin(nb_seg,price_qty,binsol=0,solver::Symbol=:Mosek)
        m=Model(optimizer_with_attributes(Mosek.Optimizer,
        "MSK_IPAR_PRESOLVE_USE" => 1,
        "MSK_IPAR_LOG" => 1,
        "MSK_IPAR_INTPNT_BASIS" => 0,
        "MSK_DPAR_MIO_TOL_ABS_RELAX_INT"=>1e-6,
        "MSK_DPAR_INTPNT_CO_TOL_NEAR_REL"=>100.0
        ))
        if solver==:Gurobi
                m=Model(optimizer_with_attributes(Gurobi.Optimizer,
                "BarHomogeneous" => 1,
                "Presolve" => -1,
                "QCPDual" => 1,
                ))
        end

        ### CONTINUOUS VARIABLES ###
        @variable(m, pg[sb=tr_seg_bid])
        @expression(m, net_pg_bid[bd=tr_bid], sum(pg[sb] for sb=bid_2_seg_bids[bd]))
        @variable(m, -pi <= theta_tr[union(tr_nodes_ext,slack_bus),times] <= pi)
        @variable(m, slackPUp[tr_nodes,times] == 0)
        @variable(m, slackPDown[tr_nodes,times] == 0)

        if binsol==0
        @variable(m, sa[tr_seg_bid], Bin)
 	@variable(m, qa[tr_bid], Bin)
 	@variable(m, qta[tr_qt_bid], Bin)
 	@variable(m, alpha[[aoi for aoi=AlphaOmegaIndex if aoi[1] in tr_nodes]], Bin)
 	@variable(m, omega[[aoi for aoi=AlphaOmegaIndex if aoi[1] in tr_nodes]], Bin)
        end

        ### FLOW DEFINITION ###
        @variable(m, fp[tr_edges,times])
        @constraint(m, DefFpTr[e=tr_edges,t=times],
        X_pu[e]*fp[e,t] == (theta_tr[e_fr[e],t]-theta_tr[e_to[e],t]))

        @variable(m,fpTo[tr_edges,times])
        @constraint(m, DefFpX0[e=tr_edges,t=times],
        fp[e,t] + fpTo[e,t] == 0)

        ### CONSTRAINTS ###
        if binsol==0
        @constraint(m, PDown[sb=tr_seg_bid], -pg[sb] <= -sa[sb]*p_min[sb])
        @constraint(m, PUp[sb=tr_seg_bid], pg[sb] <= sa[sb]*p_max[sb])
        else
        @constraint(m, PDown[sb=tr_seg_bid], -pg[sb] <= -binsol.sa[sb]*p_min[sb])
        @constraint(m, PUp[sb=tr_seg_bid], pg[sb] <= binsol.sa[sb]*p_max[sb])
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

        @variable(m, fp_n[1:nb_seg,ext_edges,times])
        @variable(m, 0 <= x_n[1:nb_seg,ext_edges,times] <= 1)
        #@variable(m, 0 <= y_n[1:nb_seg,ext_edges,times] <= 1)
        #@constraint(m, Reform1[n=1:nb_seg,e=ext_edges,t=times], y_n[n,e,t] <= x_n[n,e,t])
        #@constraint(m, Reform2[n=1:nb_seg,e=ext_edges,t=times], y_n[n,e,t] >= 2*x_n[n,e,t]-1)
        #@constraint(m, accept_one[e=ext_edges,t=times], sum(x_n[n,e,t] for n=1:nb_seg) <= 1)
        @constraint(m, DefFp_n[n=1:nb_seg,e=ext_edges,t=times],
        #fp_n[n,e,t] == y_n[n,e,t]*price_qty[e,t][n][2][2] + (x_n[n,e,t]-y_n[n,e,t])*price_qty[e,t][n][2][1])
        #fp_n[n,e,t] == x_n[n,e,t]*price_qty[e,t][n][2][2] + (1-x_n[n,e,t])*price_qty[e,t][n][2][1])
        fp_n[n,e,t] == x_n[n,e,t]*(price_qty[e,t][n][2][2]-price_qty[e,t][n][2][1]))

        @constraint(m, def_fp_n[e=ext_edges,t=times], fp[e,t] == sum(fp_n[n,e,t] for n=1:nb_seg) + price_qty[e,t][1][2][1])

        ### BINARY CONSTRAINTS ###
        if binsol==0
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
        end

        ### OBJECTIVE ###
        @objective(m, Min,
        + sum(cost[sb][1]*pg[sb] for sb=tr_seg_bid)
        + cost_ls*sum(slackPUp[i,t] + slackPDown[i,t] for i=tr_nodes for t=times)
        + sum(fp_n[n,e,t]*price_qty[e,t][n][1] for n=1:nb_seg for e=ext_edges for t=times)
        #- sum(price_qty[e,t][n][2][1]*price_qty[e,t][n][1] for n=1:nb_seg for e=ext_edges for t=times)
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
        xHist=Dict{Tuple{Int64,Int64},Any}((e,t) => [(value(x_n[n,e,t]),price_qty[e,t][n][2]) for n=1:nb_seg] for e=ext_edges for t=times)
        prices_int=Dict{Tuple{Int64,Int64},Float64}()
        dc_sol=PFsol(false,pgHist,qgHist,cHist,sHist,theta_trHist,fpHist,fqHist,
        fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
        if binsol==0
        tr_bin_sol=Binsol(
		Dict{Any,Any}(sb=>value(sa[sb]) for sb=tr_seg_bid),
		Dict{Any,Any}(bd=>value(qa[bd]) for bd=tr_bid),
		Dict{Any,Any}(qt=>value(qta[qt]) for qt=tr_qt_bid),
		Dict{Any,Any}(sb=>value(alpha[sb]) for sb=AlphaOmegaIndex if sb[1] in tr_nodes),
		Dict{Any,Any}(sb=>value(omega[sb]) for sb=AlphaOmegaIndex if sb[1] in tr_nodes)
	)
        for e=ext_edges for t=times
                cur_n=nb_seg
                while (abs(xHist[e,t][cur_n][1]) < 1e-6 && cur_n > 1) cur_n-=1 end
                prices_int[e,t]=-price_qty[e,t][cur_n][1]
        end end
        return (dc_sol,tr_bin_sol,prices_int,xHist,solve_time(m))
        else
        for e=ext_edges for t=times
                #println(dual(def_fp_n[e,t]))
                prices_int[e,t]=dual(def_fp_n[e,t])
        end end
        end
        #(dc_sol_bis,dc_dual_sol,prices_int,sol_t)=solveDC_TN(dc_sol,tr_bin_sol)
        return (dc_sol,binsol,prices_int,xHist,solve_time(m))
end


function solve_TN_res_sup_fun_bin_no_last_step(nb_seg,price_qty,binsol,solver::Symbol=:Mosek)
        m=Model(optimizer_with_attributes(Mosek.Optimizer,
        "MSK_IPAR_PRESOLVE_USE" => 1,
        "MSK_IPAR_LOG" => 1,
        "MSK_IPAR_INTPNT_BASIS" => 0,
        "MSK_DPAR_MIO_TOL_ABS_RELAX_INT"=>1e-6,
        "MSK_DPAR_INTPNT_CO_TOL_NEAR_REL"=>100.0
        ))
        if solver==:Gurobi
                m=Model(optimizer_with_attributes(Gurobi.Optimizer,
                "BarHomogeneous" => 1,
                "Presolve" => -1,
                "QCPDual" => 1,
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
        @constraint(m, PDown[sb=tr_seg_bid], -pg[sb] <= -binsol.sa[sb]*p_min[sb])
        @constraint(m, PUp[sb=tr_seg_bid], pg[sb] <= binsol.sa[sb]*p_max[sb])

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

        @variable(m, fp_n[1:nb_seg,ext_edges,times])
        @variable(m, 0 <= x_n[1:nb_seg,ext_edges,times] <= 1)
        @constraint(m, DefFp_n[n=1:nb_seg,e=ext_edges,t=times],
        fp_n[n,e,t] == x_n[n,e,t]*(price_qty[e,t][n][2][2]-price_qty[e,t][n][2][1]))

        @constraint(m, def_fp_n[e=ext_edges,t=times], fp[e,t] == sum(fp_n[n,e,t] for n=1:nb_seg) + price_qty[e,t][1][2][1])

        ### OBJECTIVE ###
        @objective(m, Min,
        + sum(cost[sb][1]*pg[sb] for sb=tr_seg_bid)
        + cost_ls*sum(slackPUp[i,t] + slackPDown[i,t] for i=tr_nodes for t=times)
        + sum(fp_n[n,e,t]*price_qty[e,t][n][1] for n=1:nb_seg for e=ext_edges for t=times)
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
        xHist=Dict{Tuple{Int64,Int64},Any}((e,t) => [(value(x_n[n,e,t]),price_qty[e,t][n][2]) for n=1:nb_seg] for e=ext_edges for t=times)
        prices_int=Dict{Tuple{Int64,Int64},Float64}()
        dc_sol=PFsol(false,pgHist,qgHist,cHist,sHist,theta_trHist,fpHist,fqHist,
        fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
        priceP=Dict{Tuple{Int64,Int64},Float64}((i,t) => dual(PowerBalanceP[i,t]) for i=tr_nodes for t=times)
        priceQ=0
        dual_sol=Dualsol(priceP,priceQ)
        for e=ext_edges for t=times
                #println(dual(def_fp_n[e,t]))
                prices_int[e,t]=dual(def_fp_n[e,t])
        end end
        return (dc_sol,dual_sol,binsol,prices_int,xHist,solve_time(m))
end

function solve_DN_res_fun_bin(sol_TN::PFsol,binsol_TN::Binsol)
	solve_times=Dict{String,Float64}(dn=>0 for dn=keys(DNs))
        st_binsol=Dict{String,Float64}()
	println("Starting decomposition...")
	sol_DN=Dict{String,PFsol}()
	dualsol_DN=Dict{String,Dualsol}()
        binsol_DN=Dict{String,Binsol}()
	println("Starting $(length(root_dist_node)) Distribution networks problems.")
	for (num,dn) in enumerate(keys(DNs))
                (sol_DN[dn],tmp_sol_dual,binsol_DN[dn],st_binsol[dn])=solve_DN_relax(dn,sol_TN)
                if sol_DN[dn].isPFsol==false
		        (sol_DN[dn],dualsol_DN[dn],solve_times[dn])=solveAC_subnet(dn,sol_TN,binsol_DN[dn])
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
        binsol=Binsol(
		Dict{Any,Any}(sb=>binsol_TN.sa[sb] for sb=tr_seg_bid),
		Dict{Any,Any}(bd=>binsol_TN.qa[bd] for bd=tr_bid),
		Dict{Any,Any}(qt=>binsol_TN.qta[qt] for qt=tr_qt_bid),
		Dict{Any,Any}(sb=>binsol_TN.alpha[sb] for sb=AlphaOmegaIndex if sb[1] in tr_nodes),
		Dict{Any,Any}(sb=>binsol_TN.omega[sb] for sb=AlphaOmegaIndex if sb[1] in tr_nodes)
	)
        for dn=keys(DNs)
                for sb=seg_bid_DN[dn] binsol.sa[sb]=binsol_DN[dn].sa[sb] end
                for bd=bid_DN[dn] binsol.qa[bd]=binsol_DN[dn].qa[bd] end
                for qt=qt_bid_DN[dn] binsol.qta[qt]=binsol_DN[dn].qta[qt] end
                for sb=AlphaOmegaIndex if sb[1] in DNs[dn]
                        binsol.alpha[sb]=binsol_DN[dn].alpha[sb]
                        binsol.omega[sb]=binsol_DN[dn].omega[sb]
                end end
        end

	println("\n\n >>> Informations about the solution <<<")
	MaxSlack=maximum(abs,union(values(slackPUpHist),values(slackQUpHist),values(slackPDownHist), values(slackQDownHist)))
	println("Max Slack: $MaxSlack")
	pfsol=PFsol(true,pgHist,qgHist,cHist,sHist,thetaHist,
	fpHist,fqHist,fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
	ComputeGaps(pfsol)
	println("SW: $(SW(pfsol))")
        # println("\nTotal bin solve time: $(sum(values(st_binsol)))")
	# if !isempty(solve_times) println("\nTotal solve time: $(sum(values(solve_times)))") end
        # println("Max bin solve time: $(findmax(st_binsol))")
	# if !isempty(solve_times) println("Max solve time: $(findmax(solve_times))") end
        for dn=keys(DNs) st_binsol[dn]+=solve_times[dn] end
	return (pfsol,binsol,st_binsol)
end

function solve_DN_res_fun_bin_no_last_step(sol_TN::PFsol,binsol_TN::Binsol)
	solve_times=Dict{String,Float64}(dn=>0 for dn=keys(DNs))
        st_binsol=Dict{String,Float64}()
	println("Starting decomposition...")
	sol_DN=Dict{String,PFsol}()
	dualsol_DN=Dict{String,Dualsol}()
        binsol_DN=Dict{String,Binsol}()
	println("Starting $(length(root_dist_node)) Distribution networks problems.")
	for (num,dn) in enumerate(keys(DNs))
                (sol_DN[dn],tmp_sol_dual,binsol_DN[dn],st_binsol[dn])=solve_DN_relax(dn,sol_TN)
                #if sol_DN[dn].isPFsol==false
		(sol_DN[dn],dualsol_DN[dn],solve_times[dn])=solveAC_subnet(dn,sol_TN,binsol_DN[dn])
                #end
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
        binsol=Binsol(
		Dict{Any,Any}(sb=>binsol_TN.sa[sb] for sb=tr_seg_bid),
		Dict{Any,Any}(bd=>binsol_TN.qa[bd] for bd=tr_bid),
		Dict{Any,Any}(qt=>binsol_TN.qta[qt] for qt=tr_qt_bid),
		Dict{Any,Any}(sb=>binsol_TN.alpha[sb] for sb=AlphaOmegaIndex if sb[1] in tr_nodes),
		Dict{Any,Any}(sb=>binsol_TN.omega[sb] for sb=AlphaOmegaIndex if sb[1] in tr_nodes)
	)
        for dn=keys(DNs)
                for sb=seg_bid_DN[dn] binsol.sa[sb]=binsol_DN[dn].sa[sb] end
                for bd=bid_DN[dn] binsol.qa[bd]=binsol_DN[dn].qa[bd] end
                for qt=qt_bid_DN[dn] binsol.qta[qt]=binsol_DN[dn].qta[qt] end
                for sb=AlphaOmegaIndex if sb[1] in DNs[dn]
                        binsol.alpha[sb]=binsol_DN[dn].alpha[sb]
                        binsol.omega[sb]=binsol_DN[dn].omega[sb]
                end end
        end
        priceP=Dict()
        priceQ=Dict()
        for dn=keys(DNs)
                for i=DNs[dn] for t=times
                        priceP[i,t] = dualsol_DN[dn].priceP[i,t]
                        priceQ[i,t] = dualsol_DN[dn].priceQ[i,t]
                end end
        end

	println("\n\n >>> Informations about the solution <<<")
	MaxSlack=maximum(abs,union(values(slackPUpHist),values(slackQUpHist),values(slackPDownHist), values(slackQDownHist)))
	println("Max Slack: $MaxSlack")
	pfsol=PFsol(true,pgHist,qgHist,cHist,sHist,thetaHist,
	fpHist,fqHist,fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
	ComputeGaps(pfsol)
	println("SW: $(SW(pfsol))")
        for dn=keys(DNs) st_binsol[dn]+=solve_times[dn] end
        dual_sol=Dualsol(priceP,priceQ)
	return (pfsol,dual_sol,binsol,st_binsol)
end

function compute_decentralized_bin(nb_points,
                exp_flow::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}())
        solve_times=Dict{String,Any}()
        println("\n##############################\n### Step 1: Build RSF ###\n##############################\n")
	(price_qty,ST_points)=compute_n_points_bin(nb_points,exp_flow)
        price_qty=correct_price_qty(price_qty)
        solve_times["points"]=ST_points
	println("\n##########################\n### Step 2: Solve Transmission problem ###\n##########################\n")
	(sol_TN,bin_TN,prices_int,x_int,ST_TN)=solve_TN_res_sup_fun_bin(nb_points+1,price_qty)
        (sol_TN,bin_TN,prices_int,x_int,ST_TN_cont)=solve_TN_res_sup_fun_bin(nb_points+1,price_qty,bin_TN)
        solve_times["TN_RSF"]=ST_TN
        solve_times["TN_RSF_cont"]=ST_TN_cont
	println("\n################################\n### Step 3: Solve Distribution problems ###\n################################\n")
	(dec_sol,dec_bin_sol,dec_ST_pr)=solve_DN_res_fun_bin(sol_TN,bin_TN)
        solve_times["DN_RSF"]=dec_ST_pr
	println("\n##############################\n### Step 4: RSF Dual ###\n##############################\n")
	(dec_sol_bis,dec_dual_sol,dec_solve_times_du)=solveDCAC_bin_decomposed(dec_sol,dec_bin_sol,prices_int)
        solve_times["dual_RSF"]=dec_solve_times_du

        total_st=sum(sum(values(solve_times[i])) for i=keys(solve_times))
        max_par_st=sum(maximum(values(solve_times[i])) for i=keys(solve_times))

        println("\n>>> Total solve time: $total_st ($max_par_st)\n")

	return (dec_sol,dec_bin_sol,dec_dual_sol,prices_int,solve_times)
end

function enhance_compute_decentralized_bin(nb_points,
                exp_flow::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}())
        solve_times=Dict{String,Any}()
        println("\n##############################\n### Step 1.1: Build RSF ###\n##############################\n")
	(price_qty,ST_points)=compute_n_points_bin(round(Int64,nb_points/2),exp_flow)
        price_qty=correct_price_qty(price_qty)
        solve_times["points"]=ST_points
        println("\n##############################\n### Step 1.2: Enhance RSF ###\n##############################\n")
        (price_qty,ST_add_points)=add_points_price_qty(round(Int64,nb_points/2),price_qty)
        price_qty=correct_price_qty(price_qty)
        solve_times["add_points"]=ST_add_points
	println("\n##########################\n### Step 2: Solve Transmission problem ###\n##########################\n")
        (sol_TN,bin_TN,prices_int,x_int,ST_TN)=solve_TN_res_sup_fun_bin(nb_points+1,price_qty)
        (sol_TN,bin_TN,prices_int,x_int,ST_TN_cont)=solve_TN_res_sup_fun_bin(nb_points+1,price_qty,bin_TN)
        solve_times["TN_RSF"]=ST_TN
        solve_times["TN_RSF_cont"]=ST_TN_cont
	println("\n################################\n### Step 3: Solve Distribution problems ###\n################################\n")
        (dec_sol,dec_bin_sol,dec_ST_pr)=solve_DN_res_fun_bin(sol_TN,bin_TN)
        solve_times["DN_RSF"]=dec_ST_pr
	println("\n##############################\n### Step 4: RSF Dual ###\n##############################\n")
        (dec_sol_bis,dec_dual_sol,dec_solve_times_du)=solveDCAC_bin_decomposed(dec_sol,dec_bin_sol,prices_int)
        solve_times["dual_RSF"]=dec_solve_times_du

        total_st=sum(sum(values(solve_times[i])) for i=keys(solve_times))
        max_par_st=sum(maximum(values(solve_times[i])) for i=keys(solve_times))

        println("\n>>> Total solve time: $total_st ($max_par_st)\n")
	return (dec_sol,dec_bin_sol,dec_dual_sol,prices_int,solve_times)
end

function enhance_compute_no_rsf_last_step(nb_points,
                exp_flow::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}())
        solve_times=Dict{String,Any}()
        println("\n##############################\n### Step 1.1: Build RSF ###\n##############################\n")
	(price_qty,ST_points)=compute_n_points_bin(round(Int64,nb_points/2),exp_flow)
        price_qty=correct_price_qty(price_qty)
        solve_times["points"]=ST_points
        println("\n##############################\n### Step 1.2: Enhance RSF ###\n##############################\n")
        (price_qty,ST_add_points)=add_points_price_qty(round(Int64,nb_points/2),price_qty)
        price_qty=correct_price_qty(price_qty)
        solve_times["add_points"]=ST_add_points
	println("\n##########################\n### Step 2: Solve Transmission problem ###\n##########################\n")
        (sol_TN,bin_TN,prices_int,x_int,ST_TN)=solve_TN_res_sup_fun_bin(nb_points+1,price_qty)
        (sol_TN,dual_sol,bin_TN,prices_int,x_int,ST_TN_cont)=solve_TN_res_sup_fun_bin_no_last_step(nb_points+1,price_qty,bin_TN)
        solve_times["TN_RSF"]=ST_TN
        solve_times["TN_RSF_cont"]=ST_TN_cont
	println("\n################################\n### Step 3: Solve Distribution problems ###\n################################\n")
        (dec_sol,dist_dual_sol,dec_bin_sol,dec_ST_pr)=solve_DN_res_fun_bin_no_last_step(sol_TN,bin_TN)
        priceP=Dict()
        for i=tr_nodes for t=times priceP[i,t]=dual_sol.priceP[i,t] end end
        for i=dist_nodes for t=times priceP[i,t]=dist_dual_sol.priceP[i,t] end end
        dec_dual_sol=Dualsol(priceP,dist_dual_sol.priceQ)
        solve_times["DN_RSF"]=dec_ST_pr

        total_st=sum(sum(values(solve_times[i])) for i=keys(solve_times))
        max_par_st=sum(maximum(values(solve_times[i])) for i=keys(solve_times))

        println("\n>>> Total solve time: $total_st ($max_par_st)\n")
	return (dec_sol,dec_bin_sol,dec_dual_sol,prices_int,solve_times)
end

function enhance_compute_decentralized_exact_bin(nb_points,
                exp_flow::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}())
        solve_times=Dict{String,Any}()
        println("\n##############################\n### Step 1.1: Build RSF single t ###\n##############################\n")
	(rsf,ST_points)=compute_n_points_bin_approx_t(round(Int64,nb_points/2),exp_flow)
        rsf=correct_price_qty(rsf)
        solve_times["points"]=ST_points
        println("\n##############################\n### Step 1.2: Enhance RSF single t###\n##############################\n")
        (rsf,ST_add_points)=add_points_price_qty_approx_t(round(Int64,nb_points/2),rsf)
        rsf=correct_price_qty(rsf)
        solve_times["add_points"]=ST_add_points
        println("\n##############################\n### Step 1.3: Build RSF ###\n##############################\n")
	(price_qty,ST_points)=compute_n_points_bin_exact_t(round(Int64,nb_points/2),rsf,exp_flow)
        price_qty=correct_price_qty(price_qty)
        solve_times["points_exact"]=ST_points
        println("\n##############################\n### Step 1.4: Enhance RSF ###\n##############################\n")
        (price_qty,ST_add_points)=add_points_price_qty_exact_t(round(Int64,nb_points/2),rsf,price_qty)
        price_qty=correct_price_qty(price_qty)
        solve_times["add_points_exact"]=ST_add_points
        for i=keys(price_qty) for k=keys(price_qty[i])
                        price_qty[i][k]=(max(price_qty[i][k][1],rsf[i][k][1]),price_qty[i][k][2])
        end end
        # println("\n##############################\n### Step 1.1: Build RSF single t ###\n##############################\n")
	# (price_qty,ST_points)=compute_n_points_bin_approx_t(round(Int64,nb_points/2),exp_flow)
        # price_qty=correct_price_qty(price_qty)
        # solve_times["points"]=ST_points
        # println("\n##############################\n### Step 1.2: Enhance RSF single t###\n##############################\n")
        # (price_qty,ST_add_points)=add_points_price_qty_approx_t(round(Int64,nb_points/2),price_qty)
        # price_qty=correct_price_qty(price_qty)
        # solve_times["add_points"]=ST_add_points

	println("\n##########################\n### Step 2: Solve Transmission problem ###\n##########################\n")
        (sol_TN,bin_TN,prices_int,x_int,ST_TN)=solve_TN_res_sup_fun_bin(nb_points+1,price_qty)
        (sol_TN,bin_TN,prices_int,x_int,ST_TN_cont)=solve_TN_res_sup_fun_bin(nb_points+1,price_qty,bin_TN)
        solve_times["TN_RSF"]=ST_TN
        solve_times["TN_RSF_cont"]=ST_TN_cont
	println("\n################################\n### Step 3: Solve Distribution problems ###\n################################\n")
        (dec_sol,dec_bin_sol,dec_ST_pr)=solve_DN_res_fun_bin(sol_TN,bin_TN)
        solve_times["DN_RSF"]=dec_ST_pr
	println("\n##############################\n### Step 4: RSF Dual ###\n##############################\n")
        (dec_sol_bis,dec_dual_sol,dec_solve_times_du)=solveDCAC_bin_decomposed(dec_sol,dec_bin_sol,prices_int)
        solve_times["dual_RSF"]=dec_solve_times_du

        total_st=sum(sum(values(solve_times[i])) for i=keys(solve_times))
        max_par_st=sum(maximum(values(solve_times[i])) for i=keys(solve_times))

        println("\n>>> Total solve time: $total_st ($max_par_st)\n")
	return (dec_sol,dec_bin_sol,dec_dual_sol,prices_int,solve_times)
end

function enhance_compute_decentralized_oneill(nb_points,
        exp_flow::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}())
        solve_times=Dict{String,Any}()
        println("\n##############################\n### Step 1.1: Build RSF ###\n##############################\n")
	(price_qty,ST_points)=compute_n_points_bin_oneill(round(Int64,nb_points/2),exp_flow)
        price_qty=correct_price_qty(price_qty)
        solve_times["points"]=ST_points
        println("\n##############################\n### Step 1.2: Enhance RSF ###\n##############################\n")
        (price_qty,ST_add_points)=add_points_price_qty_oneill(round(Int64,nb_points/2),price_qty)
        price_qty=correct_price_qty(price_qty)
        solve_times["add_points"]=ST_add_points
	println("\n##########################\n### Step 2: Solve Transmission problem ###\n##########################\n")
        (sol_TN,bin_TN,prices_int,x_int,ST_TN)=solve_TN_res_sup_fun_bin(nb_points+1,price_qty)
        (sol_TN,bin_TN,prices_int,x_int,ST_TN_cont)=solve_TN_res_sup_fun_bin(nb_points+1,price_qty,bin_TN)
        solve_times["TN_RSF"]=ST_TN
        solve_times["TN_RSF_cont"]=ST_TN_cont
	println("\n################################\n### Step 3: Solve Distribution problems ###\n################################\n")
        (dec_sol,dec_bin_sol,dec_ST_pr)=solve_DN_res_fun_bin(sol_TN,bin_TN)
        solve_times["DN_RSF"]=dec_ST_pr
	println("\n##############################\n### Step 4: RSF Dual ###\n##############################\n")
        (dec_sol_bis,dec_dual_sol,dec_solve_times_du)=solveDCAC_bin_decomposed(dec_sol,dec_bin_sol,prices_int)
        solve_times["dual_RSF"]=dec_solve_times_du

        total_st=sum(sum(values(solve_times[i])) for i=keys(solve_times))
        max_par_st=sum(maximum(values(solve_times[i])) for i=keys(solve_times))

        println("\n>>> Total solve time: $total_st ($max_par_st)\n")
	return (dec_sol,dec_bin_sol,dec_dual_sol,prices_int,solve_times)
end

function representation_price_qty(x, single_price_qty)
        for i=1:length(single_price_qty)
                if x >= single_price_qty[i][2][1] && x < single_price_qty[i][2][2]
                        return single_price_qty[i][1]
                end
        end
        return 0
end
