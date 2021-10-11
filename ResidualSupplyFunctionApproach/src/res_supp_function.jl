function compute_res_fun(dn::String,
                         exp_flow::Dict{Tuple{Int64,Int64},Float64})
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

         m=Model(optimizer_with_attributes(Mosek.Optimizer,
 	"MSK_IPAR_PRESOLVE_USE" => 1,
        "MSK_IPAR_LOG" => 0,
 	))
         ### VARIABLES ###
         @variable(m, 0 <= c[dn_c_index,times])
         @variable(m, s[dn_s_index,times])
         @variable(m, pg[dn_seg_bid])
         @expression(m, net_pg_bid[bd=dn_bid], sum(pg[sb] for sb=bid_2_seg_bids[bd]))
         @variable(m, qg[dn_seg_bid])
         @expression(m, net_qg_bid[bd=dn_bid], sum(qg[sb] for sb=bid_2_seg_bids[bd]))
         @variable(m, -pi <= theta[union(dn_nodes,dn_ext_nodes),times] <= pi)

         @variable(m, abs(net_injP[i,t]) >= slackPUp[i=dn_nodes,t=times] >= 0)
         @variable(m, abs(net_injP[i,t]) >= slackPDown[i=dn_nodes,t=times] >= 0)
         @variable(m, abs(net_injQ[i,t]) >= slackQUp[i=dn_nodes,t=times] >= 0)
         @variable(m, abs(net_injQ[i,t]) >= slackQDown[i=dn_nodes,t=times] >= 0)
         # @variable(m, slackPUp[i=dn_nodes,t=times] >= 0)
         # @variable(m, slackPDown[i=dn_nodes,t=times] >= 0)
         # @variable(m, slackQUp[i=dn_nodes,t=times] >= 0)
         # @variable(m, slackQDown[i=dn_nodes,t=times] >= 0)

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
         # if binsol!=0
         #         @constraint(m, PDown[sb=dn_seg_bid], -pg[sb] <= -binsol.sa[sb]*p_min[sb])
         #         @constraint(m, PUp[sb=dn_seg_bid], pg[sb] <= binsol.sa[sb]*p_max[sb])
         # else
                 @variable(m, 0 <= relax_sa[dn_seg_bid] <= 1)
                 @constraint(m, PDown[sb=dn_seg_bid], -pg[sb] <= -relax_sa[sb]*p_min[sb])
                 @constraint(m, PUp[sb=dn_seg_bid], pg[sb] <= relax_sa[sb]*p_max[sb])
         #end

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
         # @constraint(m, FlowLimUp[e=dn_ext_edges,t=times],
         # fp[e,t] <= rate_a[e])
         # @constraint(m, FlowLimDown[e=dn_ext_edges,t=times],
         # -fp[e,t] <= rate_a[e])
         @constraint(m, DefFpX0[e=dn_ext_edges, t=times],
         fp[e,t] + fpTo[e,t] == 0)

         # ### OBJECTIVE ###
         @objective(m, Min,
         + sum(cost[sb][1]*pg[sb] for sb=dn_seg_bid)
         + cost_ls*sum(slackPUp[i,t] + slackPDown[i,t]
                         + slackQUp[i,t] + slackQDown[i,t] for i=dn_nodes for t=times)
         )

         ### Flow at the interconnection is fixed to a given value ###
        @constraint(m,fpFix[e=dn_ext_edges,t=times], fp[e,t]-exp_flow[e,t] == 0)

         optimize!(m)
         #@assert primal_status(m)==MOI.FEASIBLE_POINT
         #print(primal_status(m),"  ")
         ## only build a partial solution for this DN ##
         # pgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(pg[j]) for j=dn_seg_bid)
         # qgHist = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(j=> value(qg[j]) for j=dn_dist_seg_bid)
         # cHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((i,t) => value(c[i,t]) for i=dn_c_index for t=times)
         # sHist = Dict{Tuple{Tuple{Int64,Int64},Int64},Float64}((l,t) => value(s[l,t]) for l=dn_s_index for t=times)
         # thetaHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(theta[i,t]) for i=dn_nodes for t=times)
         # fpHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fp[l,t]) for l=dn_edges for t=times)
         # fqHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fq[l,t]) for l=dn_edges for t=times)
         # fpToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fpTo[l,t]) for l=dn_edges for t=times)
         # fqToHist = Dict{Tuple{Int64,Int64},Float64}((l,t) => value(fqTo[l,t]) for l=dn_edges for t=times)
         # priceP=Dict{Tuple{Int64,Int64},Float64}((i,t) => dual(PowerBalanceP[i,t])
         # for i=dn_nodes for t=times)
         # priceQ=Dict{Tuple{Int64,Int64},Float64}((i,t) => dual(PowerBalanceQ[i,t]) for i=intersect(dn_nodes,dist_nodes) for t=times)
         # slackPUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPUp[i,t])
         # for i=dn_nodes for t=times)
         # slackPDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackPDown[i,t])
         # for i=dn_nodes for t=times)
         # slackQUpHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackQUp[i,t])
         # for i=intersect(dn_nodes,dist_nodes) for t=times)
         # slackQDownHist = Dict{Tuple{Int64,Int64},Float64}((i,t) => value(slackQDown[i,t])
         # for i=intersect(dn_nodes,dist_nodes) for t=times)
         # MaxSlack=(maximum(abs,union(value.(slackPUp), value.(slackQUp),value.(slackPDown), value.(slackQDown))))
         # #println("Max Slack: $MaxSlack")
         # pfsol=PFsol(MaxSlack>1e-5,pgHist,qgHist,cHist,sHist,thetaHist,
         # fpHist,fqHist,fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
         # dual_sol=Dualsol(priceP,priceQ)
         qty_price_interface=Dict{Tuple{Int64,Int64},Tuple{Float64,Float64}}((e,t) => (dual(fpFix[e,t]),exp_flow[e,t]) for e=dn_ext_edges for t=times)
        if primal_status(m)!=MOI.FEASIBLE_POINT
                for e=dn_ext_edges for t=times
                        qty_price_interface[e,t]=(10*cost_ls,exp_flow[e,t])
                end end
        end
         return (qty_price_interface,solve_time(m))
end

function compute_n_points(nb_point::Int64,exp_flow::Dict{Tuple{Int64,Int64},Float64})
        solve_times=Dict{String,Float64}()
	println("Starting decomposition...")
	price_qty_DN=Dict{Tuple{Int64,Int64},Any}((e,t)=> [] for e=ext_edges for t=times)
	println("Starting $(length(root_dist_node)) Distribution networks problems.")
	for (num,dn) in enumerate(keys(DNs))
                st_array=[]
                inter_e=ext_edges_DN[dn]
                inter_e=first(inter_e)
                copy_exp_flow=copy(exp_flow)
                for n=0:(nb_point+1)
                        if inter_e in line_limit
                                qty=-rate_a[inter_e]+n*(2*rate_a[inter_e]/(nb_point+1))
                                for t=times exp_flow[inter_e,t]=qty end
                        else
                                frac=-2+4*n/(nb_point+1)
                                for t=times exp_flow[inter_e,t]=frac*copy_exp_flow[inter_e,t] end
                        end
                        (tmp_price_qty,tmp_st)=compute_res_fun(dn,exp_flow)
                        for i=keys(tmp_price_qty)
                                push!(price_qty_DN[i],tmp_price_qty[i])
                        end
                        push!(st_array,tmp_st)
                end
                solve_times[dn]=sum(st_array)
		println(">$num. SOLVED Distribution network problem: $dn.")
	end
        price_qty=transform_price_qty(nb_point+2,price_qty_DN)
	return (price_qty,solve_times)
end

function transform_price_qty(nb_point,price_qty)
        new_price_qty=Dict{Tuple{Int64,Int64},Any}((e,t)=>[] for e=ext_edges for t=times)
        # for e=ext_edges for t=times
        #         push!(new_price_qty[e,t],(price_qty[e,t][1][1],(-rate_a[e],price_qty[e,t][1][2])))
        # end end
        for e=ext_edges for t=times for n=2:(nb_point)
                push!(new_price_qty[e,t],(price_qty[e,t][n][1],(price_qty[e,t][n-1][2],price_qty[e,t][n][2])))
        end end end
        # for e=ext_edges for t=times
        #         push!(new_price_qty[e,t],(1000.0,(price_qty[e,t][nb_point][2],rate_a[e])))
        # end end
        return new_price_qty
end

function solve_TN_res_sup_fun(nb_seg,price_qty,solver::Symbol=:Mosek)
        m=Model(optimizer_with_attributes(Mosek.Optimizer,
        "MSK_IPAR_PRESOLVE_USE" => 1,
        "MSK_IPAR_LOG" => 1,
        "MSK_IPAR_INTPNT_BASIS" => 0,
        #"MSK_IPAR_PRESOLVE_MAX_NUM_REDUCTIONS"=>0,
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
        @variable(m, 0 <= relax_sa[tr_seg_bid] <= 1)
        @constraint(m, PDown[sb=tr_seg_bid], -pg[sb] <= -relax_sa[sb]*p_min[sb])
        @constraint(m, PUp[sb=tr_seg_bid], pg[sb] <= relax_sa[sb]*p_max[sb])

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
        fp_n[n,e,t] == sum(x_n[n,e,t]*price_qty[e,t][n][2][2] + (1-x_n[n,e,t])*price_qty[e,t][n][2][1]))

        #@constraint(m, at_most_one_acc[e=ext_edges,t=times], sum(acc_n[n,e,t] for n=1:nb_seg) <= 1)
        @constraint(m, def_fp_n[e=ext_edges,t=times], fp[e,t] == sum(fp_n[n,e,t] for n=1:nb_seg))
        #@constraint(m, seg_acc[n=2:nb_seg,e=ext_edges,t=times], x_n[n-1,e,t] >= x_n[n,e,t])
        # @constraint(m, seg_fp_nDown[n=1:nb_seg,e=ext_edges,t=times],
        # fp_n[n,e,t] >= acc_n[n,e,t]*price_qty[e,t][n][2][1])
        # @constraint(m, seg_fp_nUp[n=1:nb_seg,e=ext_edges,t=times],
        # fp_n[n,e,t] <= acc_n[n,e,t]*price_qty[e,t][n][2][2])

        ### OBJECTIVE ###
        @objective(m, Min,
        + sum(cost[sb][1]*pg[sb] for sb=tr_seg_bid)
        + cost_ls*sum(slackPUp[i,t] + slackPDown[i,t] for i=tr_nodes for t=times)
        + sum(fp_n[n,e,t]*price_qty[e,t][n][1] for n=1:nb_seg for e=ext_edges for t=times)
        - sum(price_qty[e,t][n][2][1]*price_qty[e,t][n][1] for n=1:nb_seg for e=ext_edges for t=times)
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
        for e=ext_edges for t=times
                cur_n=nb_seg
                while abs(xHist[e,t][cur_n][1]) < 1e-6 cur_n-=1 end
                prices_int[e,t]=price_qty[e,t][cur_n][1]
        end end
        dc_sol=PFsol(false,pgHist,qgHist,cHist,sHist,theta_trHist,fpHist,fqHist,
        fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
        return (dc_sol,prices_int,xHist,solve_time(m))
end

function solve_DN_res_fun(sol_TN::PFsol)
	solve_times=Dict{String,Float64}()
	println("Starting decomposition...")
	sol_DN=Dict{String,PFsol}()
	dualsol_DN=Dict{String,Dualsol}()
	println("Starting $(length(root_dist_node)) Distribution networks problems.")
	for (num,dn) in enumerate(keys(DNs))
		(sol_DN[dn],dualsol_DN[dn],solve_times[dn])=solveAC_subnet(dn,sol_TN)
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
	pfsol=PFsol(true,pgHist,qgHist,cHist,sHist,thetaHist,
	fpHist,fqHist,fpToHist,fqToHist,slackPUpHist,slackPDownHist,slackQUpHist,slackQDownHist)
	ComputeGaps(pfsol)
	println("SW: $(SW(pfsol))")
	println("\nTotal solve time: $(sum(values(solve_times)))")
	println("Max solve time: $(findmax(solve_times))")
	return (pfsol,solve_times)
end

function compute_decentralized(nb_points,
                exp_flow::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}())
        solve_times=Dict{String,Any}()
        println("\n##############################\n### Step 1: Build RSF ###\n##############################\n")
	(price_qty,ST_points)=compute_n_points(nb_points,exp_flow)
        solve_times["points"]=ST_points
	println("\n##########################\n### Step 2: Solve Transmission problem ###\n##########################\n")
	(sol_TN,prices_int,x_int,ST_TN)=solve_TN_res_sup_fun(nb_points+1,price_qty)
        solve_times["TN_RSF"]=ST_TN
	println("\n################################\n### Step 3: Solve Distribution problems ###\n################################\n")
	(dec_sol,dec_ST_pr)=solve_DN_res_fun(sol_TN)
        solve_times["DN_RSF"]=dec_ST_pr
	println("\n##############################\n### Step 4: RSF Dual ###\n##############################\n")
	(dec_sol_bis,dec_dual_sol,dec_solve_times_du)=solveDCAC_bin_decomposed(dec_sol,0,prices_int)
        solve_times["dual_RSF"]=dec_ST_pr
	return (dec_sol,dec_dual_sol,solve_times)
end

function successive_levels_rsf()
        sw_sol=[]
        for n=2:2:10
                (rsf_sol,rsf_dual,rsf_ST)=compute_decentralized(n)
                push!(sw_sol,SW(rsf_sol))
        end
        for n=15:5:40
                (rsf_sol,rsf_dual,rsf_ST)=compute_decentralized(n)
                push!(sw_sol,SW(rsf_sol))
        end
        return sw_sol
end


### EXTRA ###

function compute_alternative_n_points(nb_point::Int64,exp_flow::Dict{Tuple{Int64,Int64},Float64},cur_price_qty)
        if nb_point==0 return cur_price_qty end
        println("> Starting n=$(nb_point).")
        if isempty(cur_price_qty)
                (price_qty,tmp_st)=compute_n_points(1,exp_flow)
                return compute_alternative_n_points(nb_point-1,exp_flow,price_qty)
        end
        price_qty_DN=cur_price_qty
        copy_exp_flow=copy(exp_flow)
        for (num,dn) in enumerate(keys(DNs))
                st_array=[]
                inter_e=ext_edges_DN[dn]
                inter_e=first(inter_e)
                for t=times
                        price_qty_et=cur_price_qty[inter_e,t]
                        ind_seg=1
                        min_price_diff=0
                        for n=2:(length(price_qty_et))
                                if abs(price_qty_et[n][1]) > 1e3 && abs(price_qty_et[n-1][1]) > 1e3
                                        continue
                                end
                                if abs(price_qty_et[n][1]-price_qty_et[n-1][1]) > min_price_diff
                                        min_price_diff=abs(price_qty_et[n][1]-price_qty_et[n-1][1])
                                        ind_seg=n
                                end
                        end
                        copy_exp_flow[inter_e,t]=.5*(price_qty_et[ind_seg][2][1]+price_qty_et[ind_seg][2][2])
                end
                (tmp_price_qty,tmp_st)=compute_res_fun(dn,copy_exp_flow)
                for t=times
                        cur_n=1
                        while copy_exp_flow[inter_e,t] > price_qty_DN[inter_e,t][cur_n][2][2]
                                cur_n+=1
                        end
                        insert!(price_qty_DN[inter_e,t],cur_n,(tmp_price_qty[inter_e,t][1],(price_qty_DN[inter_e,t][cur_n-1][2][2],tmp_price_qty[inter_e,t][2])))
                        price_qty_DN[inter_e,t][cur_n+1]=(price_qty_DN[inter_e,t][cur_n+1][1],(tmp_price_qty[inter_e,t][2],price_qty_DN[inter_e,t][cur_n+1][2][2]))
                end
                println(">$num. SOLVED Distribution network problem: $dn.")
	end
	return compute_alternative_n_points(nb_point-1,exp_flow,price_qty_DN)
end


function compute_alternative_dec(nb_points,
                exp_flow::Dict{Tuple{Int64,Int64},Float64}=Dict{Tuple{Int64,Int64},Float64}())
        solve_times=Dict{String,Any}()
        cur_qty_price=Dict{Tuple{Int64,Int64},Any}()
        println("\n##############################\n### Step 1: Build RSF ###\n##############################\n")
	price_qty=compute_alternative_n_points(nb_points,exp_flow,cur_qty_price)
        #solve_times["points"]=ST_points
	println("\n##########################\n### Step 2: Solve Transmission problem ###\n##########################\n")
	(sol_TN,prices_int,x_int,ST_TN)=solve_TN_res_sup_fun(nb_points+1,price_qty)
        solve_times["TN_RSF"]=ST_TN
	println("\n################################\n### Step 3: Solve Distribution problems ###\n################################\n")
	(dec_sol,dec_ST_pr)=solve_DN_res_fun(sol_TN)
        solve_times["DN_RSF"]=dec_ST_pr
	println("\n##############################\n### Step 4: RSF Dual ###\n##############################\n")
	(dec_sol_bis,dec_dual_sol,dec_solve_times_du)=solveDCAC_bin_decomposed(dec_sol,0,prices_int)
        solve_times["dual_RSF"]=dec_ST_pr
	return (dec_sol,dec_dual_sol,solve_times)
end
