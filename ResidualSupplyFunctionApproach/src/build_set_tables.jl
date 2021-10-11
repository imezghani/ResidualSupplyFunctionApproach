
function build_set_tables(cur_time::Int64,sol::PFsol,dualsol::Dualsol,loco::Dict{Int64,Float64},locg::Dict{Tuple{String,Int64},Float64})::Tuple{DataFrame,DataFrame}
	actors=union(["TN"],keys(DNs),["MO"],["Total"])
	firstdf=DataFrame(
	actor=actors,
	real_pr=zeros(length(actors)),
	real_cp=zeros(length(actors)), # =>Dict{String,Float64}("TN" => sum(net_injP[i,cur_time]*dualsol.priceP[i,cur_time] for i=tr_nodes)),
	reactive_pr=zeros(length(actors)),
	reactive_cp=zeros(length(actors)),
	loc=zeros(length(actors)),
	total=zeros(length(actors))
	)
	firstdf.real_pr[firstdf.actor .== "TN"] .= (mapreduce(sum,+,init=0,max(0,sol.pg[sb])*dualsol.priceP[sb[1],sb[5]] for sb=tr_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time])*dualsol.priceP[i,cur_time] for i=tr_nodes))
	firstdf.real_cp[firstdf.actor .== "TN"] .= (mapreduce(sum,+,init=0,min(0,sol.pg[sb])*dualsol.priceP[sb[1],sb[5]] for sb=tr_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time])*dualsol.priceP[i,cur_time] for i=tr_nodes))
	for dn=keys(DNs)
		firstdf.real_pr[firstdf.actor .== dn] .= (mapreduce(sum,+,init=0,max(0,sol.pg[sb])*dualsol.priceP[sb[1],sb[5]] for sb=dist_seg_bid if sb[5]==cur_time && sb[1] in DNs[dn])
			+ mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time])*dualsol.priceP[i,cur_time] for i=DNs[dn]))
		firstdf.real_cp[firstdf.actor .== dn] .= (mapreduce(sum,+,init=0,min(0,sol.pg[sb])*dualsol.priceP[sb[1],sb[5]] for sb=dist_seg_bid if sb[5]==cur_time && sb[1] in DNs[dn])
			+ mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time])*dualsol.priceP[i,cur_time] for i=DNs[dn]))
		firstdf.reactive_pr[firstdf.actor .== dn] .= (mapreduce(sum,+,init=0,max(0,sol.qg[sb])*dualsol.priceQ[sb[1],sb[5]] for sb=dist_seg_bid if sb[5]==cur_time && sb[1] in DNs[dn])
			+ mapreduce(sum,+,init=0,max(0,net_injQ[i,cur_time])*dualsol.priceQ[i,cur_time] for i=DNs[dn]))
		firstdf.reactive_cp[firstdf.actor .== dn] .= (mapreduce(sum,+,init=0,min(0,sol.qg[sb])*dualsol.priceQ[sb[1],sb[5]] for sb=dist_seg_bid if sb[5]==cur_time && sb[1] in DNs[dn])
			+ mapreduce(sum,+,init=0,min(0,net_injQ[i,cur_time])*dualsol.priceQ[i,cur_time] for i=DNs[dn]))
	end

	#println(typeof(firstdf.real_pr[firstdf.actor .== "TN"]),"  ,  ",firstdf.real_pr[firstdf.actor .== "TN"])
	#println(typeof(sum(firstdf.real_pr[firstdf.actor .== "TN"])),"  ,  ",sum(firstdf.real_pr[firstdf.actor .== "TN"]))
	#println(typeof(firstdf.real_pr[firstdf.actor .== dn] for dn=keys(DNs)),firstdf.real_pr[firstdf.actor .== dn] for dn=keys(DNs))
	#println(typeof(sum(firstdf.real_pr[firstdf.actor .== dn] for dn=keys(DNs))),sum(firstdf.real_pr[firstdf.actor .== dn] for dn=keys(DNs)))
	firstdf.real_pr[firstdf.actor .== "MO"] .= -(firstdf.real_pr[firstdf.actor .== "TN"] + sum(firstdf.real_pr[firstdf.actor .== dn] for dn=keys(DNs)))
	firstdf.real_cp[firstdf.actor .== "MO"] .= -(firstdf.real_cp[firstdf.actor .== "TN"] + sum(firstdf.real_cp[firstdf.actor .== dn] for dn=keys(DNs)))
	firstdf.reactive_pr[firstdf.actor .== "MO"] .= -(firstdf.reactive_pr[firstdf.actor .== "TN"] + sum(firstdf.reactive_pr[firstdf.actor .== dn] for dn=keys(DNs)))
	firstdf.reactive_cp[firstdf.actor .== "MO"] .= -(firstdf.reactive_cp[firstdf.actor .== "TN"] + sum(firstdf.reactive_cp[firstdf.actor .== dn] for dn=keys(DNs)))

	firstdf.loc[firstdf.actor .== "TN"] .= locg[("TN",cur_time)]
	for dn=keys(DNs) firstdf.loc[firstdf.actor .== dn] .= locg[(dn,cur_time)] end
	firstdf.loc[firstdf.actor .== "MO"] .= loco[cur_time]

	firstdf.real_pr[firstdf.actor .== "Total"] .= sum(firstdf.real_pr[firstdf.actor .== act] for act=setdiff(actors,["Total"]))
	firstdf.real_cp[firstdf.actor .== "Total"] .= sum(firstdf.real_cp[firstdf.actor .== act] for act=setdiff(actors,["Total"]))
	firstdf.reactive_pr[firstdf.actor .== "Total"] .= sum(firstdf.real_pr[firstdf.actor .== act] for act=setdiff(actors,["Total"]))
	firstdf.reactive_cp[firstdf.actor .== "Total"] .= sum(firstdf.real_cp[firstdf.actor .== act] for act=setdiff(actors,["Total"]))
	firstdf.loc[firstdf.actor .== "Total"] .= sum(firstdf.loc[firstdf.actor .== act] for act=setdiff(actors,["Total"]))

	firstdf.total .= firstdf.real_pr+firstdf.real_cp+firstdf.reactive_pr+firstdf.reactive_cp+firstdf.loc


	seconddf=DataFrame(
	actor=actors,
	real_pc=zeros(length(actors)),
	real_cb=zeros(length(actors)), # =>Dict{String,Float64}("TN" => sum(net_injP[i,cur_time]*dualsol.priceP[i,cur_time] for i=tr_nodes)),
	reactive_pc=zeros(length(actors)),
	reactive_cb=zeros(length(actors)),
	total_rev=firstdf.total,
	profit=zeros(length(actors))
	)

	seconddf.real_pc[seconddf.actor .== "TN"] .= -(mapreduce(sum,+,init=0,max(0,sol.pg[sb])*cost[sb][1] for sb=tr_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time])*0 for i=tr_nodes))
	seconddf.real_cb[seconddf.actor .== "TN"] .= -(mapreduce(sum,+,init=0,min(0,sol.pg[sb])*cost[sb][1] for sb=tr_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time])*0 for i=tr_nodes))
	for dn=keys(DNs)
		seconddf.real_pc[seconddf.actor .== dn] .= -(mapreduce(sum,+,init=0,max(0,sol.pg[sb])*cost[sb][1] for sb=dist_seg_bid if sb[5]==cur_time && sb[1] in DNs[dn])
			+ mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time])*0 for i=DNs[dn]))
		seconddf.real_cb[seconddf.actor .== dn] .= -(mapreduce(sum,+,init=0,min(0,sol.pg[sb])*cost[sb][1] for sb=dist_seg_bid if sb[5]==cur_time && sb[1] in DNs[dn])
			+ mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time])*0 for i=DNs[dn]))
		seconddf.reactive_pc[seconddf.actor .== dn] .= -(mapreduce(sum,+,init=0,max(0,sol.qg[sb])*0 for sb=dist_seg_bid if sb[5]==cur_time && sb[1] in DNs[dn])
			+ mapreduce(sum,+,init=0,max(0,net_injQ[i,cur_time])*0 for i=DNs[dn]))
		seconddf.reactive_cb[seconddf.actor .== dn] .= -(mapreduce(sum,+,init=0,min(0,sol.qg[sb])*0 for sb=dist_seg_bid if sb[5]==cur_time && sb[1] in DNs[dn])
			+ mapreduce(sum,+,init=0,min(0,net_injQ[i,cur_time])*0 for i=DNs[dn]))
	end

	seconddf.real_pc[seconddf.actor .== "Total"] .= -(mapreduce(sum,+,init=0,max(0,sol.pg[sb])*cost[sb][1] for sb=seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time])*0 for i=nodes))
	seconddf.real_cb[seconddf.actor .== "Total"] .= -(mapreduce(sum,+,init=0,min(0,sol.pg[sb])*cost[sb][1] for sb=seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time])*0 for i=nodes))
	seconddf.reactive_pc[seconddf.actor .== "Total"] .= -(mapreduce(sum,+,init=0,max(0,sol.qg[sb])*0 for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injQ[i,cur_time])*0 for i=dist_nodes))
	seconddf.reactive_cb[seconddf.actor .== "Total"] .= -(mapreduce(sum,+,init=0,min(0,sol.qg[sb])*0 for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injQ[i,cur_time])*0 for i=dist_nodes))

	seconddf.profit .= seconddf.real_pc + seconddf.real_cb + seconddf.reactive_pc + seconddf.reactive_cb + seconddf.total_rev
	return (firstdf,seconddf)
end

function build_solution_table(cur_time::Int64,sol::PFsol,dualsol::Dualsol)::DataFrame
	actors=union(["TN"],keys(DNs),["Total"])
	finaldf=DataFrame(
	actor=actors,
	real_prod=zeros(length(actors)),
	real_cons=zeros(length(actors)),
	reactive_prod=zeros(length(actors)),
	reactive_cons=zeros(length(actors)),
	min_real_DLMP=zeros(length(actors)),
	max_real_DLMP=zeros(length(actors)),
	av_real_DLMP=zeros(length(actors)),
	min_reactive_DLMP=zeros(length(actors)),
	max_reactive_DLMP=zeros(length(actors)),
	av_reactive_DLMP=zeros(length(actors))
	)

	finaldf.real_prod[finaldf.actor .== "TN"] .= mapreduce(sum,+,init=0,max(0,sol.pg[sb]) for sb=tr_seg_bid if sb[5]==cur_time) + mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time]) for i=tr_nodes)
	finaldf.real_cons[finaldf.actor .== "TN"] .= -(mapreduce(sum,+,init=0,min(0,sol.pg[sb]) for sb=tr_seg_bid if sb[5]==cur_time) + mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time]) for i=tr_nodes))
	for dn=keys(DNs)
		finaldf.real_prod[finaldf.actor .== dn] .= mapreduce(sum,+,init=0,max(0,sol.pg[sb]) for sb=dist_seg_bid if sb[5]==cur_time && sb[1] in DNs[dn]) + mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time]) for i=DNs[dn])
		finaldf.real_cons[finaldf.actor .== dn] .= -(mapreduce(sum,+,init=0,min(0,sol.pg[sb]) for sb=dist_seg_bid if sb[5]==cur_time && sb[1] in DNs[dn]) + mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time]) for i=DNs[dn]))
		finaldf.reactive_prod[finaldf.actor .== dn] .= mapreduce(sum,+,init=0,max(0,sol.qg[sb]) for sb=dist_seg_bid if sb[5]==cur_time && sb[1] in DNs[dn]) + mapreduce(sum,+,init=0,max(0,net_injQ[i,cur_time]) for i=DNs[dn])
		finaldf.reactive_cons[finaldf.actor .== dn] .= -(mapreduce(sum,+,init=0,min(0,sol.qg[sb]) for sb=dist_seg_bid if sb[5]==cur_time && sb[1] in DNs[dn]) + mapreduce(sum,+,init=0,min(0,net_injQ[i,cur_time]) for i=DNs[dn]))
	end

	finaldf.min_real_DLMP[finaldf.actor .== "TN"] .= minimum(dualsol.priceP[i,cur_time] for i=tr_nodes)
	finaldf.max_real_DLMP[finaldf.actor .== "TN"] .= maximum(dualsol.priceP[i,cur_time] for i=tr_nodes)
	finaldf.av_real_DLMP[finaldf.actor .== "TN"] .= mapreduce(sum,+,init=0,dualsol.priceP[i,cur_time] for i=tr_nodes)/length(tr_nodes)
	for dn=keys(DNs)
		finaldf.min_real_DLMP[finaldf.actor .== dn] .= minimum(dualsol.priceP[i,cur_time] for i=DNs[dn])
		finaldf.max_real_DLMP[finaldf.actor .== dn] .= maximum(dualsol.priceP[i,cur_time] for i=DNs[dn])
		finaldf.av_real_DLMP[finaldf.actor .== dn] .= mapreduce(sum,+,init=0,dualsol.priceP[i,cur_time] for i=DNs[dn])/length(DNs[dn])
		finaldf.min_reactive_DLMP[finaldf.actor .== dn] .= minimum(dualsol.priceQ[i,cur_time] for i=DNs[dn])
		finaldf.max_reactive_DLMP[finaldf.actor .== dn] .= maximum(dualsol.priceQ[i,cur_time] for i=DNs[dn])
		finaldf.av_reactive_DLMP[finaldf.actor .== dn] .= mapreduce(sum,+,init=0,dualsol.priceQ[i,cur_time] for i=DNs[dn])/length(DNs[dn])
	end

	finaldf.real_prod[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,max(0,sol.pg[sb]) for sb=seg_bid if sb[5]==cur_time) + mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time]) for i=nodes)
	finaldf.real_cons[finaldf.actor .== "Total"] .= -(mapreduce(sum,+,init=0,min(0,sol.pg[sb]) for sb=seg_bid if sb[5]==cur_time) + mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time]) for i=nodes))
	finaldf.reactive_prod[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,max(0,sol.qg[sb]) for sb=dist_seg_bid if sb[5]==cur_time) + mapreduce(sum,+,init=0,max(0,net_injQ[i,cur_time]) for i=dist_nodes)
	finaldf.reactive_cons[finaldf.actor .== "Total"] .= -(mapreduce(sum,+,init=0,min(0,sol.qg[sb]) for sb=dist_seg_bid if sb[5]==cur_time) + mapreduce(sum,+,init=0,min(0,net_injQ[i,cur_time]) for i=dist_nodes))
	finaldf.min_real_DLMP[finaldf.actor .== "Total"] .= minimum(dualsol.priceP[i,cur_time] for i=nodes)
	finaldf.max_real_DLMP[finaldf.actor .== "Total"] .= maximum(dualsol.priceP[i,cur_time] for i=nodes)
	finaldf.av_real_DLMP[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,dualsol.priceP[i,cur_time] for i=nodes)/length(nodes)
	finaldf.min_reactive_DLMP[finaldf.actor .== "Total"] .= minimum(dualsol.priceQ[i,cur_time] for i=dist_nodes)
	finaldf.max_reactive_DLMP[finaldf.actor .== "Total"] .= maximum(dualsol.priceQ[i,cur_time] for i=dist_nodes)
	finaldf.av_reactive_DLMP[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,dualsol.priceQ[i,cur_time] for i=dist_nodes)/length(dist_nodes)
	return finaldf
end

function set_tables_to_CSV(dir::String,sol::PFsol,dualsol::Dualsol,loco::Dict{Int64,Float64},locg::Dict{Tuple{String,Int64},Float64})
	first_table=Dict{Int64,DataFrame}(t=>DataFrame() for t=times)
	second_table=Dict{Int64,DataFrame}(t=>DataFrame() for t=times)
	solution_table=Dict{Int64,DataFrame}(t=>DataFrame() for t=times)
	for t=times
		(first_table[t],second_table[t]) = build_set_tables(t,sol,dualsol,loco,locg)
	end
	for t=times
		solution_table[t] = build_solution_table(t,sol,dualsol)
	end
	if !(dir in readdir()) mkdir(dir) end
	for t=times
		CSV.write("$dir/set_table1_T$t.csv",first_table[t])
		CSV.write("$dir/set_table2_T$t.csv",second_table[t])
		CSV.write("$dir/solution_table_T$t.csv",solution_table[t])
	end
	total_first_table = DataFrame(
	actor=first_table[first(times)].actor,
	real_pr=sum(first_table[t].real_pr for t=times),
	real_cp=sum(first_table[t].real_cp for t=times),
	reactive_pr=sum(first_table[t].reactive_pr for t=times),
	reactive_cp=sum(first_table[t].reactive_cp for t=times),
	loc=sum(first_table[t].loc for t=times),
	total=sum(first_table[t].total for t=times)
	)
	total_second_table = DataFrame(
	actor=second_table[first(times)].actor,
	real_pc=sum(second_table[t].real_pc for t=times),
	real_cb=sum(second_table[t].real_cb for t=times),
	reactive_pc=sum(second_table[t].reactive_pc for t=times),
	reactive_cb=sum(second_table[t].reactive_cb for t=times),
	total_rev=sum(second_table[t].total_rev for t=times),
	profit=sum(second_table[t].profit for t=times)
	)
	total_solution_table = DataFrame(
	actor=solution_table[first(times)].actor,
	real_prod=sum(solution_table[t].real_prod for t=times),
	real_cons=sum(solution_table[t].real_cons for t=times),
	reactive_prod=sum(solution_table[t].reactive_prod for t=times),
	reactive_cons=sum(solution_table[t].reactive_cons for t=times),
	min_real_DLMP=[minimum(solution_table[t].min_real_DLMP[i] for t=times) for i=1:length(solution_table[first(times)].actor)],
	max_real_DLMP=[maximum(solution_table[t].max_real_DLMP[i] for t=times) for i=1:length(solution_table[first(times)].actor)],
	av_real_DLMP=sum(solution_table[t].av_real_DLMP for t=times)./length(times),
	min_reactive_DLMP=[minimum(solution_table[t].min_reactive_DLMP[i] for t=times) for i=1:length(solution_table[first(times)].actor)],
	max_reactive_DLMP=[maximum(solution_table[t].max_reactive_DLMP[i] for t=times) for i=1:length(solution_table[first(times)].actor)],
	av_reactive_DLMP=sum(solution_table[t].av_reactive_DLMP for t=times)./length(times)
	)
	CSV.write("$dir/total_set_table1.csv",total_first_table)
	CSV.write("$dir/total_set_table2.csv",total_second_table)
	CSV.write("$dir/total_solution_table.csv",total_solution_table)
end

function build_aggr_set_tables(cur_time::Int64,sol::PFsol,dualsol::Dualsol,loco::Dict{Int64,Float64},locg::Dict{Tuple{String,Int64},Float64})::Tuple{DataFrame,DataFrame}
	actors=union(["TN"],["DNs"],["MO"],["Total"])
	firstdf=DataFrame(
	actor=actors,
	real_pr=zeros(length(actors)),
	real_cp=zeros(length(actors)), # =>Dict{String,Float64}("TN" => sum(net_injP[i,cur_time]*dualsol.priceP[i,cur_time] for i=tr_nodes)),
	reactive_pr=zeros(length(actors)),
	reactive_cp=zeros(length(actors)),
	loc=zeros(length(actors)),
	total=zeros(length(actors))
	)
	firstdf.real_pr[firstdf.actor .== "TN"] .= (mapreduce(sum,+,init=0,max(0,sol.pg[sb])*dualsol.priceP[sb[1],sb[5]] for sb=tr_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time])*dualsol.priceP[i,cur_time] for i=tr_nodes))
	firstdf.real_cp[firstdf.actor .== "TN"] .= (mapreduce(sum,+,init=0,min(0,sol.pg[sb])*dualsol.priceP[sb[1],sb[5]] for sb=tr_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time])*dualsol.priceP[i,cur_time] for i=tr_nodes))
	firstdf.real_pr[firstdf.actor .== "DNs"] .= (mapreduce(sum,+,init=0,max(0,sol.pg[sb])*dualsol.priceP[sb[1],sb[5]] for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time])*dualsol.priceP[i,cur_time] for i=dist_nodes))
	firstdf.real_cp[firstdf.actor .== "DNs"] .= (mapreduce(sum,+,init=0,min(0,sol.pg[sb])*dualsol.priceP[sb[1],sb[5]] for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time])*dualsol.priceP[i,cur_time] for i=dist_nodes))
	firstdf.reactive_pr[firstdf.actor .== "DNs"] .= (mapreduce(sum,+,init=0,max(0,sol.qg[sb])*dualsol.priceQ[sb[1],sb[5]] for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injQ[i,cur_time])*dualsol.priceQ[i,cur_time] for i=dist_nodes))
	firstdf.reactive_cp[firstdf.actor .== "DNs"] .= (mapreduce(sum,+,init=0,min(0,sol.qg[sb])*dualsol.priceQ[sb[1],sb[5]] for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injQ[i,cur_time])*dualsol.priceQ[i,cur_time] for i=dist_nodes))

	#println(typeof(firstdf.real_pr[firstdf.actor .== "TN"]),"  ,  ",firstdf.real_pr[firstdf.actor .== "TN"])
	#println(typeof(sum(firstdf.real_pr[firstdf.actor .== "TN"])),"  ,  ",sum(firstdf.real_pr[firstdf.actor .== "TN"]))
	#println(typeof(firstdf.real_pr[firstdf.actor .== dn] for dn=keys(DNs)),firstdf.real_pr[firstdf.actor .== dn] for dn=keys(DNs))
	#println(typeof(sum(firstdf.real_pr[firstdf.actor .== dn] for dn=keys(DNs))),sum(firstdf.real_pr[firstdf.actor .== dn] for dn=keys(DNs)))
	firstdf.real_pr[firstdf.actor .== "MO"] .= -(firstdf.real_pr[firstdf.actor .== "TN"] + (firstdf.real_pr[firstdf.actor .== "DNs"]))
	firstdf.real_cp[firstdf.actor .== "MO"] .= -(firstdf.real_cp[firstdf.actor .== "TN"] + (firstdf.real_cp[firstdf.actor .== "DNs"]))
	firstdf.reactive_pr[firstdf.actor .== "MO"] .= -(firstdf.reactive_pr[firstdf.actor .== "TN"] + (firstdf.reactive_pr[firstdf.actor .== "DNs"]))
	firstdf.reactive_cp[firstdf.actor .== "MO"] .= -(firstdf.reactive_cp[firstdf.actor .== "TN"] + (firstdf.reactive_cp[firstdf.actor .== "DNs"]))

	firstdf.loc[firstdf.actor .== "TN"] .= locg[("TN",cur_time)]
	firstdf.loc[firstdf.actor .== "DNs"] .= sum(locg[(dn,cur_time)] for dn=keys(DNs))
	firstdf.loc[firstdf.actor .== "MO"] .= loco[cur_time]

	firstdf.real_pr[firstdf.actor .== "Total"] .= sum(firstdf.real_pr[firstdf.actor .== act] for act=setdiff(actors,["Total"]))
	firstdf.real_cp[firstdf.actor .== "Total"] .= sum(firstdf.real_cp[firstdf.actor .== act] for act=setdiff(actors,["Total"]))
	firstdf.reactive_pr[firstdf.actor .== "Total"] .= sum(firstdf.real_pr[firstdf.actor .== act] for act=setdiff(actors,["Total"]))
	firstdf.reactive_cp[firstdf.actor .== "Total"] .= sum(firstdf.real_cp[firstdf.actor .== act] for act=setdiff(actors,["Total"]))
	firstdf.loc[firstdf.actor .== "Total"] .= sum(firstdf.loc[firstdf.actor .== act] for act=setdiff(actors,["Total"]))

	firstdf.total .= firstdf.real_pr+firstdf.real_cp+firstdf.reactive_pr+firstdf.reactive_cp+firstdf.loc


	seconddf=DataFrame(
	actor=actors,
	real_pc=zeros(length(actors)),
	real_cb=zeros(length(actors)), # =>Dict{String,Float64}("TN" => sum(net_injP[i,cur_time]*dualsol.priceP[i,cur_time] for i=tr_nodes)),
	reactive_pc=zeros(length(actors)),
	reactive_cb=zeros(length(actors)),
	total_rev=firstdf.total,
	profit=zeros(length(actors))
	)

	seconddf.real_pc[seconddf.actor .== "TN"] .= -(mapreduce(sum,+,init=0,max(0,sol.pg[sb])*cost[sb][1] for sb=tr_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time])*0 for i=tr_nodes))
	seconddf.real_cb[seconddf.actor .== "TN"] .= -(mapreduce(sum,+,init=0,min(0,sol.pg[sb])*cost[sb][1] for sb=tr_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time])*0 for i=tr_nodes))
	seconddf.real_pc[seconddf.actor .== "DNs"] .= -(mapreduce(sum,+,init=0,max(0,sol.pg[sb])*cost[sb][1] for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time])*0 for i=dist_nodes))
	seconddf.real_cb[seconddf.actor .== "DNs"] .= -(mapreduce(sum,+,init=0,min(0,sol.pg[sb])*cost[sb][1] for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time])*0 for i=dist_nodes))
	seconddf.reactive_pc[seconddf.actor .== "DNs"] .= -(mapreduce(sum,+,init=0,max(0,sol.qg[sb])*0 for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injQ[i,cur_time])*0 for i=dist_nodes))
	seconddf.reactive_cb[seconddf.actor .== "DNs"] .= -(mapreduce(sum,+,init=0,min(0,sol.qg[sb])*0 for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injQ[i,cur_time])*0 for i=dist_nodes))

	seconddf.real_pc[seconddf.actor .== "Total"] .= -(mapreduce(sum,+,init=0,max(0,sol.pg[sb])*cost[sb][1] for sb=seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time])*0 for i=nodes))
	seconddf.real_cb[seconddf.actor .== "Total"] .= -(mapreduce(sum,+,init=0,min(0,sol.pg[sb])*cost[sb][1] for sb=seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time])*0 for i=nodes))
	seconddf.reactive_pc[seconddf.actor .== "Total"] .= -(mapreduce(sum,+,init=0,max(0,sol.qg[sb])*0 for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injQ[i,cur_time])*0 for i=dist_nodes))
	seconddf.reactive_cb[seconddf.actor .== "Total"] .= -(mapreduce(sum,+,init=0,min(0,sol.qg[sb])*0 for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injQ[i,cur_time])*0 for i=dist_nodes))

	seconddf.profit .= seconddf.real_pc + seconddf.real_cb + seconddf.reactive_pc + seconddf.reactive_cb + seconddf.total_rev
	return (firstdf,seconddf)
end

function build_aggr_solution_table(cur_time::Int64,sol::PFsol,dualsol::Dualsol)::DataFrame
	actors=union(["TN"],["DNs"],["Total"])
	finaldf=DataFrame(
	actor=actors,
	real_prod=zeros(length(actors)),
	real_cons=zeros(length(actors)),
	reactive_prod=zeros(length(actors)),
	reactive_cons=zeros(length(actors)),
	min_real_DLMP=zeros(length(actors)),
	max_real_DLMP=zeros(length(actors)),
	av_real_DLMP=zeros(length(actors)),
	min_reactive_DLMP=zeros(length(actors)),
	max_reactive_DLMP=zeros(length(actors)),
	av_reactive_DLMP=zeros(length(actors))
	)

	finaldf.real_prod[finaldf.actor .== "TN"] .= (mapreduce(sum,+,init=0,max(0,sol.pg[sb]) for sb=tr_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time]) for i=tr_nodes))
	finaldf.real_cons[finaldf.actor .== "TN"] .= (-(mapreduce(sum,+,init=0,min(0,sol.pg[sb]) for sb=tr_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time]) for i=tr_nodes)))
	finaldf.real_prod[finaldf.actor .== "DNs"] .= (mapreduce(sum,+,init=0,max(0,sol.pg[sb]) for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time]) for i=dist_nodes))
	finaldf.real_cons[finaldf.actor .== "DNs"] .= (-(mapreduce(sum,+,init=0,min(0,sol.pg[sb]) for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time]) for i=dist_nodes)))
	finaldf.reactive_prod[finaldf.actor .== "DNs"] .= (mapreduce(sum,+,init=0,max(0,sol.qg[sb]) for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,max(0,net_injQ[i,cur_time]) for i=dist_nodes))
	finaldf.reactive_cons[finaldf.actor .== "DNs"] .= (-(mapreduce(sum,+,init=0,min(0,sol.qg[sb]) for sb=dist_seg_bid if sb[5]==cur_time)
		+ mapreduce(sum,+,init=0,min(0,net_injQ[i,cur_time]) for i=dist_nodes)))

	finaldf.min_real_DLMP[finaldf.actor .== "TN"] .= minimum(dualsol.priceP[i,cur_time] for i=tr_nodes)
	finaldf.max_real_DLMP[finaldf.actor .== "TN"] .= maximum(dualsol.priceP[i,cur_time] for i=tr_nodes)
	finaldf.av_real_DLMP[finaldf.actor .== "TN"] .= mapreduce(sum,+,init=0,dualsol.priceP[i,cur_time] for i=tr_nodes)/length(tr_nodes)
	finaldf.min_real_DLMP[finaldf.actor .== "DNs"] .= minimum(dualsol.priceP[i,cur_time] for i=dist_nodes)
	finaldf.max_real_DLMP[finaldf.actor .== "DNs"] .= maximum(dualsol.priceP[i,cur_time] for i=dist_nodes)
	finaldf.av_real_DLMP[finaldf.actor .== "DNs"] .= mapreduce(sum,+,init=0,dualsol.priceP[i,cur_time] for i=dist_nodes)/length(dist_nodes)
	finaldf.min_reactive_DLMP[finaldf.actor .== "DNs"] .= minimum(dualsol.priceQ[i,cur_time] for i=dist_nodes)
	finaldf.max_reactive_DLMP[finaldf.actor .== "DNs"] .= maximum(dualsol.priceQ[i,cur_time] for i=dist_nodes)
	finaldf.av_reactive_DLMP[finaldf.actor .== "DNs"] .= mapreduce(sum,+,init=0,dualsol.priceQ[i,cur_time] for i=dist_nodes)/length(dist_nodes)

	finaldf.real_prod[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,max(0,sol.pg[sb]) for sb=seg_bid if sb[5]==cur_time) + mapreduce(sum,+,init=0,max(0,net_injP[i,cur_time]) for i=nodes)
	finaldf.real_cons[finaldf.actor .== "Total"] .= -(mapreduce(sum,+,init=0,min(0,sol.pg[sb]) for sb=seg_bid if sb[5]==cur_time) + mapreduce(sum,+,init=0,min(0,net_injP[i,cur_time]) for i=nodes))
	finaldf.reactive_prod[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,max(0,sol.qg[sb]) for sb=dist_seg_bid if sb[5]==cur_time) + mapreduce(sum,+,init=0,max(0,net_injQ[i,cur_time]) for i=dist_nodes)
	finaldf.reactive_cons[finaldf.actor .== "Total"] .= -(mapreduce(sum,+,init=0,min(0,sol.qg[sb]) for sb=dist_seg_bid if sb[5]==cur_time) + mapreduce(sum,+,init=0,min(0,net_injQ[i,cur_time]) for i=dist_nodes))
	finaldf.min_real_DLMP[finaldf.actor .== "Total"] .= minimum(dualsol.priceP[i,cur_time] for i=nodes)
	finaldf.max_real_DLMP[finaldf.actor .== "Total"] .= maximum(dualsol.priceP[i,cur_time] for i=nodes)
	finaldf.av_real_DLMP[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,dualsol.priceP[i,cur_time] for i=nodes)/length(nodes)
	finaldf.min_reactive_DLMP[finaldf.actor .== "Total"] .= minimum(dualsol.priceQ[i,cur_time] for i=dist_nodes)
	finaldf.max_reactive_DLMP[finaldf.actor .== "Total"] .= maximum(dualsol.priceQ[i,cur_time] for i=dist_nodes)
	finaldf.av_reactive_DLMP[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,dualsol.priceQ[i,cur_time] for i=dist_nodes)/length(dist_nodes)
	return finaldf
end

function set_aggr_tables_to_CSV(dir::String,sol::PFsol,dualsol::Dualsol,loco::Dict{Int64,Float64},locg::Dict{Tuple{String,Int64},Float64})
	first_table=Dict{Int64,DataFrame}(t=>DataFrame() for t=times)
	second_table=Dict{Int64,DataFrame}(t=>DataFrame() for t=times)
	solution_table=Dict{Int64,DataFrame}(t=>DataFrame() for t=times)
	for t=times
		(first_table[t],second_table[t]) = build_aggr_set_tables(t,sol,dualsol,loco,locg)
	end
	for t=times
		solution_table[t] = build_aggr_solution_table(t,sol,dualsol)
	end
	if !(dir in readdir()) mkdir(dir) end
	for t=times
		CSV.write("$dir/set_table1_T$t.csv",first_table[t])
		CSV.write("$dir/set_table2_T$t.csv",second_table[t])
		CSV.write("$dir/solution_table_T$t.csv",solution_table[t])
	end
	total_first_table = DataFrame(
	actor=first_table[first(times)].actor,
	real_pr=sum(first_table[t].real_pr for t=times),
	real_cp=sum(first_table[t].real_cp for t=times),
	reactive_pr=sum(first_table[t].reactive_pr for t=times),
	reactive_cp=sum(first_table[t].reactive_cp for t=times),
	loc=sum(first_table[t].loc for t=times),
	total=sum(first_table[t].total for t=times)
	)
	total_second_table = DataFrame(
	actor=second_table[first(times)].actor,
	real_pc=sum(second_table[t].real_pc for t=times),
	real_cb=sum(second_table[t].real_cb for t=times),
	reactive_pc=sum(second_table[t].reactive_pc for t=times),
	reactive_cb=sum(second_table[t].reactive_cb for t=times),
	total_rev=sum(second_table[t].total_rev for t=times),
	profit=sum(second_table[t].profit for t=times)
	)
	total_solution_table = DataFrame(
	actor=solution_table[first(times)].actor,
	real_prod=sum(solution_table[t].real_prod for t=times),
	real_cons=sum(solution_table[t].real_cons for t=times),
	reactive_prod=sum(solution_table[t].reactive_prod for t=times),
	reactive_cons=sum(solution_table[t].reactive_cons for t=times),
	min_real_DLMP=[minimum(solution_table[t].min_real_DLMP[i] for t=times) for i=1:length(solution_table[first(times)].actor)],
	max_real_DLMP=[maximum(solution_table[t].max_real_DLMP[i] for t=times) for i=1:length(solution_table[first(times)].actor)],
	av_real_DLMP=sum(solution_table[t].av_real_DLMP for t=times)./length(times),
	min_reactive_DLMP=[minimum(solution_table[t].min_reactive_DLMP[i] for t=times) for i=1:length(solution_table[first(times)].actor)],
	max_reactive_DLMP=[maximum(solution_table[t].max_reactive_DLMP[i] for t=times) for i=1:length(solution_table[first(times)].actor)],
	av_reactive_DLMP=sum(solution_table[t].av_reactive_DLMP for t=times)./length(times)
	)
	CSV.write("$dir/total_set_table1.csv",total_first_table)
	CSV.write("$dir/total_set_table2.csv",total_second_table)
	CSV.write("$dir/total_solution_table.csv",total_solution_table)
end

function build_ads_set_tables(cur_time::Int64,sol::PFsol,dualsol::Dualsol,prices_int=0)::DataFrame
	actors=union(["BSPs Tr"],["BRPs Tr"],["TSO"],["ADS"],["BSPs Dist"],["BRPs Dist"],["Total"])
	firstdf=DataFrame(
	actor=actors,
	tm_bsp_set=zeros(length(actors)),
	tm_brp_set=zeros(length(actors)),
	ads_disaggregate=zeros(length(actors)),
	total=zeros(length(actors)),
	)
	firstdf.tm_bsp_set[firstdf.actor .== "BSPs Tr"] .= mapreduce(sum,+,init=0,sol.pg[sb]*dualsol.priceP[sb[1],sb[5]] for sb=tr_seg_bid if sb[5]==cur_time)
	firstdf.tm_bsp_set[firstdf.actor .== "BRPs Tr"] .= 0.0
	firstdf.tm_bsp_set[firstdf.actor .== "TSO"] .= -mapreduce(sum,+,init=0,sol.pg[sb]*dualsol.priceP[sb[1],sb[5]] for sb=tr_seg_bid if sb[5]==cur_time)
	if prices_int!=0
		#firstdf.tm_bsp_set[firstdf.actor .== "TSO"] .+=	mapreduce(sum,+,init=0,sol.fp[e,cur_time]*prices_int[e,cur_time] for e=ext_edges)
		#firstdf.tm_bsp_set[firstdf.actor .== "ADS"] .= -mapreduce(sum,+,init=0,sol.fp[e,cur_time]*prices_int[e,cur_time] for e=ext_edges)
		firstdf.tm_bsp_set[firstdf.actor .== "TSO"] .+=	mapreduce(sum,+,init=0,mapreduce(sum,+,init=0,sol.pg[sb] for sb=seg_bid_DN[dn] if sb[5]==cur_time)*prices_int[first(ext_edges_DN[dn]),cur_time] for dn=keys(DNs))
		firstdf.tm_bsp_set[firstdf.actor .== "ADS"] .= -mapreduce(sum,+,init=0,mapreduce(sum,+,init=0,sol.pg[sb] for sb=seg_bid_DN[dn] if sb[5]==cur_time)*prices_int[first(ext_edges_DN[dn]),cur_time] for dn=keys(DNs))
	end
	firstdf.tm_bsp_set[firstdf.actor .== "BSPs Dist"] .= 0.0
	firstdf.tm_bsp_set[firstdf.actor .== "BRPs Dist"] .= 0.0

	firstdf.tm_brp_set[firstdf.actor .== "BSPs Tr"] .= 0.0
	firstdf.tm_brp_set[firstdf.actor .== "BRPs Tr"] .= mapreduce(sum,+,init=0,net_injP[i,cur_time]*dualsol.priceP[i,cur_time] for i=tr_nodes)
	firstdf.tm_brp_set[firstdf.actor .== "TSO"] .= -mapreduce(sum,+,init=0,net_injP[i,cur_time]*dualsol.priceP[i,cur_time] for i=nodes)
	firstdf.tm_brp_set[firstdf.actor .== "ADS"] .= 0.0
	firstdf.tm_brp_set[firstdf.actor .== "BSPs Dist"] .= 0.0
	firstdf.tm_brp_set[firstdf.actor .== "BRPs Dist"] .= mapreduce(sum,+,init=0,net_injP[i,cur_time]*dualsol.priceP[i,cur_time] for i=dist_nodes)

	firstdf.ads_disaggregate[firstdf.actor .== "BSPs Tr"] .= 0.0
	firstdf.ads_disaggregate[firstdf.actor .== "BRPs Tr"] .= 0.0
	firstdf.ads_disaggregate[firstdf.actor .== "TSO"] .= 0.0
	firstdf.ads_disaggregate[firstdf.actor .== "ADS"] .= -mapreduce(sum,+,init=0,sol.pg[sb]*dualsol.priceP[sb[1],sb[5]] for sb=dist_seg_bid if sb[5]==cur_time)
	firstdf.ads_disaggregate[firstdf.actor .== "BSPs Dist"] .= mapreduce(sum,+,init=0,sol.pg[sb]*dualsol.priceP[sb[1],sb[5]] for sb=dist_seg_bid if sb[5]==cur_time)
	firstdf.ads_disaggregate[firstdf.actor .== "BRPs Dist"] .= 0.0

	firstdf.tm_bsp_set[firstdf.actor .== "Total"] .= sum(firstdf.tm_bsp_set[firstdf.actor .== act] for act=setdiff(actors,["Total"]))
	firstdf.tm_brp_set[firstdf.actor .== "Total"] .= sum(firstdf.tm_brp_set[firstdf.actor .== act] for act=setdiff(actors,["Total"]))
	firstdf.ads_disaggregate[firstdf.actor .== "Total"] .= sum(firstdf.ads_disaggregate[firstdf.actor .== act] for act=setdiff(actors,["Total"]))

	firstdf.total .= firstdf.tm_bsp_set +firstdf.tm_brp_set +firstdf.ads_disaggregate

	return firstdf
end

function build_ads_solution_table(cur_time::Int64,sol::PFsol,dualsol::Dualsol)::DataFrame
	actors=union(["TN"],["DNs"],["Total"])
	finaldf=DataFrame(
	actor=actors,
	real_bsp=zeros(length(actors)),
	real_brp=zeros(length(actors)),
	reactive_bsp=zeros(length(actors)),
	reactive_brp=zeros(length(actors)),
	min_real_DLMP=zeros(length(actors)),
	max_real_DLMP=zeros(length(actors)),
	av_real_DLMP=zeros(length(actors)),
	min_reactive_DLMP=zeros(length(actors)),
	max_reactive_DLMP=zeros(length(actors)),
	av_reactive_DLMP=zeros(length(actors))
	)

	finaldf.real_bsp[finaldf.actor .== "TN"] .= mapreduce(sum,+,init=0,sol.pg[sb] for sb=tr_seg_bid if sb[5]==cur_time)
	finaldf.real_brp[finaldf.actor .== "TN"] .= mapreduce(sum,+,init=0,net_injP[i,cur_time] for i=tr_nodes)
	finaldf.real_bsp[finaldf.actor .== "DNs"] .= mapreduce(sum,+,init=0,sol.pg[sb] for sb=dist_seg_bid if sb[5]==cur_time)
	finaldf.real_brp[finaldf.actor .== "DNs"] .= mapreduce(sum,+,init=0,net_injP[i,cur_time] for i=dist_nodes)
	finaldf.reactive_bsp[finaldf.actor .== "TN"] .= 0.0
	finaldf.reactive_brp[finaldf.actor .== "TN"] .= 0.0
	finaldf.reactive_bsp[finaldf.actor .== "DNs"] .= mapreduce(sum,+,init=0,sol.qg[sb] for sb=dist_seg_bid if sb[5]==cur_time)
	finaldf.reactive_brp[finaldf.actor .== "DNs"] .= mapreduce(sum,+,init=0,net_injQ[i,cur_time] for i=dist_nodes)

	finaldf.min_real_DLMP[finaldf.actor .== "TN"] .= minimum(dualsol.priceP[i,cur_time] for i=tr_nodes)
	finaldf.max_real_DLMP[finaldf.actor .== "TN"] .= maximum(dualsol.priceP[i,cur_time] for i=tr_nodes)
	finaldf.av_real_DLMP[finaldf.actor .== "TN"] .= mapreduce(sum,+,init=0,dualsol.priceP[i,cur_time] for i=tr_nodes)/length(tr_nodes)
	finaldf.min_real_DLMP[finaldf.actor .== "DNs"] .= minimum(dualsol.priceP[i,cur_time] for i=dist_nodes)
	finaldf.max_real_DLMP[finaldf.actor .== "DNs"] .= maximum(dualsol.priceP[i,cur_time] for i=dist_nodes)
	finaldf.av_real_DLMP[finaldf.actor .== "DNs"] .= mapreduce(sum,+,init=0,dualsol.priceP[i,cur_time] for i=dist_nodes)/length(dist_nodes)
	finaldf.min_reactive_DLMP[finaldf.actor .== "DNs"] .= minimum(dualsol.priceQ[i,cur_time] for i=dist_nodes)
	finaldf.max_reactive_DLMP[finaldf.actor .== "DNs"] .= maximum(dualsol.priceQ[i,cur_time] for i=dist_nodes)
	finaldf.av_reactive_DLMP[finaldf.actor .== "DNs"] .= mapreduce(sum,+,init=0,dualsol.priceQ[i,cur_time] for i=dist_nodes)/length(dist_nodes)

	finaldf.real_bsp[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,sol.pg[sb] for sb=seg_bid if sb[5]==cur_time)
	finaldf.real_brp[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,net_injP[i,cur_time] for i=nodes)
	finaldf.reactive_bsp[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,sol.qg[sb] for sb=dist_seg_bid if sb[5]==cur_time)
	finaldf.reactive_brp[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,net_injQ[i,cur_time] for i=dist_nodes)
	finaldf.min_real_DLMP[finaldf.actor .== "Total"] .= minimum(dualsol.priceP[i,cur_time] for i=nodes)
	finaldf.max_real_DLMP[finaldf.actor .== "Total"] .= maximum(dualsol.priceP[i,cur_time] for i=nodes)
	finaldf.av_real_DLMP[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,dualsol.priceP[i,cur_time] for i=nodes)/length(nodes)
	finaldf.min_reactive_DLMP[finaldf.actor .== "Total"] .= minimum(dualsol.priceQ[i,cur_time] for i=dist_nodes)
	finaldf.max_reactive_DLMP[finaldf.actor .== "Total"] .= maximum(dualsol.priceQ[i,cur_time] for i=dist_nodes)
	finaldf.av_reactive_DLMP[finaldf.actor .== "Total"] .= mapreduce(sum,+,init=0,dualsol.priceQ[i,cur_time] for i=dist_nodes)/length(dist_nodes)
	return finaldf
end

function set_ads_tables_to_CSV(dir::String,sol::PFsol,dualsol::Dualsol,prices_int=0)
	first_table=Dict{Int64,DataFrame}(t=>DataFrame() for t=times)
	solution_table=Dict{Int64,DataFrame}(t=>DataFrame() for t=times)
	for t=times
		first_table[t] = build_ads_set_tables(t,sol,dualsol,prices_int)
	end
	for t=times
		solution_table[t] = build_ads_solution_table(t,sol,dualsol)
	end
	if !(dir in readdir()) mkdir(dir) end
	for t=times
		CSV.write("$dir/set_table1_T$t.csv",first_table[t])
		CSV.write("$dir/solution_table_T$t.csv",solution_table[t])
	end
	total_first_table = DataFrame(
	actor=first_table[first(times)].actor,
	tm_bsp_set=sum(first_table[t].tm_bsp_set for t=times),
	tm_brp_set=sum(first_table[t].tm_brp_set for t=times),
	ads_disaggregate=sum(first_table[t].ads_disaggregate for t=times),
	total=sum(first_table[t].total for t=times)
	)

	total_solution_table = DataFrame(
	actor=solution_table[first(times)].actor,
	real_bsp=sum(solution_table[t].real_bsp for t=times),
	real_brp=sum(solution_table[t].real_brp for t=times),
	reactive_bsp=sum(solution_table[t].reactive_bsp for t=times),
	reactive_brp=sum(solution_table[t].reactive_brp for t=times),
	min_real_DLMP=[minimum(solution_table[t].min_real_DLMP[i] for t=times) for i=1:length(solution_table[first(times)].actor)],
	max_real_DLMP=[maximum(solution_table[t].max_real_DLMP[i] for t=times) for i=1:length(solution_table[first(times)].actor)],
	av_real_DLMP=sum(solution_table[t].av_real_DLMP for t=times)./length(times),
	min_reactive_DLMP=[minimum(solution_table[t].min_reactive_DLMP[i] for t=times) for i=1:length(solution_table[first(times)].actor)],
	max_reactive_DLMP=[maximum(solution_table[t].max_reactive_DLMP[i] for t=times) for i=1:length(solution_table[first(times)].actor)],
	av_reactive_DLMP=sum(solution_table[t].av_reactive_DLMP for t=times)./length(times)
	)
	CSV.write("$dir/total_set_table1.csv",total_first_table)
	CSV.write("$dir/total_solution_table.csv",total_solution_table)
end
