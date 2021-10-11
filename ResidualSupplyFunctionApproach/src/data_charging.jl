TrNodesF = string(TestCase,"TransmissionNodes.csv")
DistNodesF = string(TestCase,"DistributionNodes.csv")
EdgesF = string(TestCase,"Edges.csv")
NetInjF = string(TestCase,"NetInjections.csv")
BidsF = string(TestCase,"Bids.csv")
GeneralF = string(TestCase,"GeneralParam.csv")
ExclusiveQtBidsF = string(TestCase,"ExclusiveQtBids.csv")
RampConstrF = string(TestCase,"RampConstraints.csv")
MinDurPairsF = string(TestCase,"MinDurationPairs.csv")
HalfPlanesF = string(TestCase,"HalfPlanes.csv")
QPDiscF = string(TestCase,"QPDisc.csv")

TrNodesdf = DataFrame(CSV.File(TrNodesF))
DistNodesdf = DataFrame(CSV.File(DistNodesF))
Edgesdf = DataFrame(CSV.File(EdgesF))
NetInjdf = DataFrame(CSV.File(NetInjF))
Bidsdf = DataFrame(CSV.File(BidsF))
Generaldf = DataFrame(CSV.File(GeneralF))
ExQtBidsdf = DataFrame(CSV.File(ExclusiveQtBidsF))
RampConstrdf = DataFrame(CSV.File(RampConstrF))
MinDurPairsdf = DataFrame(CSV.File(MinDurPairsF))
HalfPlanesdf = DataFrame(CSV.File(HalfPlanesF))
QPDiscdf = DataFrame(CSV.File(QPDiscF))
println("DONE")

### Charging data into relevant structures (Maybe not necessary) ###
print("Charging data structures...")
### General Data structures ###
times = [i for i=Generaldf[!,:StartTime][1] : Generaldf[!,:EndTime][1]]
#times = [1]
tr_nodes = TrNodesdf[!,:ID]
ref_tr_node = TrNodesdf[!,:ID][[ii for ii=1:nrow(TrNodesdf) if TrNodesdf[!,:RefNode][ii] == 1]]
#angle_fixed=Dict{Int64,Float64}(DistNodesdf[!,:ID][i] => (DistNodesdf[!,:Root][i]-360.0)*pi/180.0
#for i=1:nrow(DistNodesdf) if DistNodesdf[!,:Root][i] != 0.0)
dist_nodes = DistNodesdf[!,:ID]
edges = Edgesdf[!,:ID]
nodes = union(tr_nodes, dist_nodes)
p_max = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(
(Bidsdf[!,:Node][i],Bidsdf[!,:QtBids][i],Bidsdf[!,:QBid][i],Bidsdf[!,:QBidSeg][i], Bidsdf[!,:forT][i])
=> Bidsdf[!,:QuantityHiQ][i] for i=1:nrow(Bidsdf))
p_min = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(
(Bidsdf[!,:Node][i],Bidsdf[!,:QtBids][i],Bidsdf[!,:QBid][i],Bidsdf[!,:QBidSeg][i], Bidsdf[!,:forT][i])
=> Bidsdf[!,:QuantityLoQ][i] for i=1:nrow(Bidsdf))
QuantityDiff=Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(
(Bidsdf[!,:Node][i],Bidsdf[!,:QtBids][i],Bidsdf[!,:QBid][i],Bidsdf[!,:QBidSeg][i], Bidsdf[!,:forT][i])
=> Bidsdf[!,:QuantityHiQ][i]-Bidsdf[!,:QuantityLoQ][i] for i=1:nrow(Bidsdf))
X_lb=Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}((Bidsdf[!,:Node][i],Bidsdf[!,:QtBids][i],Bidsdf[!,:QBid][i],
Bidsdf[!,:QBidSeg][i], Bidsdf[!,:forT][i]) => Bidsdf[!,:Xlb][i] for i=1:nrow(Bidsdf))
X_ub=Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}((Bidsdf[!,:Node][i],Bidsdf[!,:QtBids][i],Bidsdf[!,:QBid][i],
Bidsdf[!,:QBidSeg][i], Bidsdf[!,:forT][i]) => Bidsdf[!,:Xub][i] for i=1:nrow(Bidsdf))
#q_min=Dict{Int64,Float64}(DistNodesdf[!,:ID][i] => DistNodesdf[!,:Qmin][i] for i=1:nrow(DistNodesdf))
#q_max=Dict{Int64,Float64}(DistNodesdf[!,:ID][i] => DistNodesdf[!,:Qmax][i] for i=1:nrow(DistNodesdf))
PriceHiQ = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(
(Bidsdf[!,:Node][i],Bidsdf[!,:QtBids][i],Bidsdf[!,:QBid][i],Bidsdf[!,:QBidSeg][i], Bidsdf[!,:forT][i])
 => Bidsdf[!,:PriceHiQ][i] for i=1:nrow(Bidsdf))
PriceLoQ = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(
(Bidsdf[!,:Node][i],Bidsdf[!,:QtBids][i],Bidsdf[!,:QBid][i],Bidsdf[!,:QBidSeg][i], Bidsdf[!,:forT][i])
=> Bidsdf[!,:PriceLoQ][i] for i=1:nrow(Bidsdf))
MinDurIndex=Set{Tuple{Int64,Int64,Int64,Int64,Int64}}((MinDurPairsdf[!,:QtBids][i],
MinDurPairsdf[!,:QBid1][i],MinDurPairsdf[!,:QBid2][i],MinDurPairsdf[!,:Tau][i],
MinDurPairsdf[!,:ForTime][i]) for i=1:nrow(MinDurPairsdf)
)
# solution = Dict{Tuple{Int64,Int64,Int64,Int64,Int64},Float64}(
# (Bidsdf[!,:Node][i],Bidsdf[!,:QtBids][i],Bidsdf[!,:QBid][i],Bidsdf[!,:QBidSeg][i], Bidsdf[!,:forT][i])
# => Bidsdf[!,:solution][i]*(Bidsdf[!,:QuantityHiQ][i]-Bidsdf[!,:QuantityLoQ][i]) for i=1:nrow(Bidsdf))

# for BS=keys(p_max)
# 	if X_lb[BS] == 1 p_min[BS]=p_max[BS] end
# 	if p_max[BS] < p_min[BS]
# 		tmp=p_max[BS]
# 		p_max[BS]=p_min[BS]
# 		p_min[BS]=tmp
# 		PriceLoQ[BS]=-PriceLoQ[BS]
# 		PriceHiQ[BS]=-PriceHiQ[BS]
# 	end
# end
for BS=keys(p_max)
	if p_max[BS] < 0
		PriceLoQ[BS]=-PriceLoQ[BS]
		PriceHiQ[BS]=-PriceHiQ[BS]
	end
	p_max[BS]-=p_min[BS]
	if p_max[BS] < 0
		p_min[BS]=p_max[BS]
		p_max[BS]=0
	else
		p_min[BS]=0
	end
	if X_lb[BS] == 1
		p_min[BS] += p_max[BS]
		p_max[BS] = p_min[BS]
	end
end

cost = Dict{Tuple{Int64,Int64,Int64,Int64,Int64}, Tuple{Float64,Float64}}(
BS => (PriceLoQ[BS]-p_min[BS]*(PriceHiQ[BS]-PriceLoQ[BS])/(p_max[BS]-p_min[BS]),
(PriceHiQ[BS]-PriceLoQ[BS])/(2*(p_max[BS]-p_min[BS]))) for BS = keys(p_max))
for j=keys(cost)
	if (isnan(cost[j][1]) || isnan(cost[j][2]))
		cost[j]=(0.0,0.0)
	end
end
net_injP=Dict{Tuple{Int64,Int64}, Float64}((NetInjdf[!,:Node][k], NetInjdf[!,:forT][k])
 => NetInjdf[!,:P][k] for k=1:nrow(NetInjdf))
net_injQ=Dict{Tuple{Int64,Int64}, Float64}((NetInjdf[!,:Node][k], NetInjdf[!,:forT][k])
 => NetInjdf[!,:Q][k] for k=1:nrow(NetInjdf))
#valueP=Dict{Tuple{Int64,Int64}, Float64}((NetInjdf[!,:Node][k], NetInjdf[!,:forT][k])
# => NetInjdf[!,:ValueP][k] for k=1:nrow(NetInjdf))
#valueQ=Dict{Tuple{Int64,Int64}, Float64}((NetInjdf[!,:Node][k], NetInjdf[!,:forT][k])
# => NetInjdf[!,:ValueQ][k] for k=1:nrow(NetInjdf))
e_fr = Dict{Int64,Int64}(Edgesdf[!,:ID][i] => Edgesdf[!,:NodeFrom][i] for i=1:nrow(Edgesdf))
e_to = Dict{Int64,Int64}(Edgesdf[!,:ID][i] => Edgesdf[!,:NodeTo][i] for i=1:nrow(Edgesdf))
vsq_min=Dict{Int64,Float64}(DistNodesdf[!,:ID][i] => DistNodesdf[!,:Vsqmin][i] for i=1:nrow(DistNodesdf))
vsq_max=Dict{Int64,Float64}(DistNodesdf[!,:ID][i] => DistNodesdf[!,:Vsqmax][i] for i=1:nrow(DistNodesdf))
rate_a=Dict{Int64,Float64}(Edgesdf[!,:ID][i] => Edgesdf[!,:S][i] for i=1:nrow(Edgesdf))
tap=Dict{Int64,Float64}(Edgesdf[!,:ID][i] => Edgesdf[!,:Tap][i] for i=1:nrow(Edgesdf))
shift=Dict{Int64,Float64}(Edgesdf[!,:ID][i] => Edgesdf[!,:Shift][i] for i=1:nrow(Edgesdf))
tap_r=Dict{Int64,Float64}(i=> tap[i]*cos(shift[i]) for i=edges)
tap_i=Dict{Int64,Float64}(i=> tap[i]*sin(shift[i]) for i=edges)
seg_bid=Set{Tuple{Int64,Int64,Int64,Int64,Int64}}()
bid=Set{Tuple{Int64,Int64,Int64,Int64}}()
qt_bid=Set{Int64}(Bidsdf[!,:QtBids])
tr_edges = Set{Int64}()
dist_edges = Set{Int64}()
root_dist_node = Set{Int64}()
for e=edges
	if e_fr[e] in dist_nodes && e_to[e] in dist_nodes
		push!(dist_edges,e)
	elseif e_fr[e] in dist_nodes && e_to[e] in tr_nodes
		push!(tr_edges, e)
		push!(root_dist_node,e_fr[e])
	elseif e_fr[e] in tr_nodes && e_to[e] in dist_nodes
		push!(tr_edges, e)
		push!(root_dist_node,e_to[e])
	else
		push!(tr_edges,e)
	end
end
slack_bus=0
if isempty(ref_tr_node) slack_bus=root_dist_node[1] else slack_bus=ref_tr_node[1] end
for i=union(tr_nodes,root_dist_node) for t=times net_injQ[i,t]=0 end end

for i=1:nrow(Bidsdf)
	push!(seg_bid,(Bidsdf[!,:Node][i],Bidsdf[!,:QtBids][i], Bidsdf[!,:QBid][i], Bidsdf[!,:QBidSeg][i], Bidsdf[!,:forT][i]))
	push!(bid,(Bidsdf[!,:Node][i],Bidsdf[!,:QtBids][i], Bidsdf[!,:QBid][i], Bidsdf[!,:forT][i]))
end
AlphaOmegaIndex=Set{Tuple{Int64,Int64,Int64,Int64}}()
NoNewActIndex=Set{Tuple{Int64,Int64,Int64,Int64}}()
if hasproperty(Bidsdf,:AlphaOmegaSet)
	for i=1:nrow(Bidsdf)
		if Bidsdf[!,:AlphaOmegaSet][i]==1
			push!(AlphaOmegaIndex,(Bidsdf[!,:Node][i],Bidsdf[!,:QtBids][i], Bidsdf[!,:QBid][i], Bidsdf[!,:forT][i]))
		end
	end
else
	AOSF = string(TestCase,"AlphaOmegaSet.csv")
	AOSdf = DataFrame(CSV.File(AOSF))
	AOSIndex=Set{Tuple{Int64,Int64,Int64}}()
	for i=1:nrow(AOSdf)
		push!(AOSIndex,(AOSdf[!,:QtBids][i],AOSdf[!,:QBid][i],AOSdf[!,:forT][i]))
	end
	for bd=bid
		for aos=AOSIndex
			if aos[1]==bd[2] && aos[2]==bd[3] && aos[3]==bd[4]
				push!(AlphaOmegaIndex,bd)
			end
		end
	end
end
if hasproperty(Bidsdf,:NoNewAct)
	for i=1:nrow(Bidsdf)
		if Bidsdf[!,:NoNewAct][i]==1
			push!(NoNewActIndex,(Bidsdf[!,:Node][i],Bidsdf[!,:QtBids][i], Bidsdf[!,:QBid][i], Bidsdf[!,:forT][i]))
		end
	end
else
	NNAF = string(TestCase,"NoNewAct.csv")
	NNAdf = DataFrame(CSV.File(NNAF))
	NNAIndex=Set{Tuple{Int64,Int64,Int64}}()
	for i=1:nrow(NNAdf)
		push!(NNAIndex,(NNAdf[!,:QtBids][i],NNAdf[!,:QBid][i],NNAdf[!,:forT][i]))
	end
	for bd=bid
		for aos=NNAIndex
			if aos[1]==bd[2] && aos[2]==bd[3] && aos[3]==bd[4]
				push!(NoNewActIndex,bd)
			end
		end
	end
end

AlphaOmegaLate=Set{Tuple{Int64,Int64,Int64,Int64}}()
AlphaOmegaEarly=Set{Tuple{Int64,Int64,Int64,Int64}}()
if length(AlphaOmegaIndex)!=0
	AlphaOmegaLate=Set{Tuple{Int64,Int64,Int64,Int64}}(BB for BB=AlphaOmegaIndex for forT=times
	#if [BB[1],BB[2],get(QtForT2Q,[BB[2],BB[4]-1],-1),BB[4]-1] in AlphaOmegaIndex]
	if (BB[1],BB[2],BB[3],forT) in AlphaOmegaIndex && BB[4] < forT)
	AlphaOmegaEarly=Set{Tuple{Int64,Int64,Int64,Int64}}(BB for BB=AlphaOmegaIndex for forT=times
	#if !([BB[1],BB[2],get(QtForT2Q,[BB[2],BB[4]-1],-1),BB[4]-1] in AlphaOmegaIndex)]
	if (BB[1],BB[2],BB[3],forT) in AlphaOmegaIndex && BB[4]==forT)
end

QPDiscID=Set{Tuple{Int64,Int64,Int64}}((QPDiscdf[!,:QtBid][i],QPDiscdf[!,:QBids][i],QPDiscdf[!,:QPDCID][i])
for i=1:nrow(QPDiscdf))
QPDiscFrom=Dict{Tuple{Int64,Int64,Int64},Int64}((QPDiscdf[!,:QtBid][i],QPDiscdf[!,:QBids][i],QPDiscdf[!,:QPDCID][i])
=> QPDiscdf[!,:QPDFrom][i] for i=1:nrow(QPDiscdf))
QPDiscTo=Dict{Tuple{Int64,Int64,Int64},Int64}((QPDiscdf[!,:QtBid][i],QPDiscdf[!,:QBids][i],QPDiscdf[!,:QPDCID][i])
=> QPDiscdf[!,:QPDTo][i] for i=1:nrow(QPDiscdf))
QPDiscMax=Dict{Tuple{Int64,Int64,Int64},Float64}((QPDiscdf[!,:QtBid][i],QPDiscdf[!,:QBids][i],QPDiscdf[!,:QPDCID][i])
=> QPDiscdf[!,:QPDMax][i] for i=1:nrow(QPDiscdf))
QPDiscConstrIndex=Set{Tuple{Int64,Int64,Int64,Int64,Int64}}()
for bd=bid
	for qpd=QPDiscID
		if qpd[1]==bd[2] && qpd[2]==bd[3] && bd[1] in dist_nodes &&
			QPDiscTo[qpd]==-1 && QPDiscFrom[qpd]==-1
			push!(QPDiscConstrIndex, (bd[1],bd[2],bd[3],bd[4],qpd[3]))
		end
	end
end
QPHPID=Set{Tuple{Int64,Int64}}((HalfPlanesdf[!,:ID][i],HalfPlanesdf[!,:QtBid][i])
for i=1:nrow(HalfPlanesdf))
HP2Qt=Dict{Int64,Int64}(HalfPlanesdf[!,:ID][i] => HalfPlanesdf[!,:QtBid][i]
for i=1:nrow(HalfPlanesdf))
QPHPFrom=Dict{Tuple{Int64,Int64},Int64}((HalfPlanesdf[!,:ID][i],HalfPlanesdf[!,:QtBid][i])
=> HalfPlanesdf[!,:QPHPFrom][i] for i=1:nrow(HalfPlanesdf))
QPHPTo=Dict{Tuple{Int64,Int64},Int64}((HalfPlanesdf[!,:ID][i],HalfPlanesdf[!,:QtBid][i])
=> HalfPlanesdf[!,:QPHPTo][i] for i=1:nrow(HalfPlanesdf))
QPHPSlope=Dict{Tuple{Int64,Int64},Float64}((HalfPlanesdf[!,:ID][i],HalfPlanesdf[!,:QtBid][i])
=> HalfPlanesdf[!,:QPHPSlope][i] for i=1:nrow(HalfPlanesdf))
QPHPOffset=Dict{Tuple{Int64,Int64},Int64}((HalfPlanesdf[!,:ID][i],HalfPlanesdf[!,:QtBid][i])
=> HalfPlanesdf[!,:QPHPOffset][i] for i=1:nrow(HalfPlanesdf))
QPHPConstr=Dict{Tuple{Int64,Int64},Int64}((HalfPlanesdf[!,:ID][i],HalfPlanesdf[!,:QtBid][i])
=> HalfPlanesdf[!,:QPHPConstr][i] for i=1:nrow(HalfPlanesdf))
QPHPConstrIndex1=Set{Tuple{Int64,Int64,Int64,Int64,Int64}}()
QPHPConstrIndex0=Set{Tuple{Int64,Int64,Int64,Int64,Int64}}()
for bd=bid
	for hpc=QPHPID
		if hpc[2]==bd[2] && bd[1] in dist_nodes #&& QPHPTo[hpc]==-1 && QPHPFrom[hpc]==-1
			if QPHPConstr[hpc]==1
				push!(QPHPConstrIndex1, (bd[1],bd[2],bd[3],bd[4],hpc[1]))
			end
			if QPHPConstr[hpc]==0
				push!(QPHPConstrIndex0, (bd[1],bd[2],bd[3],bd[4],hpc[1]))
			end
		end
	end
end
# orig_seg_bid=copy(seg_bid)
# for sb=orig_seg_bid
# 	#if (p_min[sb] < 0 && cost[sb][1] < 0) || cost[sb][1]==0
# 	# if cost[sb][1]==0
# 	# 	delete!(seg_bid,sb)
# 	# 	net_injP[sb[1],sb[5]]+=p_max[sb]-p_min[sb]
# 	# end
# 	if p_min[sb]== 0.0 && p_max[sb] == 0.0
# 		delete!(seg_bid,sb)
# 	end
# end
# orig_seg_bid=copy(seg_bid)
# for sb1=orig_seg_bid
# 	for sb2=orig_seg_bid
# 		if sb1 == sb2 continue end
# 		if sb1 in seg_bid && sb2 in seg_bid
# 			if (sb1[1]==sb2[1] && #sb1[2]==sb2[2] && sb1[3]==sb2[3] &&
# 				sb1[5]==sb2[5] && cost[sb1] == cost[sb2])
# 				#println(sb1, sb2);  break
# 				delete!(seg_bid,sb1)
# 				p_min[sb2]=min(0,p_min[sb2])+min(0,p_min[sb1])
# 				p_max[sb2]=max(0,p_min[sb2])+max(0,p_min[sb1])
# 			end
# 		end
# 	end
# end

# orig_seg_bid=copy(seg_bid)
# for sb=orig_seg_bid
# 	if p_min[sb]== 0.0 && p_max[sb] == 0.0
# 		delete!(seg_bid,sb)
# 	end
# end

bid_2_seg_bids=Dict{Tuple{Int64,Int64,Int64,Int64},Set{Tuple{Int64,Int64,Int64,Int64,Int64}}}()
for BS=seg_bid
	get!(bid_2_seg_bids,(BS[1],BS[2],BS[3],BS[5]),
	push!(get(bid_2_seg_bids,(BS[1],BS[2],BS[3],BS[5]),Set{Tuple{Int64,Int64,Int64,Int64,Int64}}()),BS))
end

orig_bid=copy(bid)
for B=orig_bid
	if isempty(get(bid_2_seg_bids,B,Set{Tuple{Int64,Int64}}()))
		delete!(bid,B)
	end
end

orig_aoe=copy(AlphaOmegaEarly)
for aoe=orig_aoe
	if !(aoe in bid)
		delete!(AlphaOmegaEarly,aoe)
	end
end
orig_aol=copy(AlphaOmegaLate)
for aol=orig_aol
	if !(aol in bid)
		delete!(AlphaOmegaLate,aol)
	end
end
orig_aoi=copy(AlphaOmegaIndex)
for aoi=orig_aoi
	if !(aoi in bid)
		delete!(AlphaOmegaIndex,aoi)
	end
end
orig_nni=copy(NoNewActIndex)
for nni=orig_nni
	if !(nni in bid)
		delete!(NoNewActIndex,nni)
	end
end
# orig_net_injP=copy(net_injP)
# orig_net_injQ=copy(net_injQ)
# for i=keys(orig_net_injP)
# 	if orig_net_injP[i] > 0
# 		delete!(net_injP,i)
# 		#delete!(net_injQ,i)
# 	end
# end

tr_seg_bid=[sb for sb=seg_bid if sb[1] in tr_nodes]
dist_seg_bid=[sb for sb=seg_bid if sb[1] in dist_nodes]
dist_bid=[bd for bd=bid if bd[1] in dist_nodes]
tr_bid=[bd for bd=bid if bd[1] in tr_nodes]

n_fr=Dict{Int64,Array{Int64}}()
n_to=Dict{Int64,Array{Int64}}()
dist_n_fr=Dict{Int64,Array{Int64}}()
dist_n_to=Dict{Int64,Array{Int64}}()
tr_n_fr=Dict{Int64,Array{Int64}}()
tr_n_to=Dict{Int64,Array{Int64}}()
nt_2_bids=Dict{Tuple{Int64,Int64},Set{Tuple{Int64,Int64,Int64,Int64}}}()
for i=nodes
	get!(n_fr,i,[e for e=edges if e_fr[e]==i])
	get!(n_to,i,[e for e=edges if e_to[e]==i])
	get!(dist_n_fr,i,[e for e=dist_edges if e_fr[e]==i])
	get!(dist_n_to,i,[e for e=dist_edges if e_to[e]==i])
	get!(tr_n_fr,i,[e for e=tr_edges if e_fr[e]==i])
	get!(tr_n_to,i,[e for e=tr_edges if e_to[e]==i])
end
for B=bid
	get!(nt_2_bids,(B[1],B[4]),
		push!(get(nt_2_bids,(B[1],B[4]),Set{Tuple{Int64,Int64,Int64,Int64}}()),B))
end
pv_nodes=Set{Int64}()
for i=nodes
	for t=times
		if (i,t) in keys(nt_2_bids)
			push!(pv_nodes, i)
			continue
		end
	end
end
pq_nodes=setdiff(nodes,pv_nodes)
setdiff!(pv_nodes,slack_bus)
setdiff!(pq_nodes,slack_bus)
line_limit=Set{Int64}([e for e=edges if rate_a[e] > 0 && rate_a[e] < 90000])
v_limit=Set{Int64}([i for i=dist_nodes if vsq_min[i]!=0 && vsq_max[i]!=0])

s_index=Set{Tuple{Int64,Int64}}([(e_fr[e],e_to[e]) for e=dist_edges])
c_index=copy(s_index)
union!(c_index,Set{Tuple{Int64,Int64}}((i,i) for i=dist_nodes))
tr_nodes_ext=union!([e_fr[e] for e=tr_edges],[e_to[e] for e=tr_edges])
R_pu=Dict{Int64,Float64}(Edgesdf[!,:ID][i] => Edgesdf[!,:R][i] for i=1:nrow(Edgesdf))
for e=tr_edges R_pu[e]=0 end
X_pu=Dict{Int64,Float64}(Edgesdf[!,:ID][i] => Edgesdf[!,:X][i] for i=1:nrow(Edgesdf))
g=Dict{Int64,Float64}(
    i=> R_pu[i]/(R_pu[i]^2+X_pu[i]^2) for i=edges
)
for j=keys(g) if isnan(g[j]) g[j]=0.0 end end
b=Dict{Int64,Float64}(
    i=> -X_pu[i]/(R_pu[i]^2+X_pu[i]^2) for i=edges
)
for j=keys(b) if isnan(b[j]) b[j]=0.0 end end
g_to=Dict{Int64,Float64}(i=> 0.0 for i=edges)
g_fr=Dict{Int64,Float64}(i=> 0.0 for i=edges)
gs=Dict{Int64,Float64}(DistNodesdf[!,:ID][i] => DistNodesdf[!,:G][i] for i=1:nrow(DistNodesdf))
merge!(gs,Dict{Int64,Float64}(i => 0.0 for i=tr_nodes))
b_to=Dict{Int64,Float64}(Edgesdf[!,:ID][i] => Edgesdf[!,:B][i]/2 for i=1:nrow(Edgesdf))
b_fr=Dict{Int64,Float64}(Edgesdf[!,:ID][i] => Edgesdf[!,:B][i]/2 for i=1:nrow(Edgesdf))
bs=Dict{Int64,Float64}(DistNodesdf[!,:ID][i] => DistNodesdf[!,:B][i] for i=1:nrow(DistNodesdf))
Gbus=Dict{Tuple{Int64,Int64},Float64}((i,j) => 0.0 for i=nodes for j=nodes)
Bbus=Dict{Tuple{Int64,Int64},Float64}((i,j) => 0.0 for i=nodes for j=nodes)
for e=edges
	Gbus[(e_fr[e],e_to[e])] = -g[e]
	Gbus[(e_to[e],e_fr[e])] = Gbus[(e_fr[e],e_to[e])]
	Bbus[(e_fr[e],e_to[e])] = -b[e]
	Bbus[(e_to[e],e_fr[e])] = Bbus[(e_fr[e],e_to[e])]
end
for i=dist_nodes
	Gbus[(i,i)] = -sum(Gbus[(i,j)] for j=dist_nodes)
	Bbus[(i,i)] = -sum(Bbus[(i,j)] for j=dist_nodes)
end

function dist_neighbors(n::Int64)::Set{Int64}
	return Set{Int64}(union([e_fr[i] for i=dist_n_to[n]],[e_to[i] for i=dist_n_fr[n]]))
end

#ex_qt_bid_id=Set{Int64}(ExQtBidsdf[!,:ID])
ex_qt_bid=Dict{Int64,Set{Int64}}()
for i=1:nrow(ExQtBidsdf)
	if !(ExQtBidsdf[!,:QtBid][i] in qt_bid) continue end
	get!(ex_qt_bid,ExQtBidsdf[!,:ID][i],
	push!(get(ex_qt_bid,ExQtBidsdf[!,:ID][i],Set{Int64}()),ExQtBidsdf[!,:QtBid][i]))
end
ex_qt_bid_id=Set{Int64}(keys(ex_qt_bid))

RCRealPowerIncr = Dict{Tuple{Int64,Int64,Int64}, Float64}(
(RampConstrdf[!,:QtBids][i],RampConstrdf[!,:QBid][i],RampConstrdf[!,:RampConstr][i])
=> RampConstrdf[!,:RCRealPowerIncr][i] for i=1:nrow(RampConstrdf))
RCMin = Dict{Tuple{Int64,Int64,Int64}, Int64}(
(RampConstrdf[!,:QtBids][i],RampConstrdf[!,:QBid][i],RampConstrdf[!,:RampConstr][i])
=> RampConstrdf[!,:RCMin][i] for i=1:nrow(RampConstrdf))
RampIndex = Set{Tuple{Int64,Int64,Int64}}()
RampIndexMinIncr = Set{Tuple{Int64,Int64,Int64,Int64}}()
RampIndexMaxIncr = Set{Tuple{Int64,Int64,Int64,Int64}}()
for i=1:nrow(RampConstrdf)
	push!(RampIndex,(RampConstrdf[!,:QtBids][i],RampConstrdf[!,:QBid][i],RampConstrdf[!,:RampConstr][i]))
end
for B=bid
	for rc=RampIndex
		if get(RCMin,(B[2],B[3],rc[3]),-1) == 1 &&
		B[2]==rc[1] && B[3] == rc[2] && B[4] < maximum(times)
			push!(RampIndexMinIncr,(B[2],B[3],B[4],rc[3]))
		end
		if get(RCMin,(B[2],B[3],rc[3]),-1) == 0 &&
		B[2]==rc[1] && B[3] == rc[2] && B[4] < maximum(times)
			push!(RampIndexMaxIncr,(B[2],B[3],B[4],rc[3]))
		end
	end
end

DNs=Dict{String,Set{Int64}}()
edges_DN=Dict{String,Set{Int64}}()
for i=root_dist_node
	setNodes=union(Set{Int64}(i),dist_neighbors(i))
	newNodes=copy(setNodes)
	sizeN=0
	while sizeN!=length(newNodes)
		sizeN=length(newNodes)
		for j=newNodes
			union!(newNodes,dist_neighbors(j))
		end
	end
	get!(DNs,"DN-$i",newNodes)
	newEdges=Set{Int64}()
	for e=dist_edges
		if e_fr[e] in newNodes && e_to[e] in newNodes
			push!(newEdges,e)
		end
	end
	get!(edges_DN,"DN-$i",newEdges)
end
ext_edges_DN=Dict{String,Set{Int64}}(
	dn => Set{Int64}([e for e=tr_edges if (e_fr[e] in DNs[dn] && !(e_to[e] in DNs[dn]) || (e_to[e] in DNs[dn] && !(e_fr[e] in DNs[dn])))])
		for dn=keys(DNs)
)
ext_nodes_DN=Dict{String,Set{Int64}}(
	dn => Set{Int64}(union([e_fr[e] for e=ext_edges_DN[dn] if e_to[e] in DNs[dn]],
	[e_to[e] for e=ext_edges_DN[dn] if e_fr[e] in DNs[dn]]))
	for dn=keys(DNs)
)
s_index_DN=Dict{String,Set{Tuple{Int64,Int64}}}(dn => Set{Tuple{Int64,Int64}}() for dn=keys(DNs))
for dn=keys(edges_DN) for e=edges_DN[dn]
		push!(s_index_DN[dn],(e_fr[e],e_to[e]))
end end
c_index_DN=Dict{String,Set{Tuple{Int64,Int64}}}(dn => copy(s_index_DN[dn]) for dn=keys(DNs))
for dn=keys(edges_DN) for i=DNs[dn] push!(c_index_DN[dn],(i,i)) end end

seg_bid_DN=Dict{String,Set{Tuple{Int64,Int64,Int64,Int64,Int64}}}(
	dn => Set{Tuple{Int64,Int64,Int64,Int64,Int64}}([sb for sb=seg_bid if sb[1] in DNs[dn]]) for dn=keys(DNs)
)
dist_seg_bid_DN=Dict{String,Set{Tuple{Int64,Int64,Int64,Int64,Int64}}}(
	dn => Set{Tuple{Int64,Int64,Int64,Int64,Int64}}([sb for sb=dist_seg_bid if sb[1] in DNs[dn]]) for dn=keys(DNs)
)
bid_DN=Dict{String,Set{Tuple{Int64,Int64,Int64,Int64}}}(
	dn => Set{Tuple{Int64,Int64,Int64,Int64}}([(sb[1],sb[2],sb[3],sb[5]) for sb=seg_bid_DN[dn]]) for dn=keys(DNs)
)
dist_bid_DN=Dict{String,Set{Tuple{Int64,Int64,Int64,Int64}}}(
	dn => Set{Tuple{Int64,Int64,Int64,Int64}}([(sb[1],sb[2],sb[3],sb[5]) for sb=dist_seg_bid_DN[dn]]) for dn=keys(DNs)
)
qt_bid_DN=Dict{String,Set{Int64}}(
	dn => Set{Int64}([bd[2] for bd=bid_DN[dn]]) for dn=keys(DNs)
)
ex_qt_bid_id_DN=Dict{String,Set{Int64}}(
	dn => Set{Int64}([i for i=ex_qt_bid_id if first(ex_qt_bid[i]) in qt_bid_DN[dn]]) for dn=keys(DNs)
)
MinDurIndex_DN=Dict{String,Set{Tuple{Int64,Int64,Int64,Int64,Int64}}}(
	dn => Set{Tuple{Int64,Int64,Int64,Int64,Int64}}([mdi for mdi=MinDurIndex if mdi[1] in qt_bid_DN[dn]]) for dn=keys(DNs)
)
tr_qt_bid=Set{Int64}([bd[2] for bd=tr_bid])
tr_ex_qt_bid_id=Set{Int64}([i for i=ex_qt_bid_id if first(ex_qt_bid[i]) in tr_qt_bid])
tr_MinDurIndex=Set{Tuple{Int64,Int64,Int64,Int64,Int64}}([mdi for mdi=MinDurIndex if mdi[1] in tr_qt_bid])


RampIndexMinIncr_DN=Dict{String,Set{Tuple{Int64,Int64,Int64,Int64}}}()
RampIndexMaxIncr_DN=Dict{String,Set{Tuple{Int64,Int64,Int64,Int64}}}()
for dn=keys(DNs)
	newRampMin=Set{Tuple{Int64,Int64,Int64,Int64}}()
	newRampMax=Set{Tuple{Int64,Int64,Int64,Int64}}()
	for bd=bid_DN[dn]
		for rcI=RampIndexMinIncr
			if rcI[1]==bd[2] && rcI[2]==bd[3] && rcI[3]==bd[4]
				push!(newRampMin,rcI)
			end
		end
		for rcI=RampIndexMaxIncr
			if rcI[1]==bd[2] && rcI[2]==bd[3] && rcI[3]==bd[4]
				push!(newRampMax,rcI)
			end
		end
	end
	RampIndexMinIncr_DN[dn]=newRampMin
	RampIndexMaxIncr_DN[dn]=newRampMax
end
tr_RampIndexMinIncr=Set{Tuple{Int64,Int64,Int64,Int64}}()
tr_RampIndexMaxIncr=Set{Tuple{Int64,Int64,Int64,Int64}}()
for bd=tr_bid
	for rcI=RampIndexMinIncr
		if rcI[1]==bd[2] && rcI[2]==bd[3] && rcI[3]==bd[4]
			push!(tr_RampIndexMinIncr,rcI)
		end
	end
	for rcI=RampIndexMaxIncr
		if rcI[1]==bd[2] && rcI[2]==bd[3] && rcI[3]==bd[4]
			push!(tr_RampIndexMaxIncr,rcI)
		end
	end
end
ext_edges=[e for e=tr_edges if ((e_fr[e] in tr_nodes && e_to[e] in dist_nodes) || (e_to[e] in tr_nodes && e_fr[e] in dist_nodes))]
ext_nodes=Set{Int64}(union([e_fr[e] for e=ext_edges],[e_to[e] for e=ext_edges]))
bi_exclusion=Set{Tuple{Tuple{Int64,Int64,Int64,Int64,Int64},Tuple{Int64,Int64,Int64,Int64,Int64}}}()
for sb1=seg_bid
	for sb2=[sb for sb=seg_bid if sb1[1]==sb[1] && sb1[5]==sb[5] && sb[2]==sb[2] && sign(p_max[sb1]-p_min[sb1]) != sign(p_max[sb]-p_min[sb])]
		push!(bi_exclusion,(sb1,sb2))
	end
end
cost_ls=1e4
