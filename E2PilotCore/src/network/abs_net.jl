
abstract type AbsNet end
Graphs.edges(net::AbsNet) = Graphs.edges(net.g)
Graphs.inneighbors(net::AbsNet, nd::Int64) = inneighbors(net.g, nd)
Graphs.outneighbors(net::AbsNet, nd::Int64) = outneighbors(net.g, nd)
Graphs.indegree(net::AbsNet, nd::Int64) = indegree(net.g, nd)
Graphs.outdegree(net::AbsNet, nd::Int64) = outdegree(net.g, nd)
Graphs.all_neighbors(net::AbsNet, nd::Int64) = all_neighbors(net.g, nd)