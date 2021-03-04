/*
 * graph.h
 *
 *  Created on: 4 march 2021
 *      Author: Joel Nilsson
 *
 *      Graph class template
 *
 *      Dijkstras SSSP algorithm using priority
 *      queue, O((V+E)logV) complexity.
 *
 *      Time table SSSP-problems (Dijkstra modification).
 *
 *      Bellman-Ford algorithm, O(VE) complexity.
 *
 *      Floyd-Warshall algorithm, O(V^3) complexity.
 *
 *      Mximum flow computation, Edmond-Karp O(V*E^2).
 *      Requires that zero-flow is allowed (no non-zero lower limit capacity).
 *
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include<vector>
#include<queue>

template<typename T>
class graph
{
public:

	graph<T>(int nodes);
	void init_EL();
	void init_AL();
	void init_AM();
	void init_time_constraints();
	void init_AM_for_FW();

	void add_edge_to_EL(T weight, int node1, int node2);

	void add_edge_to_AL(T weight, int node1, int node2);
	void add_edge_to_temporal_AL(T weight, int node1, int node2, T start_time, T period);

	void add_edge_to_AM(T weight, int node1, int node2);
	void add_cheapest_edge_to_AM(T weight, int node1, int node2);
	void increment_edge_in_AM(T weight, int node1, int node2);

	std::vector<T> dijkstra(int start);
	std::vector<T> temporal_dijkstra(int start);
	std::vector<int> get_shortest_path(int start, int goal);
	std::vector<T>& bellman_ford(int start);

	std::vector<std::vector<T>> floyd_warshall();
	bool negative_cycle(int i, int j);

	T maxflow(int source, int sink);
	bool hasAugmentingPath(int source, int sink);

	std::vector<std::vector<std::pair<T, int>>> getAL();
	std::vector<std::vector<T>> getAM();
	T getINF();


private:

	int N;

	// Edge List (EL), Adjacency List (AL), Adjacency Matrix (AM)
	std::vector<std::pair<T, std::pair<int, int>>> EL;
	std::vector<std::vector<std::pair<T, int>>> AL;
	std::vector<std::vector<T>> AM;

	std::vector<std::vector<std::pair<T, T>>> time_constraint;
	std::vector<int> predecessor;

	T const INF = 999999999;
};

// Initialize graph with n nodes
template<typename T>
graph<T>::graph(int n)
{
	N = n;
}

/*
 * Initializes empty edge list. The weight is
 * the first element of the pair.
 */
template<typename T>
void graph<T>::init_EL()
{
	EL = std::vector<std::pair<T, std::pair<int, int>>>();
}

/*
 * Initializes adjacency list with N nodes, with no edges.
 * Weight is the first element of subvector.
 */
template<typename T>
void graph<T>::init_AL()
{
	AL = std::vector<std::vector<std::pair<T, int>>>(N, std::vector<std::pair<T, int>>());
}

/*
 * Initializes adjacency matrix to NxN with all zero entries.
 */
template<typename T>
void graph<T>::init_AM()
{
	AM = std::vector<std::vector<T>>(N, std::vector<T>(N, 0));
}

/*
 * Initializes time constraint matrix for time-table Dijkstra.
 */
template<typename T>
void graph<T>::init_time_constraints()
{
	time_constraint = std::vector<std::vector<std::pair<T, T>>>(N, std::vector<std::pair<T, T>>());
}

/*
 * Initializes all off-diagonal elements to
 * "infinity" (INF) for Floyd-Warshall algorithm.
 */
template<typename T>
void graph<T>::init_AM_for_FW()
{
	AM = std::vector<std::vector<T>>(N, std::vector<T>(N, INF));
	for(int i{0}; i < N; ++i)
	{
		AM[i][i] = 0;
	}
}

template<>
graph<long long>::graph(int nodes)
{
	N = nodes;
	INF = 1000000000000;
}

/*
 * Add edge to edge list. Does not check for potential multi-edges.
 */
template<typename T>
void graph<T>::add_edge_to_EL(T weight, int node1, int node2)
{
	EL.push_back( { weight, {node1, node2} } );
}

/*
 * Add edge to adjacency list. Does not check for potential multi-edges.
 */
template<typename T>
void graph<T>::add_edge_to_AL(T weight, int node1, int node2)
{
	AL[node1].push_back( { weight, node2 } );
}


/*
 * Add edge to adjacency matrix. Previous value is overwritten.
 */
template<typename T>
void graph<T>::add_edge_to_AM(T weight, int node1, int node2)
{
	AM[node1][node2] = weight;
}


/*
 * Increase value of position in adjacency matrix. Used for
 * e.g. capacity increase in case of multi-edges in maximum
 * flow problems.
 */
template<typename T>
void graph<T>::increment_edge_in_AM(T weight, int node1, int node2)
{
	AM[node1][node2] += weight;
}

template<typename T>
void graph<T>::add_cheapest_edge_to_AM(T weight, int node1, int node2)
{
	// Adds only cheaper edges. If edges (except for self-edges)
	// are set to INF in e.g. Floyd-Warshall, an edge can of course
	// be added without problem since (weight < INF) holds true.
	// Adding self edges with positive weights will not succeed.
	if(AM[node1][node2] > weight)
	{
		AM[node1][node2] = weight;
	}
}

/*
 * Adds constraint for time-table Dijkstra.
 */
template<typename T>
void graph<T>::add_edge_to_temporal_AL(T weight, int node1, int node2, T start_time, T period)
{
	AL[node1].push_back( { weight, node2 } );
	time_constraint[node1].push_back( { start_time, period } );
}


/*
 * Standard Dijkstra algorithm, O((V+E)logV).
 * Usable for graphs with no negative cycles.
 * Can handle negative weights.
 */
template<typename T>
std::vector<T> graph<T>::dijkstra(int start)
{
	int nodes = AL.size();
	std::vector<T> shortest_distance(nodes, INF);
	predecessor = std::vector<int>(nodes, -1);

	shortest_distance[start] = 0;

	std::priority_queue<std::pair<T, int>, std::vector<std::pair<T, int>>, std::greater<std::pair<T, int>>> q{};

	q.push( { shortest_distance[start], start } );

	while(!q.empty())
	{
		int u = q.top().second;
		T d = q.top().first;
		q.pop();

		if(d > shortest_distance[u])
		{
			continue;
		}

		int NEIGH = AL[u].size();
		for(int i{0}; i < NEIGH; ++i)
		{
			int to = AL[u][i].second;
			T dist = AL[u][i].first;
			// We can never decrease distance by going back, so those edges will not be considered
			if(shortest_distance[u] + dist >= shortest_distance[to])
			{
				continue;
			}
			shortest_distance[to] = shortest_distance[u] + dist;
			predecessor[to] = u;
			q.push( { shortest_distance[to], to } );
		}
	}

	return shortest_distance;
}

/*
 * Dijkstra variant for handling time-table SSSP-problems.
 * If time_constraints contains the start time t_0 for when the edge
 * is allowed to be used, and the interval P with which the edge
 * can be used thereafter (i.e. t = t_0 + k*P, k integer).
 */
template<typename T>
std::vector<T> graph<T>::temporal_dijkstra(int start)
{
	int nodes = AL.size();
	std::vector<T> shortest_distance(nodes, INF);
	predecessor = std::vector<int>(nodes, -1);

	shortest_distance[start] = 0;

	std::priority_queue<std::pair<T, int>, std::vector<std::pair<T, int>>, std::greater<std::pair<T, int>>> q{};

	q.push( { shortest_distance[start], start } );

	while(!q.empty())
	{
		int u = q.top().second;
		T d = q.top().first;
		q.pop();

		if(d > shortest_distance[u]) // ?
		{
			continue;
		}

		int NEIGH = AL[u].size();
		for(int i{0}; i < NEIGH; ++i)
		{
			int to = AL[u][i].second;
			T dist = AL[u][i].first;

			T t_0 = time_constraint[u][i].first;
			T period = time_constraint[u][i].second;

			if(shortest_distance[u] < t_0)
			{
				dist += (t_0 - shortest_distance[u]);
			}
			else if(shortest_distance[u] != t_0 && period == 0) // if time expired (period=0 for this node, and we are after t_0): continue
			{
				continue;
			}
			else if( ( (shortest_distance[u] - t_0) % period ) != 0 && period != 0) // Arriving to u at time shrt_dist[u], need to be at the correct time point to cont.
			{
				dist += ( period - ( (shortest_distance[u] - t_0) % period ) );
			}

			// We can never decrease distance by going back, so those edges will not be considered
			if(shortest_distance[u] + dist >= shortest_distance[to])
			{
				continue;
			}
			shortest_distance[to] = shortest_distance[u] + dist;
			predecessor[to] = u;
			q.push( { shortest_distance[to], to } );
		}
	}

	return shortest_distance;
}

/*
 * Returns shortest path for Dijsktra and Bellman-Ford algorithms
 * Note that the order of the returned path has to be reversed if
 * the reverse statement is not included.
 */
template<typename T>
std::vector<int> graph<T>::get_shortest_path(int start, int goal)
{
	std::vector<int> path{};

	for(int i{goal}; i != start; i = predecessor[i])
	{
		path.push_back(i);
	}
	path.push_back(start);
	// std::reverse(path.begin(), path.end()); // #include<algorithm>;
	return path;
}

/*
 * Bellman-Ford algorithm, O(VE). Returns the shortest distance from start node to all
 * other nodes. If such a path has arbitrarily low cost -INF is returned, and if
 * there exists no path INF is returned.
 * Predecessor will be updated even if there are negative cycles.
 */
template<typename T>
std::vector<T>& graph<T>::bellman_ford(int start)
{
	int nodes = N;
	std::vector<T> shortest_distance(nodes, INF);
	predecessor = std::vector<int>(nodes, -1);

	shortest_distance[start] = 0;

	// Go through all edges, V-1 times
	// Can add: Break if no change during a whole iteration.

	// Adjacency list version
	for(int i{1}; i <= nodes - 1; ++i)
	{
		for(int j{0}; j < nodes; ++j)
		{
			for(int m{0}; m < AL[j].size(); ++m)
			{
				// If shortest_distance[nd1] == INF then we try to go from a node that cannot (yet) be reached
				// Could otherwise get a problem where INF > INF - 2
				if(shortest_distance[AL[j][m].second] > shortest_distance[j] + AL[j][m].first && shortest_distance[j] != INF)
				{
					shortest_distance[AL[j][m].second] = shortest_distance[j] + AL[j][m].first;
					predecessor[AL[j][m].second] = j;
				}
			}
		}
	}


	std::vector<int> neg_cyc(nodes, 0);

	// Adjacency list version
	// Check for relaxation, if relaxation is possible, then a negative cycle exists.
	for(int j{0}; j < nodes; ++j)
	{
		for(int m{0}; m < AL[j].size(); ++m)
		{
			if(shortest_distance[AL[j][m].second] > shortest_distance[j] + AL[j][m].first && shortest_distance[j] != INF)
			{
				// Negative cycle from this node.
				neg_cyc[j] = 1;
				neg_cyc[AL[j][m].second] = 1;
			}
		}
	}


	// Check predecessors for negative cycles: nodes that can be reached from
	// such nodes have potentially -INF cost to be reached.
	// Find all nodes that are *reachable* from nodes with a negative cycle.
	// Queue all nodes in negative cycles.
	std::queue<int> nc{};
	for(int j{0}; j < nodes; ++j)
	{
		if(neg_cyc[j])
		{
			nc.push(j);
		}
	}

	int nd1{};
	// Mark nodes reachable from nodes in negative cycles
	while(!nc.empty())
	{
		nd1 = nc.front();
		nc.pop();
		int m = AL[nd1].size();
		for(int i{0}; i < m; ++i)
		{
			if(neg_cyc[AL[nd1][i].second])
			{
				continue;
			}
			else
			{
				neg_cyc[AL[nd1][i].second] = 1;
				nc.push(AL[nd1][i].second);
			}
		}

	}

	// Set all marked nodes to have distances -INF.
	for(int i{0}; i < nodes; ++i)
	{
		if(neg_cyc[i])
		{
			shortest_distance[i] = -INF;
		}
	}

	return shortest_distance;
}

/*
 * Floyd-Warshall algorithm, O(V^3). Finds shortest paths
 * between all pairs in a graph (with less than ~450 nodes).
 *
 * TODO: Implement specialization of floyd_warshall for real numbers with error compensation.
 */
template<typename T>
std::vector<std::vector<T>> graph<T>::floyd_warshall()
{
	int nodes = N;

	for(int k{0}; k < nodes; ++k)
	{
		for(int i{0}; i < nodes; ++i)
		{
			for(int j{0}; j < nodes; ++j)
			{
				// Check so that paths used are not infinite.
				if( ( AM[i][j] > AM[i][k] + AM[k][j] ) && ( AM[k][j] < INF ) && ( AM[i][k] < INF ) )
				{
					AM[i][j] = AM[i][k] + AM[k][j];
				}
			}
		}
	}

	// If AM[i][i] is negative then there exists a negative cycle, and i is part of it.
	// The path from this node to itself is then arbitrarily small

	/*
	 * Go through all node paths (again O(V^3)) and check if it is possible
	 * to go from i to j and visits any node k that is part of negative cycle.
	 */
	for(int i{0}; i < nodes; ++i)
	{
		for(int j{0}; j < nodes; ++j)
		{
			for(int k{0}; k < nodes; ++k)
			{
				if(AM[i][k] < INF && AM[k][j] < INF && AM[k][k] < 0)
				{
					AM[i][j] = -INF;
				}
			}
		}
	}

	return AM;
}

/* Max flow algorithm, Edmond-Karp O(V*E^2). Returns the
 * maximum flow of a given graph with sink and source.
 *
 * AM (flow graph) should contain the capacity of edges (back
 * edges are set to zero and increased when flow is sent).
 * AL is initialized to contain forward and (virtual) backward edges.
 *
 * Assumes allowed flow present. (TODO: virtual sink/source
 * for finding an allowed flow.)
 *
 * Capacity constraint, no flow can exceed capacity.
 * Conservation constraint, no vertex can have different out- and in-flow.
 */
template<typename T>
T graph<T>::maxflow(int source, int sink)
{
	// Build AL with no multiedges
	for(int i{0}; i < N; ++i)
	{
		for(int j{i}; j < N; ++j)
		{
			if(AM[i][j] > 0 || AM[j][i] > 0)
			{
				AL[i].push_back( { 0, j } );
				AL[j].push_back( { 0, i } );
			}
		}
	}

	predecessor = std::vector<int>(N, -1);

	T maximum{0};

	while(hasAugmentingPath(source, sink))
	{
		T flow = INF;

		// Find limiting flow by backtracking
		for(int current = sink; current != source; current = predecessor[current])
		{
			int parent = predecessor[current];
			if(flow > AM[parent][current]) // Look at forward capacities
			{
				flow = AM[parent][current];
			}
		}

		for(int current = sink; current != source; current = predecessor[current])
		{
			int parent = predecessor[current];
			AM[parent][current] -= flow; // Decrease forward capacity by sent flow
			AM[current][parent] += flow; // Increase capacity of back edge (Flow sent forward on edge)
		}

		maximum += flow;
	}

	return maximum;
}

/*
 * Finds an augmenting path in a maximum flow algorithm
 * iteration. Returns false if it is impossible to reach
 * the sink node with additional flow.
 */
template<typename T>
bool graph<T>::hasAugmentingPath(int source, int sink)
{
	std::vector<bool> visited(N, false);
	std::queue<int> q{};

	q.push(source);
	predecessor[source] = -1;
	visited[source] = true;

	while(!q.empty())
	{
		int from = q.front();
		q.pop();

		for(int i{0}; i < AL[from].size(); ++i)
		{
			int to = AL[from][i].second;

			// To not visited nodes with forward capacity > 0
			if(!visited[to] && AM[from][to] > 0)
			{
				q.push(to);
				visited[to] = true;
				predecessor[to] = from;

				if(to == sink)
				{
					return true;
				}
			}
		}
	}

	return false;
}

/*
 * Return the adjacency list of the graph object.
 */
template<typename T>
std::vector<std::vector<std::pair<T, int>>> graph<T>::getAL()
{
	return AL;
}

/*
 * Return the adjacency matrix of the graph object.
 */
template<typename T>
std::vector<std::vector<T>> graph<T>::getAM()
{
	return AM;
}

/*
 * Return the value of INF of the graph object.
 */
template<typename T>
T graph<T>::getINF()
{
	return INF;
}

#endif /* GRAPH_H_ */



