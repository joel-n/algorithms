/*
 * union_find.cpp
 *
 *  Created on: 25 feb. 2021
 *      Author: Joel Nilsson
 */


#include <vector>
#include"union_find.h"

/* Union find object.
 *
 * The operation of finding the node is the root of a tree to which a specific
 * node belongs can be carried out in constant time since the parent of each
 * node is automatically updated to be the root when a query for the root is
 * sent (in fact, all other nodes between the node in question and the root
 * also get updated, see function id()). When trees/sets are merged the smaller
 * tree is added to the larger tree. Merges are also constant time since the
 * find-function is called, and the tree sizes are compared.
 */
UF::UF(int N)
{
	std::vector<int> par{};
	for(int d{0}; d < N; ++d)
	{
		par.push_back(d);
	}
	parent = par;
	rank = std::vector<int>(N, 1);
}


// Finds the root of the tree to which a node n belongs.
// A node is automatically assigned to have the root node as
// a parent when searched for efficiency.
int UF::id(int n){
	if(parent[n] == n)
	{
		return n;
	}
	else
	{
		// Also setting the new parent of n to be the parent of the parent of n
		return parent[n] = id(parent[n]);
	}
}

// Union merges two disjoint sets by setting the root of
// one tree (the smaller) to be the child of the other.
void UF::U(int a, int b){
	a = id(a);
	b = id(b);
	// Do nothing if nodes have the same root
	if(b != a)
	{
		// Low rank root a added to high rank root b
		if(rank[a] <= rank[b])
		{
			parent[a] = b;
			rank[b] += rank[a]; // increase size of root
		}
		else
		{
			parent[b] = a;
			rank[a] += rank[b];
		}
	}
}


