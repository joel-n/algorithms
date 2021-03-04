/*
 * tree.h
 *
 *  Created on: 26 feb. 2021
 *      Author: Joel Nilsson
 *
 *
 *       NB: T should only be numerical data types.
 */


#ifndef TREE_H_
#define TREE_H_

#include<vector>
#include"union_find.h"

template<typename T>
class Tree
{
public:

	Tree<T>(int n);

	std::vector<std::pair<int, int>> minspantree(std::vector<std::pair<T, std::pair<int, int>>> const& EdgeList);

	T minspantreecost(std::vector<std::pair<T, std::pair<int, int>>> const& EL);

	T get_cost();

private:
	int nodes;
	T cost;
};

/*
 * Initializes tree with n nodes and zero cost.
 * NB: Tree should only be used with numerical data types.
 */
template<typename T>
Tree<T>::Tree(int n) : nodes{n}
{
	cost = 0;
}

/* Kruskals algorithm for minimum spanning tree, O(E*log V). Returns a minimum spanning tree
 * if one exists, otherwise returns an empty list. The cost of the tree is returned
 * by calling get_cost() after calling this function.
 *
 * Using union_find.h to speed up sub-tree merging.
 */
template<typename T>
std::vector<std::pair<int, int>> Tree<T>::minspantree(std::vector<std::pair<T, std::pair<int, int>>> const& EL)
{
	int EDS = EL.size();
	std::vector<std::pair<int, int>> edges{};

	UF set(nodes);

	cost = 0;
	int edgecount{};

	for(int i{0}; i < EDS; ++i)
	{
		int a = set.id(EL[i].second.first);
		int b = set.id(EL[i].second.second);
		if(a != b)
		{
			cost += EL[i].first;

			// Sorting step not necessary if edges are not
			// to be printed in some specific order.
			// Store in correct order (i,j), where i < j.
			if(EL[i].second.first < EL[i].second.second)
			{
				edges.push_back(EL[i].second);
			}
			else
			{
				edges.push_back(std::make_pair(EL[i].second.second, EL[i].second.first));
			}

			// Merge trees
			set.U(EL[i].second.first, EL[i].second.second);
		}

		if(edgecount == nodes-1)
		{
			break;
		}
	}

	int check = set.id(0);
	for(int i{0}; i < nodes; ++i)
	{
		if(set.id(i) != check)
		{
			edges.resize(0);
			break;
		}
	}

	return edges;
}

/* Returns the cost of a minimum spanning tree from Kruskals algorithm, O(E*log V).
 *
 * Function can be used only if there exists a minimum spanning tree. Otherwise
 * use minspantree(), which returns an edgelist and check if the list is empty.
 * The cost can be found by get_cost().
 */
template<typename T>
T Tree<T>::minspantreecost(std::vector<std::pair<T, std::pair<int, int>>> const& EL)
{
	int EDS = EL.size();

	UF set(nodes);

	cost = 0;
	int edgecount{};

	for(int i{0}; i < EDS; ++i)
	{
		int a = set.id(EL[i].second.first);
		int b = set.id(EL[i].second.second);
		if(a != b)
		{
			cost += EL[i].first;
			++edgecount;
			// merge trees
			set.U(EL[i].second.first, EL[i].second.second);
		}

		if(edgecount == nodes-1)
		{
			break;
		}
	}

	return cost;
}

template<typename T>
T Tree<T>::get_cost()
{
	return cost;
}


#endif /* TREE_H_ */
