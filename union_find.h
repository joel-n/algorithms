/*
 * union_find.h
 *
 *  Created on: 25 feb. 2021
 *      Author: Joel Nilsson
 */

#ifndef UNION_FIND_H_
#define UNION_FIND_H_

#include <vector>

class UF
{
public:

	/*
	 * Union find set with N vertices.
	 */
	UF(int N);

	/* Finds the root of the tree to which a node n belongs.
	 * A node is automatically assigned to have the root node as
	 * a parent when searched for efficiency.
	 */
	int id(int n);

	// Merges two disjoint sets by setting the root of
	// one tree (the smaller) to be the child of the other.
	void U(int a, int b);

private:
	std::vector<int> parent;
	std::vector<int> rank;

};



#endif /* UNION_FIND_H_ */
