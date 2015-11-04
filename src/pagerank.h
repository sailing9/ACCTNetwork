/*
Author:Liu Liu
Date: May 28, 2008 
*/

#ifndef _PAGERANK_H
#define _PAGERANK_H

#include<set>
#include<map>
#include"graph.h"

using namespace std;

class pagerank{
public:
	graph network;
	double dampening;
	map<int, double> scores;

	pagerank(graph network){
		this->network = network;
		int numberofnode = network.countnode();
		double initialscore = (double)1/(double)numberofnode;
		for(int i = 0;i<numberofnode;i++){
			this->scores.insert(std::make_pair(i,initialscore));
		}
	}

	~pagerank(){
	}

	void setdampening(double dampling);
	double getdampening();
	double getPagerankScore(int node);
	void initializePagerankScore(int node, double value);
	void computePagerankScore();
	double computeDifference(map<int, double> oldscore, map<int, double> newscore);

};

#endif