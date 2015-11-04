#ifndef _GRAPH_H
#define _GRAPH_H

#include<vector>
#include<set>
#include<map>
#include<string>
#include"dataset.h"
using namespace std;

typedef set<int> linkset;

class graph{
public:
	map<int, linkset> inlinks;
	map<int, linkset> outlinks;
	map<int,int> nametoid;
	map<int, int> idtoname;

	int nodenumber;
	
	graph(){
		nodenumber = 0;
	}

	graph(dataset * ptrndata);

	~graph(){
		if(!inlinks.empty()){
			inlinks.clear();
		}
		if(!outlinks.empty()){
			outlinks.clear();
		}
	
	}

	int fromnametoid(string name);
	string fromidtoname(int id);

	void addnode(string nodename);
	void addlink(string fromLink,string toLink);
	void addlink(int fromLink,int toLink);
	linkset inLink(int node);
	linkset outLink(int node);
	int countnode();

};
#endif