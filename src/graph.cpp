#include"graph.h"
#include"dataset.h"
#include"utils.h"
using namespace std;

graph::graph(dataset *ptrndata){
	nodenumber = 0;
	for(int i = 0;i<ptrndata->M;i++){
		string nodein = ptrndata->docs[i]->name;
		
		addnode(nodein);   //doc id is same as the node id, for the modeled docs
		for(int j = 0;j<ptrndata->docs[i]->citationcount;j++){
			string nodeout = ptrndata->docs[i]->citations[j];
			
			addnode(nodeout);
			addlink(nodein, nodeout);
		}
	}
}


void graph::addnode(string node){
	if(atoi(node.c_str())==100){
		printf("add node\n");
	}
	int id = fromnametoid(node);
	
	if(id==-1){
		id = nodenumber;
		nodenumber++;
		
		linkset ll;
//		printf("max size %d \n",ll.max_size());
		nametoid.insert(std::make_pair(atoi(node.c_str()), id));
		idtoname.insert(std::make_pair(id,atoi(node.c_str())));
		inlinks.insert(std::make_pair(id, ll));
		outlinks.insert(std::make_pair(id, ll));

	}
	
}

void graph::addlink(std::string fromLink, std::string toLink){
	int fromid = fromnametoid(fromLink);
	int toid = fromnametoid(toLink);
	if(atoi(toLink.c_str())==3368){
		printf("3368\n");
	}
	addlink(fromid, toid);
}

void graph::addlink(int fromLink, int toLink){
	map<int, linkset>::iterator it1;
	map<int, linkset>::iterator it2;
	
	it1 = outlinks.find(fromLink);
	it2 = inlinks.find(toLink);
	
	it1->second.insert(toLink); //fromLinkµÄoutlink set
	it2->second.insert(fromLink); //toLinkµÄinlink set
	//printf("maxsize %d \n", it2->second.size());
	
}

linkset graph::inLink(int node){
	map<int, linkset>::iterator it;
	it = inlinks.find(node);
	if(it!=inlinks.end())
		return it->second;
	linkset ll;
	return ll;
}

linkset graph::outLink(int node){
	map<int, linkset>::iterator it;
	it = outlinks.find(node);
	if(it!=outlinks.end())
		return it->second;
	linkset ll;
	return ll;
}

int graph::countnode(){
	return nodenumber;
}

int graph::fromnametoid(std::string name){
	map<int, int>::iterator it;
	it = nametoid.find(atoi(name.c_str()));
	if(it!=nametoid.end())
		return it->second;
	return -1;
}

string graph::fromidtoname(int id){
	map<int, int>::iterator it;
	it = idtoname.find(id);
	return utils::int2str(it->second);
}