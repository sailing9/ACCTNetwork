#include"pagerank.h"
#include<cmath>
#include<limits>
#include<set>
using namespace std;

void pagerank::setdampening(double dampening){
	this->dampening = dampening;
}

double pagerank::getdampening(){
	return this->dampening;
}

double pagerank::getPagerankScore(int node){
	map<int, double>::iterator pos;
	pos = this->scores.find(node);
	if(pos!=this->scores.end())
		return pos->second;
	else
		return -1;
}

void pagerank::initializePagerankScore(int node, double value){
	map<int, double>::iterator pos;
	pos = this->scores.find(node);
	if(pos!=this->scores.end()){
		pos->second = value;
	}
}

double pagerank::computeDifference(map<int,double> oldscore, map<int, double> newscore){
	if(oldscore.size()!=newscore.size())
		return -1;
	map<int, double>::iterator itnew;
	map<int, double>::iterator itold;

	double diff = 0.0;
	for(itold=oldscore.begin();itold!=oldscore.end();++itold){
		int key = itold->first;
		double valueold = itold->second;
		double valuenew = newscore.find(key)->second;
		diff+=fabs(valuenew-valueold);
	}
	return diff;

}


void pagerank::computePagerankScore(){
	double diff = numeric_limits<double>::max();
	int n = this->network.countnode();
	map<int, double> newscores;
	while(diff>0.5){ //to be modified the diff should be much smaller
		for(int i = 0;i<n;i++){
			linkset in = (this->network).inLink(i);  //i is the id of the node in the graph
			linkset::iterator it;
			double tmpscore = 0.0;
			/*if(in.size()==0)
				continue;*/
			for(it=in.begin();it!=in.end();++it){
				int nodekey = *it;
				int outdegree = ((this->network).outLink(nodekey)).size();
				if(outdegree!=0)
					tmpscore +=((this->scores).find(nodekey)->second)/(double)outdegree;
			}
			tmpscore = this->dampening*tmpscore + (1-this->dampening)/(double)n;
			newscores.insert(std::make_pair(i,tmpscore));
		}
		diff = this->computeDifference(this->scores, newscores);
		for(int i = 0;i<n;i++){
			this->scores.find(i)->second = newscores.find(i)->second;
		}
		newscores.clear();
	}
}