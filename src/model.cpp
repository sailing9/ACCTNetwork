/*
 * Copyright (C) 2007 by
 * 
 * 	Xuan-Hieu Phan
 *	hieuxuan@ecei.tohoku.ac.jp or pxhieu@gmail.com
 * 	Graduate School of Information Sciences
 * 	Tohoku University
 *
 * GibbsLDA++ is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation; either version 2 of the License,
 * or (at your option) any later version.
 *
 * GibbsLDA++ is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GibbsLDA++; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */

/* 
 * References:
 * + The Java code of Gregor Heinrich (gregor@arbylon.net)
 *   http://www.arbylon.net/projects/LdaGibbsSampler.java
 * + "Parameter estimation for text analysis" by Gregor Heinrich
 *   http://www.arbylon.net/publications/text-est.pdf
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <algorithm>
#include <functional>      // For greater<int>( )

#include "constants.h"
#include "strtokenizer.h"
#include "utils.h"
#include "dataset.h"
#include "graph.h"
#include "pagerank.h"
#include "gammaf.h"
#include "model.h"


#define isnan(x) ((x) != (x))
using namespace std;

model::~model() {
	/*
	if (p) {
	delete p;
    }

    if (ptrndata) {
	//delete ptrndata;
    }
    
    if (pnewdata) {
	delete pnewdata;
    }

    if (z) {
	for (int m = 0; m < M; m++) {
	    if (z[m]) {
		delete z[m];
	    }
	}
	delete z;
    }

  
    if (nw) {
	for (int w = 0; w < V; w++) {
	    if (nw[w]) {
		delete nw[w];
	    }
	}
	delete nw;
    }

    if (nd) {
	for (int m = 0; m < M; m++) {
	    if (nd[m]) {
		delete nd[m];
	    }
	}
	delete nd;
    } 
	
    
    if (nwsum) {
	delete nwsum;
    }   
    
    if (ndsum) {
	delete ndsum;
    }
    
    if (theta) {
	for (int m = 0; m < M; m++) {
	    if (theta[m]) {
		delete theta[m];
	    }
	}
	free(theta);
    }
	
    
    if (phi) {
	for (int k = 0; k < K; k++) {
	    if (phi[k]) {
		delete phi[k];
	    }
	}
	free(phi);
    }
	

    // only for inference
    if (newz) {
	for (int m = 0; m < newM; m++) {
	    if (newz[m]) {
		delete newz[m];
	    }
	}
	free(newz);
    }

    
    if (newnw) {
	for (int w = 0; w < newV; w++) {
	    if (newnw[w]) {
		delete newnw[w];
	    }
	}
	free(newnw);
    }

    if (newnd) {
	for (int m = 0; m < newM; m++) {
	    if (newnd[m]) {
		delete newnd[m];
	    }
	}
	free(newnd);
    } 
    
    if (newnwsum) {
	delete newnwsum;
    }   
    
    if (newndsum) {
	delete newndsum;
    }
    
    if (newtheta) {
	for (int m = 0; m < newM; m++) {
	    if (newtheta[m]) {
		delete newtheta[m];
	    }
	}
	free(newtheta);
    }
    
    if (newphi) {
	for (int k = 0; k < K; k++) {
	    if (newphi[k]) {
		delete newphi[k];
	    }
	}
	free(newphi);
    }
	*/
}

void model::set_default_values() {
    wordmapfile = "wordmap.txt";
	docmapfile = "docmap.txt";
citationmapfile = "citationmap.txt"; 
    trainlogfile = "trainlog.txt";
    tassign_suffix = ".tassign";
    theta_suffix = ".theta";
    phi_suffix = ".phi";
    others_suffix = ".others";
    twords_suffix = ".twords";
    
    dir = "./";
    dfile = "trndocs.dat";
    model_name = "model-final";    
    model_status = MODEL_STATUS_UNKNOWN;
    
    ptrndata = NULL;
    pnewdata = NULL;
    
		 
    M = 0;
    V = 0;
    K = 100;
    alpha = 50.0 / K;
    beta = 0.1;
    niters = 2000;
    liter = 0;
    savestep = 200;    
    twords = 0;
    withrawstrs = 0;
    
	alphaarray = NULL;
    p = NULL;
    z = NULL;
    nw = NULL;
    nd = NULL;
    nwsum = NULL;
    ndsum = NULL;
    theta = NULL;
    phi = NULL;
    
    newM = 0;
    newV = 0;
    newz = NULL;
    newnw = NULL;
    newnd = NULL;
    newnwsum = NULL;
    newndsum = NULL;
    newtheta = NULL;
    newphi = NULL;
}

int model::parse_args(int argc, char ** argv) {
    return utils::parse_args(argc, argv, this);
}

int model::init(int argc, char ** argv) {
    // call parse_args
    if (parse_args(argc, argv)) {
	return 1;
    }
    
    if (model_status == MODEL_STATUS_EST) {
	// estimating the model from scratch
	if (init_est()) {
	    return 1;
	}
	
    } else if (model_status == MODEL_STATUS_ESTC) {
	// estimating the model from a previously estimated one
	if (init_estc()) {
	    return 1;
	}
	
    } else if (model_status == MODEL_STATUS_INF) {
	// do inference
	if (init_inf()) {
	    return 1;
	}
    }else if(model_status == MODEL_STATUS_EVALUATE){ //add by liu liu
		if(init_evaluate()){
			return 1;
		}
	
	}
    
    return 0;
}

int model::load_model(string model_name) {
    int i, j;
    
    string filename = dir + model_name + tassign_suffix;
    FILE * fin = fopen(filename.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %d to load model!\n", filename.c_str());
	return 1;
    }
    
    char buff[BUFF_SIZE_LONG];
    string line;

    // allocate memory for z and ptrndata
    z = new int*[M];
    ptrndata = new dataset(M);
    ptrndata->V = V;

    for (i = 0; i < M; i++) {
	char * pointer = fgets(buff, BUFF_SIZE_LONG, fin);
	if (!pointer) {
	    printf("Invalid word-topic assignment file, check the number of docs!\n");
	    return 1;
	}
	
	line = buff;
	strtokenizer strtok(line, " \t\r\n");
	int length = strtok.count_tokens();
	
	vector<int> words;
	vector<int> topics;
	for (j = 0; j < length; j++) {
	    string token = strtok.token(j);
    
	    strtokenizer tok(token, ":");
	    if (tok.count_tokens() != 2) {
		printf("Invalid word-topic assignment line!\n");
		return 1;
	    }
	    
	    words.push_back(atoi(tok.token(0).c_str()));
	    topics.push_back(atoi(tok.token(1).c_str()));
	}
	
	// allocate and add new document to the corpus
	document * pdoc = new document(words);
	ptrndata->add_doc(pdoc, i);
	
	// assign values for z
	z[i] = new int[topics.size()];
	for (j = 0; j < topics.size(); j++) {
	    z[i][j] = topics[j];
	}
    }   
    
    fclose(fin);
    
    return 0;
}

int model::save_model(string model_name) {
    if (save_model_tassign(dir + model_name + tassign_suffix)) {
	return 1;
    }
    
    if (save_model_others(dir + model_name + others_suffix)) {
	return 1;
    }
    
    if (save_model_theta(dir + model_name + theta_suffix)) {
	return 1;
    }

	if(save_model_alpha(dir + model_name + "alpha")){
		return 1;
	}
    
    if (save_model_phi(dir + model_name + phi_suffix)) {
	return 1;
    }
    
	if(save_model_likelihood(dir + model_name + "likelihood")) {
		return 1;
	}
    if (twords > 0) {
	if (save_model_twords(dir + model_name + twords_suffix)) {
	    return 1;
	}
    }
    
    return 0;
}

int model::save_model_tassign(string filename) {
    int i, j;
    
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < ptrndata->M; i++) {    
	for (j = 0; j < ptrndata->docs[i]->length; j++) {
	    fprintf(fout, "%d:%d ", ptrndata->docs[i]->words[j], z[i][j]);
	}
	fprintf(fout, "\n");
    }

    fclose(fout);
    
    return 0;
}

int model::save_model_theta(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < M; i++) {
	for (int j = 0; j < K; j++) {
	    fprintf(fout, "%f ", theta[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_model_alpha(string filename){
	FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < M; i++) {
	for (int j = 0; j < K; j++) {
	    fprintf(fout, "%f	", alphaarray[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_model_likelihood(string filename){
	

	FILE * fout = fopen(filename.c_str(), "a");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

	double likelihood = compute_likelihood();
	fprintf(fout, "Topic = %d\n", K);
	fprintf(fout, "%f ", likelihood);
	fprintf(fout, "\n");
	  fclose(fout);
	return 0; 
}

int model::save_model_phi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < K; i++) {
	for (int j = 0; j < V; j++) {
	    fprintf(fout, "%f ", phi[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);    
    
    return 0;
}

int model::save_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    fprintf(fout, "alpha=%f\n", alpha);
    fprintf(fout, "beta=%f\n", beta);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "ndocs=%d\n", M);
    fprintf(fout, "nwords=%d\n", V);
    fprintf(fout, "liter=%d\n", liter);
    
    fclose(fout);    
    
    return 0;
}

int model::save_model_twords(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    if (twords > V) {
	twords = V;
    }
    mapid2word::iterator it;
    
    for (int k = 0; k < K; k++) {
		vector<pair<int, double> > words_probs;
		pair<int, double> word_prob;
		for (int w = 0; w < V; w++) {
			word_prob.first = w;
			word_prob.second = phi[k][w];
			words_probs.push_back(word_prob);
		}
    
        // quick sort to sort word-topic probability
		//utils::quicksort(words_probs, 0, words_probs.size() - 1);
		sort(words_probs.begin(), words_probs.end(), PairGreater);
		
		fprintf(fout, "Topic %dth:\n", k);
		for (int i = 0; i < twords; i++) {
			it = id2word.find(words_probs[i].first);
			if (it != id2word.end()) {
			fprintf(fout, "\t%s   %f\n", (it->second).c_str(), words_probs[i].second);
			}
		}
    }
    
    fclose(fout);    
    
    return 0;    
}

int model::save_inf_model(string model_name) {
    if (save_inf_model_tassign(dir + model_name + tassign_suffix)) {
	return 1;
    }
    
    if (save_inf_model_others(dir + model_name + others_suffix)) {
	return 1;
    }
    
    if (save_inf_model_newtheta(dir + model_name + theta_suffix)) {
	return 1;
    }
    
    if (save_inf_model_newphi(dir + model_name + phi_suffix)) {
	return 1;
    }

    if (twords > 0) {
	if (save_inf_model_twords(dir + model_name + twords_suffix)) {
	    return 1;
	}
    }
    
	compute_and_save_perplexity("perplexity.txt");
    return 0;
}

int model::compute_and_save_perplexity(string filename){
	double perplexity = compute_perplexity();
	FILE * fout = fopen(filename.c_str(), "a");
    if (!fout) {
		printf("Cannot open file %s to save!\n", filename.c_str());
		return 1;
    }

	fprintf(fout, "K=%d\n",K);
	fprintf(fout, "perplexity=%f\n", perplexity);
	fprintf(fout, "\n");
	fclose(fout);
	return 0;
}

int model::save_inf_model_tassign(string filename) {
    int i, j;
    
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    // wirte docs with topic assignments for words
    for (i = 0; i < pnewdata->M; i++) {    
	for (j = 0; j < pnewdata->docs[i]->length; j++) {
	    fprintf(fout, "%d:%d ", pnewdata->docs[i]->words[j], newz[i][j]);
	}
	fprintf(fout, "\n");
    }

    fclose(fout);
    
    return 0;
}

int model::save_inf_model_newtheta(string filename) {
    int i, j;

    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (i = 0; i < newM; i++) {
	for (j = 0; j < K; j++) {
	    fprintf(fout, "%f ", newtheta[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);
    
    return 0;
}

int model::save_inf_model_newphi(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    for (int i = 0; i < K; i++) {
	for (int j = 0; j < newV; j++) {
	    fprintf(fout, "%f ", newphi[i][j]);
	}
	fprintf(fout, "\n");
    }
    
    fclose(fout);    
    
    return 0;
}

int model::save_inf_model_others(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }

    fprintf(fout, "alpha=%f\n", alpha);
    fprintf(fout, "beta=%f\n", beta);
    fprintf(fout, "ntopics=%d\n", K);
    fprintf(fout, "ndocs=%d\n", newM);
    fprintf(fout, "nwords=%d\n", newV);
    fprintf(fout, "liter=%d\n", inf_liter);
    
    fclose(fout);    
    
    return 0;
}

int model::save_inf_model_twords(string filename) {
    FILE * fout = fopen(filename.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to save!\n", filename.c_str());
	return 1;
    }
    
    if (twords > newV) {
	twords = newV;
    }
    mapid2word::iterator it;
    map<int, int>::iterator _it;
    
    for (int k = 0; k < K; k++) {
	vector<pair<int, double> > words_probs;
	pair<int, double> word_prob;
	for (int w = 0; w < newV; w++) {
	    word_prob.first = w;
	    word_prob.second = newphi[k][w];
	    words_probs.push_back(word_prob);
	}
    
        // quick sort to sort word-topic probability
	//utils::quicksort(words_probs, 0, words_probs.size() - 1);
	sort(words_probs.begin(), words_probs.end(), PairGreater);
	
	fprintf(fout, "Topic %dth:\n", k);
	for (int i = 0; i < twords; i++) {
	    _it = pnewdata->_id2id.find(words_probs[i].first);
	    if (_it == pnewdata->_id2id.end()) {
		continue;
	    }
	    it = id2word.find(_it->second);
	    if (it != id2word.end()) {
		fprintf(fout, "\t%s   %f\n", (it->second).c_str(), words_probs[i].second);
	    }
	}
    }
    
    fclose(fout);    
    
    return 0;    
}


int model::init_est() {
    int m, n, w, k;

    p = new double[K];

    // + read training data
    ptrndata = new dataset;
    if (ptrndata->read_trndata(dir + dfile, dir + wordmapfile, dir+docmapfile, dir+citationmapfile)) {
        printf("Fail to read training data!\n");
        return 1;
    }
		
    // + allocate memory and assign values for variables
    M = ptrndata->M;
    V = ptrndata->V;
    // K: from command line or default value
    // alpha, beta: from command line or default values
    // niters, savestep: from command line or default values

	alphaarray = new double*[M];
	for(m = 0;m<M;m++){
		alphaarray[m] = new double[K];
		for(k = 0;k<K;k++){
			alphaarray[m][k] = alpha;
		}
	}

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }
	
    nwsum = new int[K];
    for (k = 0; k < K; k++) {
	nwsum[k] = 0;
    }
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	ndsum[m] = 0;
    }

    //srandom(time(0)); // initialize for random number generation
	srand(time(0));
    z = new int*[M];
    for (m = 0; m < ptrndata->M; m++) {
	int N = ptrndata->docs[m]->length;
	z[m] = new int[N];
	
        // initialize for z
        for (n = 0; n < N; n++) {
    	    //int topic = (int)(((double)random() / RAND_MAX) * K);
			int topic = (int)(((double)rand() / RAND_MAX) * K);
    	    z[m][n] = topic;
    	    
    	    // number of instances of word i assigned to topic j
    	    nw[ptrndata->docs[m]->words[n]][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        } 
        // total number of words in document i
        ndsum[m] = N;
    }
    
    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
    }
	
    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }

    return 0;
}

int model::init_estc() {
    // estimating the model from a previously estimated one
    int m, n, w, k;

    p = new double[K];

    // load moel, i.e., read z and ptrndata
    if (load_model(model_name)) {
	printf("Fail to load word-topic assignmetn file of the model!\n");
	return 1;
    }

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }
	
    nwsum = new int[K];
    for (k = 0; k < K; k++) {
	nwsum[k] = 0;
    }
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	ndsum[m] = 0;
    }

    for (m = 0; m < ptrndata->M; m++) {
	int N = ptrndata->docs[m]->length;

	// assign values for nw, nd, nwsum, and ndsum	
        for (n = 0; n < N; n++) {
    	    int w = ptrndata->docs[m]->words[n];
    	    int topic = z[m][n];
    	    
    	    // number of instances of word i assigned to topic j
    	    nw[w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        } 
        // total number of words in document i
        ndsum[m] = N;      
    }
	
    theta = new double*[M];
    for (m = 0; m < M; m++) {
        theta[m] = new double[K];
    }
	
    phi = new double*[K];
    for (k = 0; k < K; k++) {
        phi[k] = new double[V];
    }    

    return 0;        
}

void model::estimate() {
    if (twords > 0) {
	// print out top words per topic
	dataset::read_wordmap(dir + wordmapfile, &id2word);
	dataset::read_id2citationmap(dir + citationmapfile, &id2citation);
    }

	//add by liu liu
	//step 1
	//build the graph and compute the pagerank
	graph directedGraph = graph(ptrndata);
	pagerank directedPagerank = pagerank(directedGraph);
	directedPagerank.setdampening(0.85);
	directedPagerank.computePagerankScore();

	/*double * pagerankscore = new double[M];
	for(int i = 0;i<M;i++){
		pagerankscore[i] = directedPagerank.scores.find(i)->second;
	}*/
	//数组pagerankscore保存了各个点的pagerank scores

    printf("Sampling %d iterations!\n", niters);

	//step 2: Gibbs sampling 200 iterations

	printf("Sampling the first 200 steps");

	 for (int i = 0; i<200;i++) {
		printf("Iteration %d ...\n", i);
	
		// for all z_i
		for (int m = 0; m < M; m++) {
			for (int n = 0; n < ptrndata->docs[m]->length; n++) {
				
			// (z_i = z[m][n])
			// sample from p(z_i|z_-i, w)
			int topic = sampling(m, n);
			z[m][n] = topic;
			}
		}
	
    }

	//step 3: begin the E(Gibbes samping)M(smooth the objective function) step:

	 int gibbsperstep = 100;  //这个数值是可以自己设的
	 int smoothperstep = 1;  //这个数值也是可以自己设的
	 int emstep = niters/gibbsperstep; //保证gibbs sampling的次数和输入的要求基本相同

	 for(int em = 0;em<emstep;em++){
		  printf("EM step %d ...\n", em);
		//first: Gibbs sampling
		 for(int gibbs = 0;gibbs<gibbsperstep;gibbs++){
			 printf("Iteration %d ...\n", gibbs);
			 for (int m = 0; m < M; m++) {
				for (int n = 0; n < ptrndata->docs[m]->length; n++) {
					
				// (z_i = z[m][n])
				// sample from p(z_i|z_-i, w)
				int topic = sampling(m, n);
				z[m][n] = topic;
				}
			}
		 }

		//second: compute the doc-topic distribution
		 compute_theta();

		//third: update distribution and alpha
		  printf("Smooth -----------%d ...\n", em);
		 double distributionfactor = 0.5;
		 double alphafactor = 0.1;
		 double objectivefactor = 0.1;
		 double newobject = 0.0;
		 double oldobject = 0.0;
		 newobject = compute_objective(objectivefactor, directedPagerank);
		 for(int smooth = 0;smooth<smoothperstep;smooth++){
			oldobject = newobject;
			smooth_distribution(distributionfactor,directedPagerank);
			smooth_alpha(alphafactor, directedPagerank);
			newobject = compute_objective(objectivefactor, directedPagerank);
			if(newobject<oldobject)
				break;
		 }
		 
	 }

//old code for gibbs sampling
	 /*
	 int last_iter = liter;
    for (liter = last_iter + 1; liter <= niters + last_iter; liter++) {
		printf("Iteration %d ...\n", liter);
	
		// for all z_i
		for (int m = 0; m < M; m++) {
			for (int n = 0; n < ptrndata->docs[m]->length; n++) {
				
			// (z_i = z[m][n])
			// sample from p(z_i|z_-i, w)
			int topic = sampling(m, n);
			z[m][n] = topic;
			}
		}
	
		if (savestep > 0) {
			if (liter % savestep == 0) {
			// saving the model
			printf("Saving the model at iteration %d ...\n", liter);
			compute_theta();
			compute_phi();
			save_model(utils::generate_model_name(liter));
			}
		}
    }
    */
    printf("Gibbs sampling completed!\n");
    printf("Saving the final model!\n");
    compute_theta();
    compute_phi();
    liter--;
   save_model(utils::generate_model_name(-1));
}

//to be modified with alpha
int model::sampling(int m, int n) {
    // remove z_i from the count variables
    int topic = z[m][n];
    int w = ptrndata->docs[m]->words[n];
    nw[w][topic] -= 1;
    nd[m][topic] -= 1;
    nwsum[topic] -= 1;
    ndsum[m] -= 1;

    double Vbeta = V * beta;
    //double Kalpha = K * alpha;
	double Kalpha = 0.0;
	for(int i = 0;i<K;i++){
		Kalpha += alphaarray[m][i];
	}
    // do multinomial sampling via cumulative method
    for (int k = 0; k < K; k++) {
	p[k] = (nw[w][k] + beta) / (nwsum[k] + Vbeta) *
		    (nd[m][k] + alphaarray[m][k]) / (ndsum[m] + Kalpha);
    }
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
	p[k] += p[k - 1];
    }
    // scaled sample because of unnormalized p[]
    //double u = ((double)random() / RAND_MAX) * p[K - 1];
	double u = ((double)rand() / RAND_MAX) * p[K - 1];
    
    for (topic = 0; topic < K; topic++) {
	if (p[topic] > u) {
	    break;
	}
    }
    
    // add newly estimated z_i to count variables
    nw[w][topic] += 1;
    nd[m][topic] += 1;
    nwsum[topic] += 1;
    ndsum[m] += 1;    
    
    return topic;
}

//to be modified with alpha 
void model::compute_theta() {
    for (int m = 0; m < M; m++) {
		double Kalpha = 0.0;
		for(int k = 0;k<K;k++){
			Kalpha +=alphaarray[m][k];
		}
	for (int k = 0; k < K; k++) {
	    theta[m][k] = (nd[m][k] + alphaarray[m][k]) / (ndsum[m] + Kalpha);
	}
    }
}

void model::compute_phi() {
    for (int k = 0; k < K; k++) {
	for (int w = 0; w < V; w++) {
	    phi[k][w] = (nw[w][k] + beta) / (nwsum[k] + V * beta);
	}
    }
}

int model::init_inf() {
    // estimating the model from a previously estimated one
    int m, n, w, k;

    p = new double[K];

    // load moel, i.e., read z and ptrndata
    if (load_model(model_name)) {
	printf("Fail to load word-topic assignmetn file of the model!\n");
	return 1;
    }

    nw = new int*[V];
    for (w = 0; w < V; w++) {
        nw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    nw[w][k] = 0;
        }
    }
	
    nd = new int*[M];
    for (m = 0; m < M; m++) {
        nd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    nd[m][k] = 0;
        }
    }
	
    nwsum = new int[K];
    for (k = 0; k < K; k++) {
	nwsum[k] = 0;
    }
    
    ndsum = new int[M];
    for (m = 0; m < M; m++) {
	ndsum[m] = 0;
    }

    for (m = 0; m < ptrndata->M; m++) {
	int N = ptrndata->docs[m]->length;

	// assign values for nw, nd, nwsum, and ndsum	
        for (n = 0; n < N; n++) {
    	    int w = ptrndata->docs[m]->words[n];
    	    int topic = z[m][n];
    	    
    	    // number of instances of word i assigned to topic j
    	    nw[w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    nd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    nwsum[topic] += 1;
        } 
        // total number of words in document i
        ndsum[m] = N;      
    }
    
    // read new data for inference
    pnewdata = new dataset;
    if (withrawstrs) {
	if (pnewdata->read_newdata_withrawstrs(dir + dfile, dir + wordmapfile)) {
    	    printf("Fail to read new data!\n");
    	    return 1;
	}    
    } else {
	if (pnewdata->read_newdata(dir + dfile, dir + wordmapfile, dir+citationmapfile)) {
    	    printf("Fail to read new data!\n");
    	    return 1;
	}    
    }
    
    newM = pnewdata->M;
    newV = pnewdata->V;
    
    newnw = new int*[newV];
    for (w = 0; w < newV; w++) {
        newnw[w] = new int[K];
        for (k = 0; k < K; k++) {
    	    newnw[w][k] = 0;
        }
    }
	
    newnd = new int*[newM];
    for (m = 0; m < newM; m++) {
        newnd[m] = new int[K];
        for (k = 0; k < K; k++) {
    	    newnd[m][k] = 0;
        }
    }
	
    newnwsum = new int[K];
    for (k = 0; k < K; k++) {
	newnwsum[k] = 0;
    }
    
    newndsum = new int[newM];
    for (m = 0; m < newM; m++) {
	newndsum[m] = 0;
    }

    //srandom(time(0)); // initialize for random number generation
	srand(time(0)); // initialize for random number generation
    newz = new int*[newM];
    for (m = 0; m < pnewdata->M; m++) {
	int N = pnewdata->docs[m]->length;
	newz[m] = new int[N];

	// assign values for nw, nd, nwsum, and ndsum	
        for (n = 0; n < N; n++) {
    	    int w = pnewdata->docs[m]->words[n];
    	    int _w = pnewdata->_docs[m]->words[n];
    	    //int topic = (int)(((double)random() / RAND_MAX) * K);
			int topic = (int)(((double)rand() / RAND_MAX) * K);
    	    newz[m][n] = topic;
    	    
    	    // number of instances of word i assigned to topic j
    	    newnw[_w][topic] += 1;
    	    // number of words in document i assigned to topic j
    	    newnd[m][topic] += 1;
    	    // total number of words assigned to topic j
    	    newnwsum[topic] += 1;
        } 
        // total number of words in document i
        newndsum[m] = N;      
    }    
    
    newtheta = new double*[newM];
    for (m = 0; m < newM; m++) {
        newtheta[m] = new double[K];
    }
	
    newphi = new double*[K];
    for (k = 0; k < K; k++) {
        newphi[k] = new double[newV];
    }    
    
    return 0;        
}

int model::load_model_evaluate(std::string model_name){
	int i=0, j;
    
	
	
	string filename_phi = dir + model_name +phi_suffix;
	string filename_theta = dir + model_name + theta_suffix;
  
	
	FILE * fin_phi = fopen(filename_phi.c_str(), "r");
	FILE * fin_theta = fopen(filename_theta.c_str(), "r");
	

	phi = new double*[K];
	for(int k = 0;k<K;k++){
		phi[k] = new double[V];
	}
	
	theta = new double*[M];
	for(int m=0;m<M;m++){
		theta[m] = new double[K];
	}

	char buff[BUFF_SIZE_LONG];
    string line;
	int length;


	if (!fin_phi) {
		printf("Cannot open file %d to load model!\n", filename_phi.c_str());
		return 1;
    }
    
	for(int k = 0; k<K;k++){
		if(fgets(buff, BUFF_SIZE_LONG - 1, fin_phi)==NULL)
			break;
				
		line = buff;
		line = utils::trimstring(line);
		strtokenizer strtok(line, " \t\r\n");
		length = strtok.count_tokens();

		if(length!=V){
			printf("Wrong with the file %d to load model!\n", filename_phi.c_str());
			return 1;
		}
		for(j=0; j<length; j++){
			phi[k][j] = atof(strtok.token(j).c_str());
		}
	}
	fclose(fin_phi);

	if (!fin_theta) {
		printf("Cannot open file %d to load model!\n", filename_theta.c_str());
		return 1;
    }
    
	for(int m = 0; m<M;m++){
		if(fgets(buff, BUFF_SIZE_LONG - 1, fin_theta)==NULL)
			break;
				
		line = buff;
		line = utils::trimstring(line);
		strtokenizer strtok(line, " \t\r\n");
		length = strtok.count_tokens();

		if(length!=K){
			printf("Wrong with the file %d to load model!\n", filename_theta.c_str());
			return 1;
		}
		for(j=0; j<length; j++){
			theta[m][j] = atof(strtok.token(j).c_str());
		}
		
	}
	fclose(fin_theta);


    return 0;
}

int model::init_evaluate(){
	    // estimating the model from a previously estimated one
    int m, n, w, k, a, c, r, t;

    // load moel, i.e., read probability distribution
    if (load_model_evaluate(model_name)) {
		printf("Fail to load probability distribution file of the model!\n");
		return 1;
    }

	
    // read new data for evaluation
    pnewdata = new dataset;
    if (withrawstrs) {
	if (pnewdata->read_newdata_withrawstrs(dir + dfile, dir + wordmapfile)) {
    	    printf("Fail to read new data!\n");
    	    return 1;
	}    
    } else {
	if (pnewdata->read_newdata(dir + dfile, dir + wordmapfile, dir+citationmapfile)) {
    	    printf("Fail to read new data!\n");
    	    return 1;
	}    
    }
    
    newM = pnewdata->M;
    newV = pnewdata->V;
   return 0;
}


void model::evaluate(){
	map<int, string> id2docname;
	map<string, int> docname2id;
	if (twords > 0) {
		// print out top words per topic
		dataset::read_wordmap(dir + wordmapfile, &id2word);
		dataset::read_id2docmap(dir+docmapfile, &id2docname);
		dataset::read_doc2idmap(dir+docmapfile, &docname2id );
		dataset::read_id2citationmap(dir+citationmapfile, &id2citation); //add by liu liu
	}
	map<int, string>::iterator id2nameit;
	mapid2word::iterator wit;
	
	mapword2id word2id;
	mapword2id::iterator w2idit;

	dataset::read_wordmap(dir+wordmapfile, &word2id);

	double * rank = new double[3];

	double p5c = 0.0;
	double p5r = 0.0;
	double p10c = 0.0;
	double p10r = 0.0;
	double rprec = 0.0;
	double rprer = 0.0;
	double meanapc = 0.0;
	double meanapr = 0.0;
    
	double recall5r = 0.0;
	 
	double recall10r = 0.0;
	 
	double mrrr = 0.0;
 
	double bpref10r = 0.0;

	double *rklr = new double[3];

	

	for(int i = 0;i<3;i++){
		
		rklr[i] = 0.0;
		rank[i] = 0.0;
	}
		

	string filecitation = "citation_result.txt";
	FILE * foutcitation = fopen(filecitation.c_str(), "w");
	
	 if (!foutcitation) {
		printf("Cannot open file %s to save!\n", filecitation.c_str());
		return;
    }

	 int R = id2citation.size();
	double p5tmp, p10tmp, rpretmp, meanaptmp, recall5tmp, recall10tmp, bpref10tmp;
		for(int m = 0;m<newM;m++){

			
			//conf_evaluation
			vector<pair<string, double>> doc_probs;
			pair<string, double> doc_prob;

			for(int c = 0;c<R; c++){
				double currentprob = 1.0;
				string docname = utils::int2str(id2citation.find(c)->second);
				int citation = docname2id.find(docname)->second;
				for(int w = 0;w<pnewdata->docs[m]->length;w++){
					double probability = 0.0;
					for(int k = 0;k<K;k++){
						probability +=theta[citation][k]*phi[k][pnewdata->docs[m]->words[w]];
					}
					currentprob *=probability;
				}
				//string name = id2docname.find(c)->second;
				doc_prob.first = docname;
				doc_prob.second = currentprob;
				doc_probs.push_back(doc_prob);

			}

			std::sort(doc_probs.begin(), doc_probs.end(), PairStringGreater);


			//output the result to a file
			fprintf(foutcitation, "============================================================\n");
			fprintf(foutcitation, "Doc %d\n", m);
			//fprintf(foutcitation, "Content: %s \n", pnewdata->docs[m]->title.c_str());
			fprintf(foutcitation, "Actural Citations : \n");
			for(int i = 0;i<pnewdata->docs[m]->citationcount;i++){
				fprintf(foutcitation, "%s \n", pnewdata->docs[m]->citations[i].c_str());
				printf("%s \n", pnewdata->docs[m]->citations[i].c_str());
			}			
			fprintf(foutcitation, "\n");

			
			p5tmp = utils::calcu_precisiontop(5, doc_probs, pnewdata->docs[m]->citations, pnewdata->docs[m]->citationcount);
			p10tmp = utils::calcu_precisiontop(10, doc_probs, pnewdata->docs[m]->citations, pnewdata->docs[m]->citationcount);
			rpretmp = utils::calcu_rpre(doc_probs, pnewdata->docs[m]->citations, pnewdata->docs[m]->citationcount);
			meanaptmp = utils::calcu_map(doc_probs, pnewdata->docs[m]->citations, pnewdata->docs[m]->citationcount);
			recall5tmp = utils::calcu_recalltop(5,doc_probs,pnewdata->docs[m]->citations,pnewdata->docs[m]->citationcount);
			recall10tmp =utils::calcu_recalltop(10,doc_probs,pnewdata->docs[m]->citations,pnewdata->docs[m]->citationcount);
			bpref10tmp = utils::calcu_bpref10(doc_probs, pnewdata->docs[m]->citations, pnewdata->docs[m]->citationcount);

			p5r +=p5tmp/newM;
			p10r +=p10tmp/newM;
			rprer +=rpretmp/newM;
			meanapr +=meanaptmp/newM;
			recall5r +=recall5tmp/newM;
			recall10r +=recall10tmp/newM;
			//mrrr[i] += utils::calcu_mrr(citation_probs[i],pnewdata->docs[m]->citations, pnewdata->docs[m]->citationcount)/newM;
			bpref10r += bpref10tmp/newM;

			rank= utils::calcu_rkl(doc_probs, pnewdata->docs[m]->citations, pnewdata->docs[m]->citationcount);
			for(int j = 0;j<3;j++)
				rklr[j] +=rank[j]/(double)newM;

			fprintf(foutcitation, "P@5=%f   P@10=%f  Rpre=%f  MAP=%f  R@5=%f  R@10=%f  Bpref10=%f  Rank_Best=%f  Rank_Median=%f  Rank_Last=%f \n", p5tmp, p10tmp, rpretmp,  meanaptmp,recall5tmp, recall10tmp, bpref10tmp, rank[0], rank[1], rank[2]);
			//fprintf(foutcitation, "Top10 predict citations: \n");
			for(int j = 0;j<10;j++){
				fprintf(foutcitation, "%s \n", doc_probs[j].first.c_str());
			}
			fprintf(foutcitation, "\n");
	

		}

		string filename = "evaluation.txt";
		FILE * fout = fopen(filename.c_str(), "a");
		if (!fout) {
			printf("Cannot open file %s to save!\n", filename.c_str());
			return;
		}

		fprintf(fout, "=========================================================\n");
		fprintf(fout, "LDA Model \n");
		fprintf(fout, "策略一: no inference. train conf&citation&author \n");
		fprintf(fout, "Citation prediction \n");
		fprintf(fout, "%f\t", p5r);
		fprintf(fout, "%f\t", p10r);
		fprintf(fout, "%f\t", rprer);
		fprintf(fout, "%f\t", meanapr);
		fprintf(fout, "%f\t", recall5r);
		fprintf(fout, "%f\t", recall10r);
		fprintf(fout, "%f\t", bpref10r);
		fprintf(fout, "%f\t", rklr[0]);
		fprintf(fout, "%f\t", rklr[1]);
		fprintf(fout, "%f\t", rklr[2]);
		fprintf(fout, "\n");
		fprintf(fout, "\n");

		fclose(fout);
		fclose(foutcitation);
	
}


//

void model::inference() {
    if (twords > 0) {
	// print out top words per topic
	dataset::read_wordmap(dir + wordmapfile, &id2word);
dataset::read_id2citationmap(dir+citationmapfile, &id2citation); //add by liu liu
    }

    printf("Sampling %d iterations for inference!\n", niters);
    
    for (inf_liter = 1; inf_liter <= niters; inf_liter++) {
	printf("Iteration %d ...\n", inf_liter);
	
	// for all newz_i
		for (int m = 0; m < newM; m++) {
			for (int n = 0; n < pnewdata->docs[m]->length; n++) {
			// (newz_i = newz[m][n])
			// sample from p(z_i|z_-i, w)
			int topic = inf_sampling(m, n);
			newz[m][n] = topic;
			}
		}
    }
    
    printf("Gibbs sampling for inference completed!\n");
    printf("Saving the inference outputs!\n");
    compute_newtheta();
    compute_newphi();
    inf_liter--;
    save_inf_model(dfile);
}

int model::inf_sampling(int m, int n) {
    // remove z_i from the count variables
    int topic = newz[m][n];
    int w = pnewdata->docs[m]->words[n];
    int _w = pnewdata->_docs[m]->words[n];
    newnw[_w][topic] -= 1;
    newnd[m][topic] -= 1;
    newnwsum[topic] -= 1;
    newndsum[m] -= 1;
    
    double Vbeta = V * beta;
    double Kalpha = K * alpha;
    // do multinomial sampling via cumulative method
    for (int k = 0; k < K; k++) {
		p[k] = (nw[w][k] + newnw[_w][k] + beta) / (nwsum[k] + newnwsum[k] + Vbeta) *
		    (newnd[m][k] + alpha) / (newndsum[m] + Kalpha);
    }
    // cumulate multinomial parameters
    for (int k = 1; k < K; k++) {
		p[k] += p[k - 1];
    }
    // scaled sample because of unnormalized p[]
    //double u = ((double)random() / RAND_MAX) * p[K - 1];
	double u = ((double)rand() / RAND_MAX) * p[K - 1];
    
    for (topic = 0; topic < K; topic++) {
	if (p[topic] > u) {
	    break;
	}
    }
    
    // add newly estimated z_i to count variables
    newnw[_w][topic] += 1;
    newnd[m][topic] += 1;
    newnwsum[topic] += 1;
    newndsum[m] += 1;    
    
    return topic;
}

void model::compute_newtheta() {
    for (int m = 0; m < newM; m++) {
	for (int k = 0; k < K; k++) {
	    newtheta[m][k] = (newnd[m][k] + alpha) / (newndsum[m] + K * alpha);
	}
    }
}

void model::compute_newphi() {
    map<int, int>::iterator it;
    for (int k = 0; k < K; k++) {
	for (int w = 0; w < newV; w++) {
	    it = pnewdata->_id2id.find(w);
	    if (it != pnewdata->_id2id.end()) {
		newphi[k][w] = (nw[it->second][k] + newnw[w][k] + beta) / (nwsum[k] + newnwsum[k] + V * beta);
	    }
	}
    }
}

double model::compute_perplexity(){
	double perplexity = 0.0;
	double denominator = 0.0;
	double numerator = 0.0;
	double perp = 0.0;

	for(int m=0; m<newM; m++){
		int Nd = pnewdata->docs[m]->length;		

		perp = 0.0;
		double p = 0.0;
		for(int n=0; n<Nd; n++){
			int w = pnewdata->_docs[m]->words[n];
			for(int k=0; k<K; k++){				
				//p += newphi[k][w];
				p += newphi[k][w] * newtheta[m][k];
			}
			if(p==0)
				printf("fds");
			if(isnan(p))
				printf("nan ");
			perp += log(p);
		}
		numerator = perp;
		denominator = Nd;
		perplexity += exp(-numerator/denominator);
	}

	//perplexity = exp(-numerator/denominator);
	return perplexity;
}

//strategy 1: 只考虑trndata中的model的data
double model::compute_energy(pagerank page){
	double energy = 0.0;
	//int nodenumber = page.network.countnode();
	int nodenumber = this->M;
	for(int i = 0;i<nodenumber;i++){
		//int u = i;
		string uname = ptrndata->docs[i]->name;
		int u = page.network.fromnametoid(uname);
		linkset out = page.network.outLink(u);
		int udegree = out.size();
		//int udegree = 0;
		linkset::iterator it;
		double temp_energy = 0;
		for(it = out.begin();it!=out.end();++it){
			int v = *it;
			
			string nodename = page.network.fromidtoname(v);
			map<string, int>::iterator nit;
			nit = this->ptrndata->name2id.find(nodename);
			double diff = 0;
			if(nit!=this->ptrndata->name2id.end()){
				int vdocid = nit->second; 
				for(int j = 0;j<K;j++){
					double tmp = theta[i][j]/sqrt(page.getPagerankScore(u))-theta[vdocid][j]/sqrt(page.getPagerankScore(v));
					diff += pow(tmp, 2);
				}
				//udegree++;
			}
			temp_energy +=page.getPagerankScore(u)*diff;
			
		}
		if(udegree!=0){
			temp_energy = temp_energy/(double)udegree; 
		}
		
		energy+=temp_energy;
	}
	energy = energy/2; 
	return energy;
}

double model::compute_likelihood(){
	double doctopic = 0.0;
	double topicword = 0.0;
	
	//这样算不行，因为调用的gamma函数范围有限，不能算超过171.6的数。所以应采用求ln(|Gamma(x)|)的方式
	//for(int m = 0;m<M;m++){
	//	double numerator_up = 1.0;
	//	double numerator_down = 0.0;
	//	double denominator_up = 1.0;
	//	double denominator_down = 0.0;
	//	for(int k = 0;k<K;k++){
	//		numerator_up *= mathlib::lagam(nd[m][k]+alphaarray[m][k]);
	//		numerator_down +=nd[m][k]+alphaarray[m][k];
	//		denominator_up *= mathlib::lagam(alphaarray[m][k]);
	//		denominator_down += alphaarray[m][k];
	//	}
	//	numerator_down = mathlib::lagam(numerator_down);
	//	denominator_down = mathlib::lagam(denominator_down);
	//	
	//	double tmp = (numerator_up/numerator_down)/(denominator_up/denominator_down);
	//	doctopic +=log(tmp);
	//}

	//for(int k = 0;k<K;k++){
	//	double numerator_up = 1.0;
	//	double numerator_down = 0.0;
	//	double denominator_up = 1.0;
	//	double denominator_down = 0.0;
	//	for(int v = 0;v<V;v++){
	//		numerator_up *= mathlib::lagam(nw[v][k]+beta);
	//		numerator_down +=nw[v][k]+beta;
	//		denominator_up *= mathlib::lagam(beta);
	//		denominator_down += beta;
	//	}
	//	numerator_down = mathlib::lagam(numerator_down);
	//	denominator_down = mathlib::lagam(denominator_down);
	//	
	//	double tmp = (numerator_up/numerator_down)/(denominator_up/denominator_down);
	//	topicword +=log(tmp);
	//	
	//}

	double sign;
	double lngammabeta = lngamma(beta,sign);
	double Vlngammabeta = lngamma(V*beta, sign);
	
	for(int m = 0;m<M;m++){
		double alphasum = 0.0;
		double tmp = 0.0;
		for(int k = 0;k<K;k++){
			alphasum +=alphaarray[m][k];
			tmp +=lngamma(nd[m][k]+alphaarray[m][k], sign)-lngamma(alphaarray[m][k], sign);			
		}
		doctopic +=tmp + lngamma(alphasum,sign) - lngamma(ndsum[m]+alphasum, sign);
	}

	for(int k = 0;k<K;k++){
		double tmp = 0.0;
		for(int v = 0;v<V;v++){
			tmp +=lngamma(nw[v][k]+beta, sign)-lngammabeta;
		}
		topicword +=tmp + Vlngammabeta - lngamma(nwsum[k]+V*beta, sign);
	}
	return (doctopic+topicword);
}



double model::compute_objective(double factor, pagerank page){
	double result = 0.0;
	double energy = compute_energy(page);
	double likelihood = compute_likelihood();
	result = (1-factor)* energy - factor * likelihood;
	return result;
}

void model::smooth_distribution(double factor, pagerank page){
	double ** theta_tmp = new double*[M];
	for(int i = 0;i<M;i++){
		//update theta[m][0-k]
		theta_tmp[i] = new double[K];
		string nodename = ptrndata->docs[i]->name;
		int graphid = page.network.fromnametoid(nodename);
		linkset out = page.network.outLink(graphid);
		linkset in = page.network.inLink(graphid);
		linkset::iterator it; 
		map<string, int>::iterator nit;
		int vdegree  = out.size();

		for(int j = 0; j<K;j++){
			theta_tmp[i][j] = 0.0;
			double change = 0.0;
			
			for(it=in.begin();it!=in.end();++it){
				int u = *it;
				string uname = page.network.fromidtoname(u);
				nit = ptrndata->name2id.find(uname);
				if(nit!=ptrndata->name2id.end()){
					int udocid = nit->second;
					int udegree = page.network.outLink(u).size(); 
					change += (double)page.getPagerankScore(u)*theta[udocid][j]/(double)udegree/sqrt(page.getPagerankScore(u)*page.getPagerankScore(graphid));
				}
				
			}

			for(it=out.begin();it!=out.end();++it){
				int u = *it;
				string uname = page.network.fromidtoname(u);
				nit = ptrndata->name2id.find(uname);
				if(nit!=ptrndata->name2id.end()){
					int udocid = nit->second;
					change +=(double)page.getPagerankScore(graphid)*theta[udocid][j]/(double)vdegree/sqrt(page.getPagerankScore(u)*page.getPagerankScore(graphid));
				}
				
			}

			theta_tmp[i][j] = (1-factor)*theta[i][j]+factor*0.5*change;
		}

	}

	for(int i = 0;i<M;i++){
		for(int j = 0;j<K;j++){
			theta[i][j] = theta_tmp[i][j];
		}
	}
	delete theta_tmp;

}

void model::smooth_alpha(double factor, pagerank page){
	double ** alpha_tmp = new double*[M];
	for(int i = 0;i<M;i++){
		//update theta[m][0-k]
		alpha_tmp[i] = new double[K];
		string nodename = ptrndata->docs[i]->name;

		int graphid = page.network.fromnametoid(nodename);
		linkset out = page.network.outLink(graphid);
		linkset in = page.network.inLink(graphid);
		linkset::iterator it; 
		map<string, int>::iterator nit;
		int vdegree  = out.size();

		
		for(int j = 0; j<K;j++){
			alpha_tmp[i][j] = 0.0;
			double change = 0.0;
			
			
			for(it=in.begin();it!=in.end();++it){
				int u = *it;
				string uname = page.network.fromidtoname(u);
				nit = ptrndata->name2id.find(uname);
				if(nit!=ptrndata->name2id.end()){
					int udocid = nit->second;
					int udegree = page.network.outLink(u).size(); 
					change += (double)page.getPagerankScore(u)*theta[udocid][j]/(double)udegree/sqrt(page.getPagerankScore(u)*page.getPagerankScore(graphid));
				}
				
			}

			for(it=out.begin();it!=out.end();++it){
				int u = *it;
				string uname = page.network.fromidtoname(u);
				nit = ptrndata->name2id.find(uname);
				if(nit!=ptrndata->name2id.end()){
					int udocid = nit->second;
					change +=(double)page.getPagerankScore(graphid)*theta[udocid][j]/(double)vdegree/sqrt(page.getPagerankScore(u)*page.getPagerankScore(graphid));
				}
				
			}

			alpha_tmp[i][j] = (1-factor)*alphaarray[i][j]+factor*0.5*change;
		}

	}

	for(int i = 0;i<M;i++){
		for(int j = 0;j<K;j++){
			alphaarray[i][j] = alpha_tmp[i][j];
		}
	}
	delete alpha_tmp;

}

