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

#ifndef	_DATASET_H
#define	_DATASET_H

#include <string>
#include <vector>
#include <map>

using namespace std;

// map of words/terms [string => int]
typedef map<string, int> mapword2id;
// map of words/terms [int => string]
typedef map<int, string> mapid2word;

//map of citation id to citation [int=>string]
typedef map<int, int> mapid2citation;
//map of citation to citation id [string => int]
typedef map<int, int> mapcitation2id;

class document {
public:
    int * words;
    string rawstr;
    int length;
    
	//add by liuliu
	string * citations;
	int citationcount;

	//add by liuliu
	string name;

    document() {
	words = NULL;
	rawstr = "";
	length = 0;	
	//add by liuliu
	citations = NULL;
	citationcount = 0;
	name = "";
    }
    
    document(int length) {
	this->length = length;
	rawstr = "";
	words = new int[length];	
	//add by liuliu
	citations = NULL;
	citationcount = 0;
    name = "";
    }
    
    document(int length, int * words) {
	this->length = length;
	rawstr = "";
	this->words = new int[length];
	for (int i = 0; i < length; i++) {
	    this->words[i] = words[i];
	}
    
	//add by liuliu
	citations = NULL;
	citationcount = 0;
	name = "";
    }

    document(int length, int * words, string rawstr) {
	this->length = length;
	this->rawstr = rawstr;
	this->words = new int[length];
	for (int i = 0; i < length; i++) {
	    this->words[i] = words[i];
	}
    //add by liuliu
	citations = NULL;
	citationcount = 0;
	name = "";
	}
    
    document(vector<int> & doc) {
	this->length = doc.size();
	rawstr = "";
	this->words = new int[length];
	
	for (int i = 0; i < length; i++) {
	    this->words[i] = doc[i];
	}
    //add by liuliu
	citations = NULL;
	citationcount = 0;
	name = "";
    }

    document(vector<int> & doc, string rawstr) {
	this->length = doc.size();
	this->rawstr = rawstr;
	this->words = new int[length];
	for (int i = 0; i < length; i++) {
	    this->words[i] = doc[i];
	}
    //add by liuliu
	citations = NULL;
	citationcount = 0;

	name = "";
    }
  

    ~document() {
	if (words) {
	    delete words;
	}
	if(citations){
		delete citations;
	}
    }

	//add by liuliu
	void setcitationcount(int count){
		if(count == 0)
			this->citations = NULL;
		else
			this->citations = new string[count];
		this->citationcount = count;
	}

	void setdocname(string name){
		this->name = name;
	}
};

class dataset {
public:
    document ** docs;
    document ** _docs; // used only for inference
    map<int, int> _id2id; // also used only for inference

	map<int, string> id2name;
	map<string, int> name2id;

    int M; // number of documents
    int V; // number of words
    
    dataset() {
	docs = NULL;
	_docs = NULL;
	M = 0;
	V = 0;
    }
    
    dataset(int M) {
	this->M = M;
	this->V = 0;
	docs = new document*[M];	
	_docs = NULL;
    }   
    
    ~dataset() {
	if (docs) {
	    for (int i = 0; i < M; i++) {
		delete docs[i];
	    }
		delete docs;
	}
	
	
	if (_docs) {
	    for (int i = 0; i < M; i++) {
		delete _docs[i];		
	    }
		delete _docs;
	}
	
    }
    
    void deallocate() {
	if (docs) {
	    for (int i = 0; i < M; i++) {
		delete docs[i];
	    }
	}
	delete docs;
	docs = NULL;

	if (_docs) {
	    for (int i = 0; i < M; i++) {
		delete _docs[i];
	    }
	}
	delete _docs;
	_docs = NULL;
    }
    
    void add_doc(document * doc, int idx) {
	if (0 <= idx && idx < M) {
	    docs[idx] = doc;
	}
    }   
    
    void _add_doc(document * doc, int idx) {
	if (0 <= idx && idx < M) {
	    _docs[idx] = doc;
	}
    }       

    static int write_wordmap(string wordmapfile, mapword2id * pword2id);
    static int read_wordmap(string wordmapfile, mapword2id * pword2id);
    static int read_wordmap(string wordmapfile, mapid2word * pid2word);
    
	static int write_docmap(string docmapfile, map<string, int> * docname2id);
	static int read_id2docmap(string docmapfile, map<int, string> * id2docname);
	static int read_doc2idmap(string docmapfile, map<string, int> * docname2id);

		static int write_citationmap(string citationmapfile, mapcitation2id * pcitation2id);
	static int read_id2citationmap(string citationmapfile, mapid2citation * pid2citation);
	static int read_citation2idmap(string citationmapfile, mapcitation2id * pcitation2id);

   int read_trndata(string dfile, string wordmapfile, string docmapfile, string citationmapfile);
    int read_newdata(string dfile, string wordmapfile,string citationmapfile);
    int read_newdata_withrawstrs(string dfile, string wordmapfile);
};

#endif

