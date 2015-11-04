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

#include <stdio.h>
#include <stdlib.h>
#include "constants.h"
#include "strtokenizer.h"
#include "utils.h"
#include "dataset.h"

using namespace std;

int dataset::write_wordmap(string wordmapfile, mapword2id * pword2id) {
    FILE * fout = fopen(wordmapfile.c_str(), "w");
    if (!fout) {
	printf("Cannot open file %s to write!\n", wordmapfile.c_str());
	return 1;
    }    
    
    mapword2id::iterator it;
    fprintf(fout, "%d\n", pword2id->size());
    for (it = pword2id->begin(); it != pword2id->end(); it++) {
	fprintf(fout, "%s %d\n", (it->first).c_str(), it->second);
    }
    
    fclose(fout);
    
    return 0;
}

int dataset::read_wordmap(string wordmapfile, mapword2id * pword2id) {
    pword2id->clear();
    
    FILE * fin = fopen(wordmapfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", wordmapfile.c_str());
	return 1;
    }    
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int nwords = atoi(buff);
    
    for (int i = 0; i < nwords; i++) {
	fgets(buff, BUFF_SIZE_SHORT - 1, fin);
	line = buff;
	
	strtokenizer strtok(line, " \t\r\n");
	if (strtok.count_tokens() != 2) {
	    continue;
	}
	
	pword2id->insert(pair<string, int>(strtok.token(0), atoi(strtok.token(1).c_str())));
    }
    
    fclose(fin);
    
    return 0;
}

int dataset::read_wordmap(string wordmapfile, mapid2word * pid2word) {
    pid2word->clear();
    
    FILE * fin = fopen(wordmapfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", wordmapfile.c_str());
	return 1;
    }    
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int nwords = atoi(buff);
    
    for (int i = 0; i < nwords; i++) {
	fgets(buff, BUFF_SIZE_SHORT - 1, fin);
	line = buff;
	
	strtokenizer strtok(line, " \t\r\n");
	if (strtok.count_tokens() != 2) {
	    continue;
	}
	
	pid2word->insert(pair<int, string>(atoi(strtok.token(1).c_str()), strtok.token(0)));
    }
    
    fclose(fin);
    
    return 0;
}

int dataset::read_trndata(string dfile, string wordmapfile, string docmapfile, string citationmapfile) {
    mapword2id word2id;
	 mapcitation2id citation2id; //add by liu liu

    FILE * fin = fopen(dfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", dfile.c_str());
	return 1;
    }   
    
    mapword2id::iterator it; 

	map<string, int>::iterator nit;
	map<int,string>::iterator iit;

	mapcitation2id::iterator cit;//add by liu liu

	char buff[BUFF_SIZE_LONG];
    string line, citations; //add by liu liu
	string docname;
	string token;
    
    // get the number of documents
    fgets(buff, BUFF_SIZE_LONG - 1, fin);
    M = atoi(buff);
    if (M <= 0) {
	printf("No document available!\n");
	return 1;
    }
    
    // allocate memory for corpus
    if (docs) {
	deallocate();
    } else {
	docs = new document*[M];
    }
    
    // set number of words to zero
    V = 0;
    
    for (int i = 0; i < M; i++) { //modified by liu liu
	fgets(buff, BUFF_SIZE_LONG - 1, fin);
	line = buff;
	
	//先分成name,words和citations两个部分
	strtokenizer word_citation(line, ";");
	string docname = utils::trimstring(word_citation.token(0));
	string words = utils::trimstring(word_citation.token(1));
	string citations = word_citation.token(2);
	
	strtokenizer strtok(words, " \t\r\n");
	int length = strtok.count_tokens();

	if (length <= 0) {
	    printf("Invalid (empty) document!\n");
	    deallocate();
	    M = V = 0;
	    return 1;
	}
	
	// allocate new document
	document * pdoc = new document(length);
	
	for (int j = 0; j < length; j++) {
	    it = word2id.find(strtok.token(j));
	    if (it == word2id.end()) {
		// word not found, i.e., new word
		pdoc->words[j] = word2id.size();
		word2id.insert(pair<string, int>(strtok.token(j), word2id.size()));
	    } else {
			pdoc->words[j] = it->second;
	    }
	}
	
	/*if(i==1300){
		printf("haha\n");
		printf("lala %d",atoi(citations.c_str()));
	}*/
	if(atoi(citations.c_str())==0){
		pdoc->setcitationcount(0);
	}
	else{
		citations = utils::trimstring(word_citation.token(2));
		strtokenizer strtokc(citations, ",");
		int citationcount = strtokc.count_tokens();
		pdoc->setcitationcount(citationcount);

		for (int j = 0;j<citationcount;j++){
			token  = utils::trimstring(strtokc.token(j));
			//pdoc->citations[j] = atoi(token.c_str());
			
			pdoc->citations[j] = token.c_str();
			cit = citation2id.find(atoi(token.c_str()));
			if(cit==citation2id.end()){
			citation2id.insert(pair<int, int>( atoi(token.c_str()), citation2id.size()));
			}
		}
	}
	

	pdoc->name = docname;
	
	// add new doc to the corpus
	add_doc(pdoc, i);
	nit = name2id.find(docname);
	if(nit==name2id.end()){
		name2id.insert(pair<string, int>(docname, name2id.size()));
		id2name.insert(pair<int, string>(name2id.size(), docname));
	}
   
	}
    fclose(fin);
    
    // write word map to file
    if (write_wordmap(wordmapfile, &word2id)) {
	return 1;
    }
    
	if(write_docmap(docmapfile, &name2id)){
	return 1;
	}

	if(write_citationmap(citationmapfile, &citation2id)){
 		return 1;
 	}

    // update number of words
    V = word2id.size();
    
    return 0;
}

int dataset::read_newdata(string dfile, string wordmapfile, string citationmapfile) {
    mapword2id word2id;
    map<int, int> id2_id;
	  mapcitation2id citation2id; //add by liu liu

    
    read_wordmap(wordmapfile, &word2id);
    if (word2id.size() <= 0) {
	printf("No word map available!\n");
	return 1;
    }

	read_citation2idmap(citationmapfile, &citation2id);
	if(citation2id.size()<=0){
		printf("No citations map available!\n");
		return 1;
	}


    FILE * fin = fopen(dfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", dfile.c_str());
	return 1;
    }   

    mapword2id::iterator it;
    map<int, int>::iterator _it;
	mapcitation2id::iterator rit; //add by liu liu

	map<string, int>::iterator nit;
	map<int,string>::iterator iit;


    char buff[BUFF_SIZE_LONG];
    string line, token;
    
    // get number of new documents
    fgets(buff, BUFF_SIZE_LONG - 1, fin);
    M = atoi(buff);
    if (M <= 0) {
	printf("No document available!\n");
	return 1;
    }
    
    // allocate memory for corpus
    if (docs) {
	deallocate();
    } else {
	docs = new document*[M];
    }
    _docs = new document*[M];
    
    // set number of words to zero
    V = 0;
    
	int i = 0;
	while(i<M) {
	fgets(buff, BUFF_SIZE_LONG - 1, fin);
	line = buff;

//先分成name,words和citations两个部分
	strtokenizer word_citation(line, ";");
	string docname = utils::trimstring(word_citation.token(0));
	string words = utils::trimstring(word_citation.token(1));
	string citations = utils::trimstring(word_citation.token(2));

	//words
	strtokenizer strtok(words, " \t\r\n");
	int length = strtok.count_tokens();
	
	vector<int> doc;
	vector<int> _doc;
	for (int j = 0; j < length; j++) {
	    it = word2id.find(strtok.token(j));
	    if (it == word2id.end()) {
		// word not found, i.e., word unseen in training data
		// do anything? (future decision)
	    } else {
		int _id;
		_it = id2_id.find(it->second);
		if (_it == id2_id.end()) {
		    _id = id2_id.size();
		    id2_id.insert(pair<int, int>(it->second, _id));
		    _id2id.insert(pair<int, int>(_id, it->second));
		} else {
		    _id = _it->second;
		}
		
		doc.push_back(it->second);
		_doc.push_back(_id);
	    }
	}
	
	// allocate memory for new doc
	document * pdoc = new document(doc);
	document * _pdoc = new document(_doc);
	
	strtokenizer strtokc(citations, ",");
	int citationcount = strtokc.count_tokens();
	pdoc->setcitationcount(citationcount);
	_pdoc->setcitationcount(citationcount);
	int ll = 0;
	for (int j = 0;j<citationcount;j++){
		token  = utils::trimstring(strtokc.token(j));
		rit = citation2id.find(atoi(token.c_str()));
		if(rit==citation2id.end()){
			printf("reference not exist\n");
			 pdoc->citationcount--;
			 _pdoc->citationcount--;
			 continue;
	 	} else {
		pdoc->citations[ll] = token.c_str();
		_pdoc->citations[ll] = token.c_str();
		ll++;
		}
	}

	pdoc->name = docname;
	_pdoc->name = docname;


	// add new doc
	if(pdoc->length>0&&pdoc->citationcount>0){
	add_doc(pdoc, i);
	_add_doc(_pdoc, i);
	i++;
	}else{
		M--;
	}
	

	nit = name2id.find(docname);
	if(nit==name2id.end()){
		name2id.insert(pair<string, int>(docname, name2id.size()));
		id2name.insert(pair<int, string>(name2id.size(), docname));
	}

    }
    
    fclose(fin);
    
    // update number of new words
    V = id2_id.size();
    
    return 0;
}

int dataset::read_newdata_withrawstrs(string dfile, string wordmapfile) {
    mapword2id word2id;
    map<int, int> id2_id;
    
    read_wordmap(wordmapfile, &word2id);
    if (word2id.size() <= 0) {
	printf("No word map available!\n");
	return 1;
    }

    FILE * fin = fopen(dfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", dfile.c_str());
	return 1;
    }   

    mapword2id::iterator it;
    map<int, int>::iterator _it;
    char buff[BUFF_SIZE_LONG];
    string line;
    
    // get number of new documents
    fgets(buff, BUFF_SIZE_LONG - 1, fin);
    M = atoi(buff);
    if (M <= 0) {
	printf("No document available!\n");
	return 1;
    }
    
    // allocate memory for corpus
    if (docs) {
	deallocate();
    } else {
	docs = new document*[M];
    }
    _docs = new document*[M];
    
    // set number of words to zero
    V = 0;
    
    for (int i = 0; i < M; i++) {
	fgets(buff, BUFF_SIZE_LONG - 1, fin);
	line = buff;
	strtokenizer strtok(line, " \t\r\n");
	int length = strtok.count_tokens();
	
	vector<int> doc;
	vector<int> _doc;
	for (int j = 0; j < length - 1; j++) {
	    it = word2id.find(strtok.token(j));
	    if (it == word2id.end()) {
		// word not found, i.e., word unseen in training data
		// do anything? (future decision)
	    } else {
		int _id;
		_it = id2_id.find(it->second);
		if (_it == id2_id.end()) {
		    _id = id2_id.size();
		    id2_id.insert(pair<int, int>(it->second, _id));
		    _id2id.insert(pair<int, int>(_id, it->second));
		} else {
		    _id = _it->second;
		}
		
		doc.push_back(it->second);
		_doc.push_back(_id);
	    }
	}
	
	// allocate memory for new doc
	document * pdoc = new document(doc, line);
	document * _pdoc = new document(_doc, line);
	
	// add new doc
	add_doc(pdoc, i);
	_add_doc(_pdoc, i);
    }
    
    fclose(fin);
    
    // update number of new words
    V = id2_id.size();
    
    return 0;
}

int dataset::write_docmap(string docmapfile, map<string, int> * docname2id) {
    FILE * fout = fopen(docmapfile.c_str(), "w");
    if (!fout) {
		printf("Cannot open file %s to write!\n", docmapfile.c_str());
		return 1;
    }    
    
    map<string, int>::iterator it;
	fprintf(fout, "%d\n", docname2id->size());
    for (it = docname2id->begin(); it != docname2id->end(); it++) {
		fprintf(fout, "%s\t%d\n", (it->first).c_str(), it->second);
    }
    
    fclose(fout);
    
    return 0;
}

//add by liu liu
int dataset::read_id2docmap(string docmapfile, map<int, string> * id2docname){
	id2docname->clear();
    
    FILE * fin = fopen(docmapfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", docmapfile.c_str());
	return 1;
    }    
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int ndocs = atoi(buff);
    
    for (int i = 0; i < ndocs; i++) {
	fgets(buff, BUFF_SIZE_SHORT - 1, fin);
	line = buff;
	
	strtokenizer strtok(line, "\t\r\n");
	if (strtok.count_tokens() != 2) {
	    continue;
	}

	string id = utils::trimstring(strtok.token(1));
	string docname = utils::trimstring(strtok.token(0));
	id2docname->insert(pair<int, string>(atoi(id.c_str()), docname));
    }
    
    fclose(fin);
    
    return 0;

}

//add by liu liu
int dataset::read_doc2idmap(string docmapfile, map<string, int> * docname2id){
	docname2id->clear();
    
    FILE * fin = fopen(docmapfile.c_str(), "r");
    if (!fin) {
		printf("Cannot open file %s to read!\n", docmapfile.c_str());
		return 1;
    }    
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int ndocs = atoi(buff);
    
    for (int i = 0; i < ndocs; i++) {
	fgets(buff, BUFF_SIZE_SHORT - 1, fin);
	line = buff;
	
	strtokenizer strtok(line, "\t\r\n");
	if (strtok.count_tokens() != 2) {
	    continue;
	}
	
	string id = utils::trimstring(strtok.token(1));
	string docname = utils::trimstring(strtok.token(0));	
	docname2id->insert(pair<string, int>(docname, atoi(id.c_str())));
    }
    
    fclose(fin);
    
    return 0;
}



//add by liu liu
int dataset::read_id2citationmap(string citationmapfile, mapid2citation * pid2citation){
	pid2citation->clear();
    
    FILE * fin = fopen(citationmapfile.c_str(), "r");
    if (!fin) {
	printf("Cannot open file %s to read!\n", citationmapfile.c_str());
	return 1;
    }    
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int ncitations = atoi(buff);
    
    for (int i = 0; i < ncitations; i++) {
	fgets(buff, BUFF_SIZE_SHORT - 1, fin);
	line = buff;
	
	strtokenizer strtok(line, "\t\r\n");
	if (strtok.count_tokens() != 2) {
	    continue;
	}

	string id = utils::trimstring(strtok.token(1));
	string citation = utils::trimstring(strtok.token(0));
	pid2citation->insert(pair<int, int>(atoi(id.c_str()), atoi(citation.c_str())));
    }
    
    fclose(fin);
    
    return 0;

}

//add by liu liu
int dataset::read_citation2idmap(string citationmapfile, mapcitation2id * pcitation2id){
	pcitation2id->clear();
    
    FILE * fin = fopen(citationmapfile.c_str(), "r");
    if (!fin) {
		printf("Cannot open file %s to read!\n", citationmapfile.c_str());
		return 1;
    }    
    
    char buff[BUFF_SIZE_SHORT];
    string line;
    
    fgets(buff, BUFF_SIZE_SHORT - 1, fin);
    int ncitations = atoi(buff);
    
    for (int i = 0; i < ncitations; i++) {
	fgets(buff, BUFF_SIZE_SHORT - 1, fin);
	line = buff;
	
	strtokenizer strtok(line, "\t\r\n");
	if (strtok.count_tokens() != 2) {
	    continue;
	}
	
	string id = utils::trimstring(strtok.token(1));
	string citation = utils::trimstring(strtok.token(0));	
	pcitation2id->insert(pair<int, int>(atoi(citation.c_str()), atoi(id.c_str())));
    }
    
    fclose(fin);
    
    return 0;
}



int dataset::write_citationmap(string citationmapfile, mapcitation2id * pcitation2id) {
    FILE * fout = fopen(citationmapfile.c_str(), "w");
    if (!fout) {
		printf("Cannot open file %s to write!\n", citationmapfile.c_str());
		return 1;
    }    
    
    mapcitation2id::iterator it;
    fprintf(fout, "%d\n", pcitation2id->size());
    for (it = pcitation2id->begin(); it != pcitation2id->end(); it++) {
		fprintf(fout, "%d\t%d\n", it->first, it->second);
    }
    
    fclose(fout);
    
    return 0;
}