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

#include "model.h"
#include "utils.h"


using namespace std;

void show_help();


int main(int argc, char ** argv) {
//LDA
	vector<string> vec;
	if(utils::isbatch(argc, argv, vec)){
		int ntopics[9] = {2, 5, 10, 20, 40, 60, 80, 100, 120};		

		for(int i=0; i<9; i++){
			char a[10];
			itoa(ntopics[i], a, 10);
			{
				int nparam = 9;
				char **params = new char*[nparam];
				params[0] = "-est";
				params[1] = "-ntopics";
				params[2] = a;
				params[3] = "-savestep";
				params[4] = "1000";
				params[5] = "-twords";
				params[6] = "10";
				params[7] = "-dfile";
				params[8] = (char*)vec[0].c_str();//"lda-train.txt";

				model lda;
				if (lda.init(nparam, params)) {
					show_help();
					return 1;
				}
			    
				if (lda.model_status == MODEL_STATUS_EST || lda.model_status == MODEL_STATUS_ESTC) {
					// parameter estimation
					lda.estimate();
				}
				delete params;
			}
			{
				int nparam = 9;
				char **params = new char*[nparam];
				params[0] = "-inf";
				params[1] = "-dir";
				params[2] = "./";
				params[3] = "-model";
				params[4] = "model-final";
				params[5] = "-twords";
				params[6] = "10";
				params[7] = "-dfile";
				params[8] = (char*)vec[1].c_str();//"lda-test.txt";

				model lda;
				if (lda.init(nparam, params)) {
					show_help();
					return 1;
				}

				if (lda.model_status == MODEL_STATUS_INF) {
					// do inference
					lda.inference();
				}
				delete params;
			}
		}
	}else{	
		model lda;

		if (lda.init(argc, argv)) {
			show_help();
			return 1;
		}
	    
		if (lda.model_status == MODEL_STATUS_EST || lda.model_status == MODEL_STATUS_ESTC) {
			// parameter estimation
			//-est -ntopics 2 -savestep 1000 -twords 20 -dfile nips-lda.txt
			lda.estimate();
		}
	    
		if (lda.model_status == MODEL_STATUS_INF) {
			// do inference

			lda.inference();
		}

		if(lda.model_status == MODEL_STATUS_EVALUATE) {
			lda.evaluate();
		}
	}    
    return 0;
	
}

void show_help() {
    printf("Command line usage:\n");
    printf("\tlda -est -alpha <double> -beta <double> -ntopics <int> -niters <int> -savestep <int> -twords <int> -dfile <string>\n");
    printf("\tlda -estc -dir <string> -model <string> -niters <int> -savestep <int> -twords <int>\n");
    printf("\tlda -inf -dir <string> -model <string> -niters <int> -twords <int> -dfile <string>\n");
    // printf("\tlda -inf -dir <string> -model <string> -niters <int> -twords <int> -dfile <string> -withrawdata\n");
}

