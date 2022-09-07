//
//  Trun_full_simu.cpp
//  Trun_int
//
//  Created by Fei Gao on 3/31/16.
//  Copyright © 2016 Fei Gao. All rights reserved.
//
//input --n or --nsub sample size
//      --nrep #replicates for this file
//      --out output name for file
//      --seed seed
// Optional:
//      --hn hn for profile-likelihood
//      --r for the model
//      --ngrid number of grid for Gaussian quadrature

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <random>
#include <string>
#include <stdlib.h>
#include "TwoPhase.h"
using namespace Eigen;
using namespace std;

string int2string(int x){
    char temp[64];
    string str;
    sprintf(temp, "%d", x);
    string xs(temp);
    return(xs);
}
bool ParseCommandLineArgs(const int argc, char* argv[],int* seed, int* nsub, int* nrep, double* hn, double* r, int* ngrid, string* out_file) {
    // Initialize with bad parameter, so we can check at end value was updated.
    *seed = -1;
    *nsub = -1;
    *nrep = -1;
    string inputstring;
    // Loop through command line arguments, starting at 1 (first arg is the command
    // to run the program).
    for (int i = 1; i < argc; ++i) {
        string arg = string(argv[i]);
        if (arg == "--seed") {
            if (i == argc - 1) {
                cout << "\nERROR Reading Command: Expected argument after '--seed'.\n"
                << "Aborting.\n";
                return false;
            }
            ++i; inputstring=string(argv[i]); *seed = atoi(inputstring.c_str());
        } else if (arg == "--nsub" || arg == "--n") {
            if (i == argc - 1) {
                cout << "\nERROR Reading Command: Expected argument after '--n'.\n"
                << "Aborting.\n";
                return false;
            }
            ++i; inputstring=string(argv[i]); *nsub = atoi(inputstring.c_str());
        } else if (arg == "--nrep") {
            if (i == argc - 1) {
                cout << "\nERROR Reading Command: Expected argument after '--nrep'.\n"
                << "Aborting.\n";
                return false;
            }
            ++i; inputstring=string(argv[i]); *nrep = atoi(inputstring.c_str());
        } else if (arg == "--hn") {
            if (i == argc - 1) {
                cout << "\nERROR Reading Command: Expected argument after '--hn'.\n"
                << "Aborting.\n";
                return false;
            }
            ++i; inputstring=string(argv[i]); *hn = atoi(inputstring.c_str());
        } else if (arg == "--r") {
            if (i == argc - 1) {
                cout << "\nERROR Reading Command: Expected argument after '--r'.\n"
                << "Aborting.\n";
                return false;
            }
            ++i; inputstring=string(argv[i]); *r = atof(inputstring.c_str());
        } else if (arg == "--ngrid") {
            if (i == argc - 1) {
                cout << "\nERROR Reading Command: Expected argument after '--ngrid'.\n"
                << "Aborting.\n";
                return false;
            }
            ++i; inputstring=string(argv[i]); *ngrid = atoi(inputstring.c_str());
        } else if (arg == "--out") {
            if (i == argc - 1) {
                cout << "\nERROR Reading Command: Expected argument after '--out'.\n";
                return false;
            }
            ++i; *out_file = string(argv[i]);
        } else {
            cout << "\nERROR Reading Command: unrecognized argument: " << arg << endl;
            return false;
        }
    }
    
    return *nsub != -1 && *nrep != -1 && *seed != -1;
}


int main(int argc, char* argv[]){
    // Parse simulation parameters from command line.
    int seed,nsub,nrep; double hn=0.0; double r=0; int ngrid = 20; string out_file;
    if (!ParseCommandLineArgs(argc, argv, &seed, &nsub, &nrep, &hn, &r, &ngrid, &out_file)) { return -1;}
    
    default_random_engine generator(seed);
    VectorXd beta(1); beta << -0.5; VectorXd gamma(2); gamma << 0.5, 0.8;
    
    string out_file_beta = out_file + "_beta.dat";
    string out_file_lambda = out_file + "_lambda.dat";
    ofstream myfile_beta (out_file_beta);
    ofstream myfile_lambda (out_file_lambda);
    int nt = 400; VectorXd t (nt); VectorXd Lambda(nt);
    for (int l=0;l<nt;++l) t(l) = (double) l / (double) nt * 3;
    if (myfile_beta.is_open()) {
        for (int k=0;k<nrep;++k){
            TwoPhase data=simudata_depZ(nsub, beta, gamma, generator, r, ngrid);
            data.solve();
            myfile_beta<<data.beta_.transpose()<<" "<<data.gamma_.transpose()<<" "<<data.iter_<<" ";
           
            VectorXi tline = gettline(data.t_,t);
            VectorXd dataLambda = cumsum(data.lambda_);
            for (int l=0;l<nt;++l) { if (tline(l)>=0) Lambda(l) = dataLambda(tline(l)); else Lambda(l) = 0;}
            myfile_lambda<<Lambda.transpose()<<endl;
            
            data.est_sd(10.0); myfile_beta<<data.sd_score.transpose()<<endl;
        }
        myfile_beta.close(); myfile_lambda.close();
        cout << "File written!";
    }
    else cout << "Unable to open file";
    return(0);
}
