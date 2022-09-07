//
//  main.cpp
//  Trun_int
//
//  Created by Fei Gao on 11/08/17.
//  Copyright © 2017 Fei Gao. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <vector>
#include <random>
#include <string>
#include <stdlib.h>
#include "TwoPhase_discZ_int.h"
using namespace Eigen;
using namespace std;

string int2string(int x){
    char temp[64];
    string str;
    sprintf(temp, "%d", x);
    string xs(temp);
    return(xs);
}
bool ParseCommandLineArgs(const int argc, char* argv[], double* hn, double* r, int* ngrid, string* out_file) {
    // Initialize with bad parameter, so we can check at end value was updated.
    
    string inputstring;
    // Loop through command line arguments, starting at 1 (first arg is the command
    // to run the program).
    for (int i = 1; i < argc; ++i) {
        string arg = string(argv[i]);
        if (arg == "--hn") {
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
    
    return true;
}

int main(int argc, char* argv[]){
    
    // Parse parameters from command line.
    double hn=0.0; double r=0; int ngrid = 20; string out_file;
    if (!ParseCommandLineArgs(argc, argv, &hn, &r, &ngrid, &out_file)) { return -1;}
    
    int pX = 1; int pZ = 3; int pW = 11;
    int p = 4 + pX + pZ + pW;
    
    string inputfile = "data.dat";
    string out_beta_file=out_file + "_beta.dat";
    string out_lambda_file=out_file + "_lambda.dat";
    string out_file_res = out_file + "_res.dat";
    ofstream myfile_beta (out_beta_file);
    ofstream myfile_lambda (out_lambda_file);
    ofstream myfile (out_file_res);
    if (myfile.is_open())
    {   ifstream data_file(inputfile);
        if (! data_file.is_open()){
            cout << "Unable to open input file" << endl;
            exit(1);
        }
        int nsub=1706;
        MatrixXd dat(nsub,p);
        // read in the data
        for (int i=0; i<nsub;++i) { for (int j=0; j<p && data_file >> dat(i,j); ++j) {}}
        VectorXd L=dat.col(0); VectorXd R=dat.col(1);
        MatrixXd X = dat.block(0,4,nsub,pX);
        MatrixXd Z = dat.block(0,4+pX,nsub,pZ);
        MatrixXd W = dat.block(0,4+pX+pZ,nsub,pW);
        VectorXi Rinf(nsub); VectorXi Q(nsub);
        for (int i=0;i<nsub;++i) { Rinf(i)=dat(i,2); Q(i)=dat(i,3);}
        TwoPhase data(L, R, Rinf, Q, X, Z, W, r, ngrid);
        
        data.solve();
        myfile_beta<<data.beta_.transpose()<<" "<<data.gamma_.transpose()<<endl;
        for (int j=0;j<data.m_;++j) myfile_lambda<<data.t_(j)<<" "<<data.lambda_(j)<<endl;
        myfile<<data.iter_<<" "<<data.logLik_MLE_i.sum()<<endl;
        
        data.est_sd(hn);
        myfile_beta<<data.Cov_score<<endl;
        myfile_beta.close(); myfile_lambda.close(); myfile.close();
        cout << "File written!";
    }
    else cout << "Unable to open file";
    return(0);
}
