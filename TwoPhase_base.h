//
//  TwoPhase_base.h
//  TwoPhase
//
//  Created by Fei Gao on 1/11/18.
//  Copyright © 2018 Fei Gao. All rights reserved.
//

#ifndef TwoPhase_base_h
#define TwoPhase_base_h
#include "Eigen/Dense"
using namespace Eigen;
using namespace std;

double K_epan (VectorXd zdiff){
    if (zdiff.norm()<1) return(0.75*(1-pow(zdiff.norm(),2)));
    else return 0;
}
// Sort VectorXd increasingly and with only unique values
VectorXd VecSortUniq(VectorXd x){
    long p=x.size(); double x1;
    for (int i=1;i<p;++i){
        for(int j=0;j<p-i;++j) {
            if (x(j)>=x(j+1)) {
                x1=x(j); x(j)=x(j+1); x(j+1)=x1;
            }
        }
    }
    VectorXd z(1); z(0)=x(0); int linez=0;
    for(int i=1;i<p;++i){
        if (x(i)>z(linez)){
            VectorXd z1=z; linez=linez+1; z.resize(linez+1);
            for (int i=0;i<linez;++i) z(i)=z1(i);
            z(linez)=x(i);
        }
    }
    return(z);
}

VectorXd get_distinct_time(VectorXd L, VectorXd R, VectorXi Rinf){
    //VectorXd Luniq = VecSortUniq(L);
    VectorXd Luniq = L;
    long n1 = Luniq.size(); long n2 = R.size();
    VectorXd x(n1+n2+1); x.head(n1)=Luniq;
    int line=(int) n1;
    for (int i=0;i<n2;++i) {
        if (Rinf(i)==0){//&(R(i)<Luniq(n1-1))) {
            x(line)=R(i); line++;
        }
    }
    //x(line) = 0; line++;//add zero to the jump set
    VectorXd x1=x.head(line); return(VecSortUniq(x1));
}

//Lambda -> lambda
VectorXd cumdiff(VectorXd x){ long n=x.size(); for (int i=(int)n-1;i>0;--i) x(i) -= x(i-1); return(x);}
//lambda -> Lambda
VectorXd cumsum(VectorXd x){ long n=x.size(); for (int i=1;i<n;++i) x(i) += x(i-1); return(x);}

//get line
VectorXi gettline(VectorXd t, VectorXd Y){
    long m=t.size(); long n = Y.size(); VectorXi Yline(n);
    for (int i=0;i<n;++i){
        int l=0;
        while (l<m){
            if (t(l)<=Y(i)) l++;
            else break;
        }
        Yline(i)=l-1;
    }
    return(Yline);
}
VectorXi gettline(VectorXd t, VectorXd Y, VectorXi start){
    long m=t.size(); long n = Y.size(); VectorXi Yline(n);
    for (int i=0;i<n;++i){
        int l=start(i);
        while (l<m){
            if (t(l)<=Y(i)) l++;
            else break;
        }
        Yline(i)=l-1;
    }
    return(Yline);
}

double Glog(double x, double r){ double Gx=x; if (r>0) Gx=log(1.0+r*x)/r; return(Gx);}
double Glogp(double x, double r){ double Gx=1; if (r>0) Gx=1/(1.0+r*x); return(Gx);}
double Gloginv(double x, double r){ double Ginvx=x; if (r>0) Ginvx=(exp(r*x)-1)/r; return(Ginvx);}

/* Compute the grid and weight for Hermite integration*/
// Matrix: ngrid*2: weight, grid
MatrixXd Hermite_grid (long ngrid){
    MatrixXd J(ngrid,ngrid); J.setZero();
    for (int j=0;j<ngrid-1;++j){
        J(j,j+1)=pow(((double)j+1.0)/2.0,0.5);
        J(j+1,j)=J(j,j+1);
    }
    EigenSolver<MatrixXd> es(J);
    VectorXcd x=es.eigenvalues();
    MatrixXcd y=es.eigenvectors();
    
    MatrixXd grid(ngrid,2);
    for (int j=0;j<ngrid;++j){
        grid(j,0)=pow(y(0,j).real(),2)*pow(M_PI,0.5);
        grid(j,1)=x(j).real();
    }
    return(grid);
}

/* Compute the grid and weight for Laguerre integration*/
// Matrix: ngrid*2: weight, grid
MatrixXd Laguerre_grid (long ngrid, double alpha){
    MatrixXd J(ngrid,ngrid); J.setZero();
    for (int j=0;j<ngrid-1;++j){
        J(j,j)= 2*(double)j+alpha+1.0;
        J(j,j+1)=pow(((double)j+alpha+1.0)*((double)j+1.0),0.5);
        J(j+1,j)=J(j,j+1);
    }
    J(ngrid-1,ngrid-1)=2*(double)ngrid+alpha-1.0;
    EigenSolver<MatrixXd> es(J);
    VectorXcd x=es.eigenvalues();
    MatrixXcd y=es.eigenvectors();
    
    MatrixXd grid(ngrid,2);
    for (int j=0;j<ngrid;++j){
        grid(j,0)=pow(y(0,j).real(),2)*tgamma(alpha+1.0);
        grid(j,1)=x(j).real();
    }
    return(grid);
}

/* generate a sequence of U with Uniform(c1,c2)*/
vector<double> genC(double A, double tau, double c1, double c2, default_random_engine & generator,int& nint){
    uniform_real_distribution<double> gapdist(c1,c2);
    double C1 = A; double C2 = C1 + 0.1 + gapdist(generator);
    vector<double> C;
    if (C2<tau){
        while(C2<tau){
            C.push_back(C2); nint++;
            C1 = C2; C2 = C1 + 0.1 + gapdist(generator);
        }
    }
    return(C);
}
/* generate failure time from Lambda = log(1+at) and Gr = log(1+rx)/r */
double genlog1aT(double a, double r, double expbetaX, default_random_engine & generator){
    uniform_real_distribution<double> udist(0.0,1.0);
    double u=udist(generator); double Lambda = Gloginv(-log(u),r);
    return((exp(Lambda/expbetaX)-1.0)/a);
}
/* generate failure time from Lambda = at */
double genaT(double a, double r, double expbetaX, default_random_engine & generator){
    uniform_real_distribution<double> udist(0.0,1.0);
    double u=udist(generator); double Lambda = Gloginv(-log(u),r);
    return(Lambda/(expbetaX * a));
}


#endif /* TwoPhase_base_h */