//
//  TwoPhase.h
//  TwoPhase
//
//  Created by Fei Gao on 1/11/18.
//  Copyright © 2018 Fei Gao. All rights reserved.
//

#ifndef TwoPhase_h
#define TwoPhase_h

#include <vector>
#include "TwoPhase_base.h"

class TwoPhase{
private:
    VectorXd beta0_, gamma0_, lambda0_; MatrixXd puv0_; double lik0_;
    VectorXd Score_; MatrixXd Hessian_;
    VectorXd xiexpXW_; vector<VectorXd> xiXexpXW_; vector<MatrixXd> xiX2expXW_;
    MatrixXd Vil_, wiv_; vector<MatrixXd> VilX_;
    VectorXd logLik_i_;
    VectorXi Lline_,Rline_,Rstar_;
    MatrixXd pI_score, pI_infor;
    MatrixXd w; double an; // bandwidth
    MatrixXd x_, z_; long nx_,nz_;
    VectorXd Xline_, Zline_;
    MatrixXd Gxi_; long ngrid_; double r_;
    
    
    const double epsilon_lambda=1E-3; const double epsilon_beta=1E-4; const double epsilon_lik=1E-6;
    const double maxbound = 1E3; const int maxiter_=1E5;
    
    void class_setup(VectorXd L, VectorXd R, VectorXi Rinf, VectorXi P2, MatrixXd X, MatrixXd Z, MatrixXd W, double r, int ngrid){
        L_ = L; R_ = R; Rinf_ = Rinf; P2_ = P2; X_ = X; Z_ = Z; W_ = W;
        n_ = L_.size(); pX_ = X_.cols(); pZ_ = Z_.cols(); pW_ = W.cols(); p_ = pX_ + pW_+1;
        r_ = r;
        if (r_>0)  { Gxi_= Laguerre_grid(ngrid,1.0/r_-1.0); Gxi_.col(1) = Gxi_.col(1) * r_;}
        else { Gxi_.resize(1,2); Gxi_ << 1,1;}
        ngrid_ = Gxi_.rows();
        
        beta_.resize(pX_+1); beta_.setZero();
        gamma_.resize(pW_); gamma_.setZero();
        t_ = get_distinct_time(L_, R_, Rinf_); m_ = t_.size();
        lambda_.resize(m_); lambda_.setOnes(); lambda_ = lambda_/(double)(m_);
        getline(); uniqueXZ();
        
        //cout<<"nx="<<nx_<<", nz="<<nz_<<endl;
        puv_.resize(nx_,nz_); puv_.setOnes(); puv_ /= (double)nx_;
    
        logLik_i_.resize(n_);
        xiexpXW_.resize(n_); xiXexpXW_.resize(n_); xiX2expXW_.resize(n_);
        Vil_.resize(n_,m_); wiv_.resize(n_,nz_); VilX_.resize(m_);
        for (int i=0;i<n_;++i) { xiXexpXW_[i].resize(pX_+1); xiX2expXW_[i].resize(pX_+1,pX_+1);}
        for (int l=0;l<m_;++l) VilX_[l].resize(n_,pX_+1);
        cal_w();
        logLik_i_.setZero();
        Score_.resize(p_); Hessian_.resize(p_,p_);
    }
    
    //Get the lines in t and s for L, R, and U
    void getline(){
        Lline_ = gettline(t_,L_); Rline_ = gettline(t_,R_,Lline_); Rstar_ = Rline_;
        for (int i=0;i<n_;++i){
            if (R_(i) > t_(m_-1)) Rinf_(i)=1;
            if (Rinf_(i)==1) { Rline_(i) = (int) m_; Rstar_(i) = Lline_(i);}
        }
    }
    
    // Obtain unique values and corresponding lines for X and Z
    void uniqueXZ(){
        x_.resize(n_,pX_); z_.resize(n_,pZ_); Xline_.resize(n_); Zline_.resize(n_);
        int linex = 0; int linez = 0; bool match;
        for (int i=0;i<n_;++i){
            if (P2_(i)==1){
                if (linex ==0) { x_.row(linex) = X_.row(i); Xline_(i) = linex; linex++;}
                else{
                    match = false;
                    for (int line=0;line<linex;++line){
                        int j=0;
                        while(j<pX_){ if (X_(i,j)!=x_(line,j)) break; else ++j;}
                        if (j==pX_) { match = true; Xline_(i) = line; break;}
                    }
                    if (!match) { x_.row(linex) = X_.row(i); Xline_(i) = linex; linex++;}
                }
            }
            if (linez ==0) { z_.row(linez) = Z_.row(i); Zline_(i) = linez; linez++;}
            else{
                match = false;
                for (int line=0;line<linez;++line){
                    int j=0;
                    while(j<pZ_){ if (Z_(i,j)!=z_(line,j)) break; else ++j;}
                    if (j==pZ_) { match = true; Zline_(i) = line; break;}
                }
                if (!match) { z_.row(linez) = Z_.row(i); Zline_(i) = linez; linez++;}
            }
        }
        nx_ = (long)linex; nz_ = (long)linez;
        x_ = x_.block(0,0,linex,pX_); z_ = z_.block(0,0,linez,pZ_);
        
        cout<<"nx="<<nx_<<", nz= "<<nz_<<endl;
    }
    
    void cal_w(){
        wiv_.setZero();
        for (int i=0;i<n_;++i) wiv_(i,Zline_(i)) = 1;
    }
    
    MatrixXd calint(VectorXd beta, VectorXd gamma, VectorXd lambda, MatrixXd puv) {
        VectorXd Lambda = cumsum(lambda);
        logLik_i_.setZero(); xiexpXW_.setZero(); Vil_.setZero();
        double lik_xi, lik_xiu, lik_u, mus, xi, liki, denom, expXW;
        MatrixXd puv_new(nx_,nz_); puv_new.setZero();
        MatrixXd IXZ(nx_,nz_); VectorXd puvwiv(nz_);
        for (int l=0;l<m_;++l) VilX_[l].setZero();
        
        for (int i=0;i<n_;++i){
            liki=0;
            xiXexpXW_[i].setZero(); xiX2expXW_[i].setZero(); IXZ.setZero();
            if (P2_(i)==1){
                expXW = exp(W_.row(i)*gamma) * exp(X_.row(i)*beta.head(pX_)) * exp(X_(i,0)*W_(i,0)*beta(pX_));
                for (int s=0;s<ngrid_;++s){
                    mus = Gxi_(s,0); xi = Gxi_(s,1);
                    lik_xi = exp(-xi * Lambda(Lline_(i))* expXW);
                    if (Rinf_(i)==0) {
                        lik_xi -= exp(-xi * Lambda(Rline_(i))* expXW);
                        denom = 1-exp(-xi * expXW * (Lambda(Rline_(i))-Lambda(Lline_(i))));
                        for (int l=(Lline_(i)+1);l<=Rline_(i);++l){
                            Vil_(i,l) += mus * lik_xi * xi * lambda(l) * expXW / denom;
                            VilX_[l].block(i,0,1,pX_) += mus * lik_xi * xi * lambda(l) * expXW / denom * X_.row(i);
                            VilX_[l](i,pX_) += mus * lik_xi * xi * lambda(l) * expXW / denom * X_(i,0)*W_(i,0);
                        }
                    }
                    xiexpXW_(i) += mus * lik_xi * xi * expXW;
                    xiXexpXW_[i].head(pX_) += mus * lik_xi * xi * expXW * X_.row(i).transpose();
                    xiXexpXW_[i](pX_) += mus * lik_xi * xi * expXW * X_(i,0)*W_(i,0);
                    xiX2expXW_[i].block(0,0,pX_,pX_) += mus * lik_xi * xi * expXW * X_.row(i).transpose() * X_.row(i);
                    xiX2expXW_[i].block(pX_,0,1,pX_) += mus * lik_xi * xi * expXW * X_.row(i)* X_(i,0)*W_(i,0);
                    xiX2expXW_[i].block(0,pX_,pX_,1) = xiX2expXW_[i].block(pX_,0,1,pX_).transpose();
                    xiX2expXW_[i](pX_,pX_) += mus * lik_xi * xi * expXW * pow(X_(i,0)*W_(i,0),2);
                }
                liki = exp(-Glog(Lambda(Lline_(i))*expXW,r_));
                if (Rinf_(i)==0) liki -= exp(-Glog(Lambda(Rline_(i))*expXW,r_));
                
                for (int v=0;v<nz_;++v) {
                    IXZ(Xline_(i),v) = wiv_(i,v);
                    if (wiv_(i,v)>0) logLik_i_(i) += wiv_(i,v) * log(puv(Xline_(i),v));
                }
            } else {
                for (int u=0;u<nx_;++u){
                    expXW = exp(W_.row(i)*gamma) * exp(x_.row(u)*beta.head(pX_)) * exp(x_(u,0)*W_(i,0)*beta(pX_));
                    for (int v=0;v<nz_;++v) puvwiv(v) = puv_(u,v) * wiv_(i,v);
                    for (int s=0;s<ngrid_;++s){
                        mus = Gxi_(s,0); xi = Gxi_(s,1);
                        lik_xiu = exp(-xi * Lambda(Lline_(i)) * expXW);
                        if (Rinf_(i)==0) lik_xiu -= exp(-xi * Lambda(Rline_(i))*expXW);
                        lik_xiu *= puvwiv.sum();
                        
                        xiexpXW_(i) += mus * lik_xiu * xi * expXW;
                        xiXexpXW_[i].head(pX_) += mus * lik_xiu * xi * expXW * x_.row(u).transpose();
                        xiXexpXW_[i](pX_) += mus * lik_xiu * xi * expXW * x_(u,0)*W_(i,0);
                        xiX2expXW_[i].block(0,0,pX_,pX_) += mus * lik_xiu * xi * expXW * x_.row(u).transpose() * x_.row(u);
                        xiX2expXW_[i].block(pX_,0,1,pX_) += mus * lik_xiu * xi * expXW * x_.row(u)* x_(u,0)*W_(i,0);
                        xiX2expXW_[i].block(0,pX_,pX_,1) = xiX2expXW_[i].block(pX_,0,1,pX_).transpose();
                        xiX2expXW_[i](pX_,pX_) += mus * lik_xiu * xi * expXW * pow(x_(u,0)*W_(i,0),2);
                        
                        if (Rinf_(i)==0) {
                            denom = 1-exp(-xi * expXW * (Lambda(Rline_(i))-Lambda(Lline_(i))));
                            for (int l=(Lline_(i)+1);l<=Rline_(i);++l){
                                Vil_(i,l) += mus * lik_xiu * xi * lambda(l) * expXW /denom;
                                VilX_[l].block(i,0,1,pX_) += mus * lik_xiu * xi * lambda(l) * expXW /denom * x_.row(u);
                                VilX_[l](i,pX_) += mus * lik_xiu * xi * lambda(l) * expXW /denom * x_(u,0)*W_(i,0);
                            }
                        }
                    }
                    lik_u = exp(-Glog(Lambda(Lline_(i))*expXW,r_));
                    if (Rinf_(i)==0) lik_u -= exp(-Glog(Lambda(Rline_(i))*expXW,r_));
                    for (int v=0;v<nz_;++v) IXZ(u,v) = lik_u * puv(u,v) * wiv_(i,v); 
                    liki += IXZ.row(u).sum();
                }
                IXZ /= liki;
            }
            xiexpXW_(i) /= liki; xiXexpXW_[i] /= liki; xiX2expXW_[i] /= liki; Vil_.row(i) /= liki;
            for (int l=0;l<m_;++l) VilX_[l].row(i)/= liki;
            puv_new += IXZ; logLik_i_(i) += log(liki);
        }
        for (int v=0;v<nz_;++v) puv_new.col(v) /= puv_new.col(v).sum();
        return(puv_new);
    }
    
    VectorXd Mstep(VectorXd& beta, VectorXd& gamma){
        Score_.setZero(); Hessian_.setZero(); VectorXd lambda(m_);
        double sumtheta=0; VectorXd sumthetaZ(p_); MatrixXd sumthetaZZ(p_,p_);
        sumthetaZ.setZero(); sumthetaZZ.setZero();
        for (int l=(int)m_-1;l>=0;--l){
            for (int i=0;i<n_;++i){
                if (l==Rstar_(i)){
                    sumtheta += xiexpXW_(i);
                    sumthetaZ.head(pX_+1) += xiXexpXW_[i];
                    sumthetaZ.tail(pW_) += xiexpXW_(i) * W_.row(i).transpose();
                    sumthetaZZ.block(0,0,pX_+1,pX_+1) += xiX2expXW_[i];
                    sumthetaZZ.block(pX_+1,0,pW_,pX_+1) += W_.row(i).transpose() * xiXexpXW_[i].transpose();
                    sumthetaZZ.block(0,pX_+1,pX_+1,pW_) = sumthetaZZ.block(pX_+1,0,pW_,pX_+1).transpose();
                    sumthetaZZ.block(pX_+1,pX_+1,pW_,pW_) += xiexpXW_(i) * W_.row(i).transpose() * W_.row(i);
                }
            }
            lambda(l) = Vil_.col(l).sum()/sumtheta;
            for (int p=0;p<pX_+1;++p) Score_(p) += VilX_[l].col(p).sum();
            Score_.tail(pW_) += W_.transpose() * Vil_.col(l);
            Score_ -= lambda(l) * sumthetaZ;
            Hessian_ += lambda(l) * (sumthetaZZ-sumthetaZ*sumthetaZ.transpose()/sumtheta);
        }
        VectorXd step = Hessian_.inverse()*Score_;
        beta += step.head(pX_+1); gamma += step.tail(pW_);
        return(lambda);
    }
    
    VectorXd Msteplambda() {
        double denom = 0; VectorXd lambda(m_);
        for (int l=(int)m_-1;l>=0;--l){
            for (int i=0;i<n_;++i) { if (l==Rstar_(i)) denom += xiexpXW_(i);}
            lambda(l) = Vil_.col(l).sum()/denom;
        }
        return(lambda);
    }
    
    double diff(VectorXd v1,VectorXd v2){
        VectorXd diff=v1-v2;
        if (v2.norm()>=0.01) return(diff.norm());
        else return(diff.norm()/v2.norm());
    }
    
    double diff_mat(MatrixXd v1,MatrixXd v2){
        long p = v1.cols(); double x = 0;
        for (int j=0;j<p;++j) x = max(diff(v1.col(j),v2.col(j)),x);
        return(x);
    }
    
    VectorXd profilelogLik_i(VectorXd beta,VectorXd gamma){
        iter_=0;
        cout<<"beta="<<beta.transpose()<<", gamma="<<gamma.transpose()<<endl;
        VectorXd lambda=lambda_; VectorXd lambda0;
        MatrixXd puv = puv_; MatrixXd puv0;
        puv = calint(beta,gamma,lambda,puv);
        while (iter_<=maxiter_){
            lambda0=lambda; puv0 = puv; lik0_ = logLik_i_.mean();
            lambda=Msteplambda();
            puv = calint(beta, gamma, lambda, puv0);
            cout<<"lambda = "<<lambda_.head(5).transpose()<<" "<<lambda_.tail(5).transpose()<<endl;
            cout<<"iter="<<iter_<<", lambdadiff="<<diff(lambda,lambda0)<<",diff_lik = "<<abs(logLik_i_.mean()-lik0_)<<endl;
            cout<<" logLik="<<logLik_i_.mean()<<endl;
            iter_++;
            if ((diff(lambda,lambda0)<epsilon_lambda)|(abs(logLik_i_.mean()-lik0_)<epsilon_lik)) break;
        }
        return(logLik_i_);
    }
    
    void Profile_I(double h){
        VectorXd beta,betas,gamma,gammas; MatrixXd pL_i(n_, p_);
        /* calculate the value of pli(beta+h e)*/
        for (int j=0;j<pX_+1;++j){ beta=beta_; gamma = gamma_; beta(j)+=h; pL_i.col(j)=profilelogLik_i(beta,gamma);}
        for (int j=0;j<pW_;++j){ beta=beta_; gamma = gamma_; gamma(j)+=h; pL_i.col(pX_+1+j)=profilelogLik_i(beta,gamma);}
        pI_score.resize(p_,p_); pI_score.setZero(); MatrixXd pL_i_c=pL_i;
        for (int j=0;j<p_;++j) pL_i_c.col(j) = (pL_i_c.col(j)-logLik_MLE_i)/h;
        for (int i=0;i<n_;++i) pI_score += pL_i_c.row(i).transpose()*pL_i_c.row(i);
    }
    
public:
    VectorXd L_, R_; VectorXi Rinf_, P2_; MatrixXd X_, Z_, W_;
    long n_, pX_, pZ_, pW_, p_, m_;
    VectorXd beta_, gamma_, lambda_, t_;
    MatrixXd puv_;
    VectorXd logLik_MLE_i; int iter_; bool cov_;
    MatrixXd Cov_score, Cov_infor; VectorXd sd_score, sd_infor;
    
    TwoPhase() = default;
    TwoPhase(VectorXd L, VectorXd R, VectorXi Rinf, VectorXi P2, MatrixXd X, MatrixXd Z, MatrixXd W, double r, int ngrid){ class_setup(L, R, Rinf, P2, X, Z, W, r, ngrid);}
    
    void solve(){
        iter_=0;
        while (iter_<=maxiter_){
            beta0_ = beta_; gamma0_ = gamma_; lambda0_ = lambda_; puv0_ = puv_; lik0_ = logLik_i_.mean();
            puv_ = calint(beta_, gamma_, lambda_, puv_); lambda_= Mstep(beta_,gamma_);
            iter_++;
            if (((diff(beta_,beta0_)+diff(gamma_,gamma0_)<epsilon_beta)&(diff(lambda_,lambda0_)<epsilon_lambda)&(diff_mat(puv_,puv0_)<epsilon_lambda))|(abs(logLik_i_.mean()-lik0_)<epsilon_lik)) {cov_=true; break;}
            else if (lambda_.cwiseAbs().maxCoeff()>10){
                if (((diff(beta_,beta0_)+diff(gamma_,gamma0_)<epsilon_beta)&(diff_mat(puv_,puv0_)<epsilon_lambda))|(abs(logLik_i_.mean()-lik0_)<epsilon_lik)){cov_=true; break;}
            }
            if (beta_.cwiseAbs().maxCoeff()>maxbound) {cov_=false; break;}
            if (std::isnan(puv_(0,0))==1){
                cout<<"beta= "<<beta0_.transpose()<<endl;
                cout<<"gamma= "<<gamma0_.transpose()<<endl;
                cout<<"lambda = "<<lambda0_.transpose()<<endl;
                cout<<"puv = "<<puv0_.transpose()<<endl;
                cov_=false; break;
            }
            cout<<"iter= "<<iter_<<" beta = "<<beta_.transpose()<<" gamma = "<<gamma_.transpose()<<endl;
            cout<<"lambda = "<<lambda_.head(5).transpose()<<" "<<lambda_.tail(5).transpose()<<endl;
            cout<<", diff_beta="<<diff(beta_,beta0_)+diff(gamma_,gamma0_)<<", diff_lambda="<<diff(lambda_,lambda0_)<<", diff_puv="<<diff_mat(puv_,puv0_)<<",diff_lik = "<<abs(logLik_i_.mean()-lik0_)<<", logLik="<<logLik_i_.mean()<<endl;
        }
        logLik_MLE_i=logLik_i_;
        if (cov_){
            cout<<"The estimate for beta is: "<<beta_.transpose()<<endl;
        }
        else {cout<<"Maximum Iteration Times ("<<maxiter_<<") Reached!"<<endl;}
    }
    
    void est_sd(double h){
        h *= pow((double)n_,-0.5); Profile_I(h);
        Cov_score=pI_score.inverse(); sd_score.resize(p_);
        for (int j=0;j<p_;++j) sd_score(j)=sqrt(Cov_score(j,j));
    }
};

TwoPhase simudata(long nsub, VectorXd beta, VectorXd gamma, default_random_engine & generator, double r, int ngrid){
    //Distributions
    uniform_real_distribution<double> Xdist(-0.5,0.5);
    uniform_real_distribution<double> Udist(0,1);
    uniform_real_distribution<double> Cdist(0,3);
    bernoulli_distribution P2_dist(0.3);
    
    MatrixXd X(nsub,1); MatrixXd Z(nsub,1); MatrixXd W(nsub,2);
    VectorXd L(nsub); VectorXd R(nsub); VectorXi Rinf(nsub);
    VectorXi P2(nsub); P2.setOnes();
    
    double nL=0; double nR=0; double nI=0; int nint=0;
    
    double T, C, expXW, Z2, Z3;
    for (int i=0;i<nsub;++i){
        X(i,0) = Xdist(generator);
        Z2 = Udist(generator); Z3 = floor(Udist(generator) * sin(X(i,0))* 10 + 0.5)/10;
        Z(i,0) = Z3; W.row(i) << Z3, Z2;
        
        expXW = exp(X.row(i) * beta) * exp(W.row(i) * gamma);
        T = genaT(0.1, r, expXW, generator); //Lambda(t) = 0.1t
        C = Cdist(generator);
        vector<double> CC=genC(0, C, 0, 0.4, generator, nint);
        long nC = CC.size();
        if (nC==0){ L(i)=0; R(i)=1; Rinf(i)=1; nR++; }
        else {
            if (T<CC[0]) {L(i)=0; R(i)=CC[0]; Rinf(i)=0;}
            else {
                int j=0;
                while (j<nC){
                    if (T>CC[j]) {L(i)=CC[j]; R(i)=1; Rinf(i)=1;}
                    else {L(i)=CC[j-1]; R(i)=CC[j]; Rinf(i)=0; break;}
                    j++;
                }
            }
            if (Rinf(i)==1) nR++;
            else if (L(i)==0) nL++;
            else nI++;
        }
        if (Rinf(i)==1) P2(i) = P2_dist(generator);
    }
    cout<<"left-censored: "<<nL/(double)nsub<<endl;
    cout<<"interval-censored: "<<nI/(double)nsub<<endl;
    cout<<"right-censored: "<<nR/(double)nsub<<endl;
    
    cout<<"Measured X: "<<(double)P2.sum()/(double)nsub<<endl;
    cout<<"Average #intervals"<<(double)nint / (double)nsub<<endl;
    TwoPhase data(L, R, Rinf, P2, X, Z, W, r, ngrid);
    return(data);
}
#endif /* TwoPhase_h */
