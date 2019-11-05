#include "qmc.h"
#include "sqrmtx.h"
#include <vector>
#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

// Information Printout
std::ostream& operator<<(std::ostream& stream, const SimQMC& s) {
    stream <<    "L=" << s.L            << '\t';
    stream <<   "NN=" << s.NN           << '\t';
    stream << "beta=" << s.beta         << '\t';
    stream <<    "u=" << s.u            << '\t';
    stream <<    "M=" << s.m            << '\t';
    stream << std::fixed << setprecision(4);
    stream <<  "acc=" << 1.*s.acc/s.all << endl;
    stream <<   "nu=" << 1.*s.n[0]/s.m  << " +/- " << 1.*s.nn[0][0]/s.m/sqrt(s.m) << '\t';
    stream <<   "nd=" << 1.*s.n[1]/s.m  << " +/- " << 1.*s.nn[1][1]/s.m/sqrt(s.m) << '\t';
    stream <<  "nu2=" << 1.*s.nn[0][0]/s.m    << '\t';
    stream <<  "nd2=" << 1.*s.nn[1][1]/s.m    << '\t';
    stream <<   "nn=" << 1.*s.nn[0][1]/s.m    << '\t';
    stream << setprecision(6) << endl;
    stream.unsetf(ios_base::floatfield);
    return stream;
}

// Initialization
SimQMC::SimQMC(int L, int N, double beta, double u) : 
    L(L), N(N), NN(N*N), l(0),
    beta(beta), dt(beta/L), u(u), 
    lmbd(acosh(exp(dt*u/2))),
    s(L, vector<sgt>(NN)),
    exp_v(2, vector<dmx>(L, dmx(NN))),
    g(2, mtx(NN)), 
    n(2, 0), nn(2, vector<double>(2,0)),
    exp_k(_exp_k(dt)), 
    exp_k_(_exp_k(-dt)) 
{
    srand(time(nullptr));
    #pragma omp parallel for 
    for (int l=0; l<L; l++)
        for (sgt& o : s[l])
            o = rand()%2*2-1;
    exp_v_();
    g_();
}

SimQMC::~SimQMC() {
    std::ofstream file;
    file.open("data.dat", ios::app);
    double m2 = (n[0]+n[1]-2*nn[0][1])/m;
    file << setprecision(6);
    file << beta << '\t' << u << '\t' << m2 << '\n';
    file.close();
    std::cout << "1/T=" << beta << '\t' << "u=" << u << '\t' << "m2=" << m2 << endl;
}

// Calculation functions
mtx SimQMC::_k() const {
   mtx r(NN);
   #pragma omp parallel for collapse(2) 
   for (int i=0; i<NN; i++)
       for (int j=0; j<NN; j++) {
           double sum = 0;
           for (int k=0; k<NN; k++)
               sum += (1./NN) * 2.*(cos((2*M_PI/N)*(k/N))+cos((2*M_PI/N)*(k%N))) 
                   * cos((2*M_PI/N)*( (i/N-j/N)*(k/N) + (i%N-j%N)*(k%N) ));
           r[i][j] = sum;
       }
   return r.chp();
}

mtx SimQMC::_exp_k(double dt) const {
    mtx r(NN);
    #pragma omp parallel for collapse(2) 
    for (int i=0; i<NN; i++)
        for (int j=0; j<NN; j++) {
            double s = 0;
            for (int k=0; k<NN; k++)
                s += (1./NN) * exp(2.*dt*(cos((2*M_PI/N)*(k/N))+cos((2*M_PI/N)*(k%N))))
                    * cos((2*M_PI/N)*( (i/N-j/N)*(k/N) + (i%N-j%N)*(k%N) ));
            r[i][j] = s;
        }
    return r.chp();
}

dmx SimQMC::_exp_v(spt q, int l) const {
    sgt o = sgn(q);
    dmx exp_v(NN);
    for (int i=0; i<NN; i++)
        exp_v[i] = exp(o*lmbd*s[l][i]);
    return exp_v;
}

mtx SimQMC::_g(spt q) const {
    mtx g(NN,1);
    for (int l=L-1; l>=0; l--) {
        g *= exp_k_*exp_v[q][l];
    }
    g += 1;
    return g.inv();
}

// Update Functions
void SimQMC::exp_v_(int i0) {
    if (i0 < 0) {
        for (spt q : {0,1})
            #pragma omp parallel for collapse(2) 
            for (int l=0; l<L; l++)
                for (int i=0; i<NN;i++)
                    exp_v[q][l][i] = exp(sgn(q)*lmbd*s[l][i]);
        return;
    }
    #pragma omp parallel for num_threads(2)
    for (spt q=0; q<2; q++)
        exp_v[q][l][i0] *= exp(-2*sgn(q)*lmbd*s[l][i0]);
}

void SimQMC::g_(int i0) {
    if (i0 < 0) {
        #pragma omp parallel for num_threads(2)
        for (spt q=0; q<2; q++)
            g[q] = _g(q);
        return;
    }
    #pragma omp parallel for num_threads(2)
    for (spt q=0; q<2; q++) {
        sgt o = sgn(q);
        mtx dg(NN);
        double D = (exp(-2*lmbd*o*s[l][i0])-1)/(1+(1-g[q][i0][i0])*(exp(-2*lmbd*o*s[l][i0])-1));
        #pragma omp parallel for collapse(2)
        for (int i=0; i<NN; i++)
            for (int j=0; j<NN; j++)
                dg[i][j] = D*(KronekerDelta(i,i0) - g[q][i][i0])*g[q][i0][j];
        g[q] -= dg;
    }
}

void SimQMC::n_() {
    for (int i=0; i<NN; i++) {
        #pragma omp parallel for num_threads(2)
        for (spt q=0; q<2; q++) {
            n[q] += (1.-g[q][i][i])/NN;
            #pragma omp parallel for num_threads(2)
            for (spt p=0; p<2; p++)
                nn[q][p] += (1.-g[q][i][i])*(1.-g[p][i][i])/NN;
        }
    }
    m++;
}

void SimQMC::l_() {
    #pragma omp parallel for num_threads(2)
    for (spt q=0; q<2; q++)
        g[q] = exp_k_*exp_v[q][l]*g[q]*exp_v[(q+1)%2][l]*exp_k;
    l = (l+1)%L;
}

// Monte Carlo sweep
void SimQMC::operator++() {
    int i = rand()%NN;
    double r = 1.*rand()/RAND_MAX;
    double p = (1+(1-g[0][i][i])*(exp(-2*lmbd*s[l][i])-1))*(1+(1-g[1][i][i])*(exp(2*lmbd*s[l][i])-1));
    if (r < p) {
        exp_v_(i);
        g_();
        s[l][i] *= -1; 
        acc++;
    }
    all++;
}


