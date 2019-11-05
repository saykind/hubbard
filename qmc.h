#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include "sqrmtx.h"
#define sgn(q) (q==0 ? +1 : -1)
#define spn(o) (o==1 ?  0 :  1)
#define KronekerDelta(i,j) (i==j ? 1 : 0)

typedef int sgt;
typedef int spt;

class SimQMC {
    private:
        int L = 0;
        int N = 0;
        int NN = 0;
        int l = 0;
        int acc = 0;
        int all = 0;
        int m = 0;
        double beta = 0.;
        double dt = 0.;
        double u = 0.;
        double lmbd = 0.;
        std::vector<std::vector<sgt>> s;
        std::vector<std::vector<dmx>> exp_v;
        std::vector<mtx> g;
        std::vector<double> n;
        std::vector<std::vector<double>> nn;
        mtx exp_k;
        mtx exp_k_;
        std::vector<mtx> G;
    public:
            SimQMC(int L=8, int N=6, double beta=1., double u=4.);
            ~SimQMC();
        mtx _k() const;                             // calculate kinetic energy matrix
        mtx _exp_k(double dt = 1.) const;           // calculate exp of kinetic energy matrix
        dmx _exp_v(spt q = 0, int l = 0) const;     // calculate exp of potential energy matrix
        mtx _g(spt q = 0) const;                    // calculate green's function from the scratch
        void exp_v_(int i0 = -1);                   // update exp_v when i-th spin is flipped
        void g_(int i0 = -1);                       // update green's function when i-th spin is flipped
        void n_();                                  // update occupancy using current configuration s
        void l_();                                  // wrap green's functions and shift l
        void operator++();                          // sweep
        friend std::ostream& operator<<(std::ostream&, const SimQMC&);
        friend std::ofstream& operator<<(std::ofstream&, const SimQMC&);
};

