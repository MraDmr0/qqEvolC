#include "potentials.h"
#include <iostream>
#include <complex>


const double hbar = 6.582119569e-13;

//update potential matrix for qbmode = off
void UpdatePotential(const json& input, int D, double t, double dt, double* wl, std::complex<double>** wr, EnvelopeFunction envelope, std::complex<double>** Vmatrices , double* env, double* env2)
{
    //env should be a vector with envelope function at times [t, t + 0.5dt, t+dt] 
    assert(env != nullptr);
    assert(Vmatrices != nullptr);
    assert(wl != nullptr);
    assert(wr != nullptr);

    double w1 = input["w1"];
    std::complex<double> im;
    im = std::complex<double>(0.0, 1.0);
    
    double tvec[3];
    tvec[0] = t;
    tvec[1] = t + 0.5*dt;
    tvec[2] = t + dt;


    envelope(input, tvec, env, env2);
    
    for (int k = 0; k < 3; k++) 
    {
        for (int i = 0; i < D; i++) 
        {
            for (int j = 0; j < D; j++) 
            {
                Vmatrices[k][i*D + j] = -im*env[k]*wr[i][j] * std::exp((im/hbar)*(wl[i]-wl[j])*tvec[k])*std::cos(w1*tvec[k]);
            }
        }
    }

    
}

void UpdatePotential2(const json& input, int D, double t, double dt, double* wl, std::complex<double>** wr, EnvelopeFunction envelope, std::complex<double>** Vmatrices , double* env, double* env2)
{
    //env should be a vector with envelope function at times [t, t + 0.5dt, t+dt] 
    assert(env != nullptr);
    assert(Vmatrices != nullptr);
    assert(wl != nullptr);
    assert(wr != nullptr);

    double w1 = input["w1"];
    double w2 = input["w2"];
    std::complex<double> im;
    im = std::complex<double>(0.0, 1.0);

    
    std::complex<double>** Vmatrices1 = new std::complex<double>*[3];
    
    for (int i = 0; i < 3; i++) 
    {
        Vmatrices1[i] = new std::complex<double>[D * D];
    }
    
    double tvec[3];
    tvec[0] = t;
    tvec[1] = t + 0.5*dt;
    tvec[2] = t + dt;

    envelope(input, tvec, env, env2);

    for (int k = 0; k < 3; k++) 
    {
        for (int i = 0; i < D; i++) 
        {
            for (int j = 0; j < D; j++) 
            {
                Vmatrices1[k][i*D + j] = -im*env[k]*wr[i][j] * std::exp((im/hbar)*(wl[i]-wl[j])*tvec[k])*std::cos(w1*tvec[k]);
                Vmatrices[k][i*D + j] = -im*env2[k]*wr[i][j] * std::exp((im/hbar)*(wl[i]-wl[j])*tvec[k])*std::cos(w2*tvec[k])+ Vmatrices1[k][i*D+j];
            }
        }
    }

    env[2] += env2[2];
    for (int i = 0; i < 3; ++i) 
    {
        delete[] Vmatrices1[i];
    }
    delete[] Vmatrices1;
 

}
