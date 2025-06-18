#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <iostream>
#include "algorithms.h"

//executes RK4 simulation for qbmode = off
void EvolveRK4(const json& input, PotentialFunction potential, EnvelopeFunction envelope) 
{
    //Assign base input data to local variables
    int Nstep    = input["Nstep"];
    int D        = input["Dstates"];
    int Nprint   = input["Nprint"];
    double ti    = input["ti"];
    double tf    = input["tf"];
    //dynamically allocate memory for structured input data 
    double* wl = new double[D];
    std::complex<double>* psi0 = new std::complex<double>[D];
    
    std::complex<double>** wr = new std::complex<double>*[D];
    for (int k = 0; k < D ; k++)
    {
        wr[k] = new std::complex<double>[D];
    }

    for (int k = 0; k < D ; k++)
    {
        psi0[k] = std::complex<double>(input["psi"][k],0.0);
        wl[k]   = input["wl"][k];

        for (int j = 0; j < D ; j++)
        {
            wr[k][j] = std::complex<double>(input["wr"][k][j],0.0);
        }
    }
    //allocate vector for time array
    double* t   = new double[Nstep+1];
    //compute time step
    double dt = (tf-ti)/(double)(Nstep);
    //fill the time array
    t[0] = ti;
    for (int i = 1; i < Nstep+1; i++)
    {
        t[i] = i*dt;
    }
    //dynamically define local variables
    std::complex<double>* psiPrev = new std::complex<double>[D];
    std::complex<double>* psiCurr = new std::complex<double>[D];
    std::complex<double>** psiOut = new std::complex<double>*[D];
    std::complex<double>* K0 = new std::complex<double>[D];
    std::complex<double>* K1 = new std::complex<double>[D];
    std::complex<double>* K2 = new std::complex<double>[D];
    std::complex<double>* K3 = new std::complex<double>[D];
    
    //allocate array for envelope function
    double* env = new double[3];
    double* env2 = new double[3];

    for (int k = 0; k < 3; k++)
    {
        env[k] = 0.0;
        env2[k] = 0.0;

    }
    //compute how many points are printed in output
    int Nsave;
    if (Nstep%Nprint != 0.0) 
    {
        Nsave = Nstep/Nprint + 2;
    }
    else
    {
        Nsave = Nstep/Nprint + 1;
    }
    //allocate potential matrix
    std::complex<double>** Vmatrices = new std::complex<double>*[3];
    for (int i = 0; i < 3; i++) 
    {
        Vmatrices[i] = new std::complex<double>[D * D];
    }

    //allocate output time and envelope array
    double* tOut = new double[Nsave];
    double* envOut = new double [Nsave];
    for (int k = 0; k < D; k++)
    {
        psiOut[k] = new std::complex<double>[Nsave];
    }
    
    //initialize arrays
    for (int i = 0; i < Nsave; i++)
    {
        envOut[i] = 0.0;
        tOut[i]   = 0.0;
    }
    tOut[0] = ti;

    for (int k = 0; k < D; k++)
    {
        psiPrev[k]   = psi0[k];
        psiCurr[k]   = psi0[k];
        psiOut[k][0] = psi0[k];  
    }
    //compute envelope at initial time for output
    potential(input, D, t[0], dt,  wl,  wr,  envelope,  Vmatrices , env, env2);
    envOut[0] = env[2] + env2[2];
    //index for output data
    int idx = 1;

    for (int i = 1; i < Nstep+1 ; i++)
    {
        //update potential
        potential(input, D, t[i-1], dt,  wl,  wr,  envelope,  Vmatrices , env, env2);

        //compute K0
        for (int j = 0; j < D; ++j) 
        {
            K0[j] = std::complex<double>(0.0, 0.0);
            for (int k = 0; k < D; ++k) 
            {
                K0[j] += Vmatrices[0][j*D + k] * psiPrev[k];
            }
        }

        //compute K1
        for (int j = 0; j < D; ++j) 
        {
            psiCurr[j] = psiPrev[j] + 0.5 * dt * K0[j];
        }
        for (int j = 0; j < D; ++j) {

            K1[j] = std::complex<double>(0.0, 0.0);
            for (int k = 0; k < D; ++k) 
            {
                K1[j] += Vmatrices[1][j * D + k] * psiCurr[k];
            }
        }

        //compute K2
        for (int j = 0; j < D; ++j) 
        {
            psiCurr[j] = psiPrev[j] + 0.5 * dt * K1[j];
        }
        for (int j = 0; j < D; ++j) 
        {
            K2[j] = std::complex<double>(0.0, 0.0);
            for (int k = 0; k < D; ++k) 
            {
                K2[j] += Vmatrices[1][j * D + k] * psiCurr[k];
            }
        }

        //compute K3
        for (int j = 0; j < D; ++j) 
        {
            psiCurr[j] = psiPrev[j] + dt * K2[j];
        }
        for (int j = 0; j < D; ++j) 
        {
            K3[j] = std::complex<double>(0.0, 0.0);
            for (int k = 0; k < D; ++k) 
            {
                K3[j] += Vmatrices[2][j * D + k] * psiCurr[k];
            }
        }

        // New psi state
        for (int j = 0; j < D; ++j) 
        {
            psiPrev[j] += (dt / 6.0) * (K0[j] + 2.0 * K1[j] + 2.0 * K2[j] + K3[j]);
        }

        // Normalization
        std::complex<double> norm = 0.0;
        for (int j = 0; j < D; ++j) 
        {
            norm += std::conj(psiPrev[j]) * psiPrev[j];
        }
        norm = std::sqrt(norm);
        for (int j = 0; j < D; ++j) 
        {
            psiPrev[j] /= norm;
        }

        // Save every Nprint
        if (i % Nprint == 0 || i == Nstep) 
        {
            for (int j = 0; j < D; ++j) 
            {
                psiOut[j][idx] = psiPrev[j];
            }
            envOut[idx]    = env[2];
            tOut[idx] = t[i];
            idx += 1;
        }

    }
    //print out result
    for (int i = 0; i < Nsave; i++)
    {
        std::cout << tOut[i]<<" " << envOut[i] << " ";
        for (int k = 0; k < D; k++)
        {
            std::cout << real( psiOut[k][i]) << "+" << imag( psiOut[k][i]) << "j"<< " ";
        }
        std::cout << "\n";
    }

    //Erase dynamically allocated memory
    delete [] t;
    delete [] tOut;
    delete [] psiPrev;
    delete [] psiCurr;
    delete [] K0;
    delete [] K1;
    delete [] K2;
    delete [] K3;
    for (int k = 0; k < D; k++)
    {
        delete [] psiOut[k];
        delete[] wr[k];
    }
    delete [] psiOut;
    delete [] psi0;
    delete [] envOut;
    delete [] wr;
    delete [] wl; 
    for (int i = 0; i < 3; ++i) 
    {
        delete[] Vmatrices[i];
    }
    delete[] Vmatrices;
    delete[] env;
    delete[] env2;

}




