#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <iostream>
#include <fstream>
#include <string>
#include "algorithms.h"

//executes RK4 simulation for qbmode = off
void EvolveRK4(const json& input, PotentialFunction potential, EnvelopeFunction envelope) 
{
    //Assign base input data to local variables
    std::string prefix = input["prefix"];
    int Nstep          = input["Nstep"];
    int D              = input["Dstates"];
    int Nprint         = input["Nprint"];
    double ti          = input["ti"];
    double tf          = input["tf"];
    //dynamically allocate memory for structured input data 
 
    std::vector<double> wl(D);
    std::vector<std::complex<double>> psi0(D);
    std::vector<std::complex<double>> wr(D*D);

    //dynamically define local variables
    std::vector<std::complex<double>> psiPrev(D, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> psiCurr(D, std::complex<double>(0.0, 0.0));
    std::vector<std::complex<double>> K0(D, 0.0);
    std::vector<std::complex<double>> K1(D, 0.0);
    std::vector<std::complex<double>> K2(D, 0.0);
    std::vector<std::complex<double>> K3(D, 0.0);
        
    //allocate array for envelope function
    std::vector<double> env(3,0);
    std::vector<double> env2(3,0);

    for (int k = 0; k < D ; k++)
    {
        psi0[k] = std::complex<double>(input["psi"][k],0.0);
        psiPrev[k]   = psi0[k];
        psiCurr[k]   = psi0[k]; 

        wl[k]   = input["wl"][k];

        for (int j = 0; j < D ; j++)
        {
            wr[k*D + j] = std::complex<double>(input["wr"][k][j],0.0);
        }
    }

    //allocate vector for time array
    std::vector<double> t(Nstep+1);
    //compute time step
    double dt = (tf-ti)/(double)(Nstep);
    //fill the time array
    for (int i = 0; i < (int)(t.size()); i++)
    {
        t[i] = ti + i*dt;
    }

    std::vector<std::vector<std::complex<double>>> psiOut;

    //allocate potential matrix
    std::vector<std::vector<std::complex<double>>> Vmatrices(3, std::vector<std::complex<double>>(D*D));

    //allocate output time and envelope array
    std::vector<double> tOut;
    std::vector<double> envOut;
    
    //initialize arrays
   
    tOut.push_back(ti);
    psiOut.push_back(psi0);
    //compute envelope at initial time for output
    potential(input, D, t[0], dt,  wl,  wr,  envelope,  Vmatrices , env, env2);
    envOut.push_back(env[0] + env2[0]);

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

        for (int j = 0; j < D; ++j) {

            K1[j] = std::complex<double>(0.0, 0.0);
            for (int k = 0; k < D; ++k) 
            {
                psiCurr[k] = psiPrev[k] + 0.5 * dt * K0[k];
                K1[j] += Vmatrices[1][j * D + k] * psiCurr[k];
            }
        }

        for (int j = 0; j < D; ++j) 
        {
            K2[j] = std::complex<double>(0.0, 0.0);
            for (int k = 0; k < D; ++k) 
            {
                psiCurr[k] = psiPrev[k] + 0.5 * dt * K1[k];
                K2[j] += Vmatrices[1][j * D + k] * psiCurr[k];
            }
        }

        //compute K3
    
        for (int j = 0; j < D; ++j) 
        {
            K3[j] = std::complex<double>(0.0, 0.0);
            for (int k = 0; k < D; ++k) 
            {
                psiCurr[k] = psiPrev[k] + dt * K2[k];
                K3[j] += Vmatrices[2][j * D + k] * psiCurr[k];
            }
        }

        // New psi state
        for (int j = 0; j < D; ++j) 
        {
            psiPrev[j] += (dt / 6.0) * (K0[j] + 2.0 * K1[j] + 2.0 * K2[j] + K3[j]);
        }

        // Normalization
        
        double norm = std::sqrt(
        std::accumulate(
                psiPrev.begin(), psiPrev.end(), 0.0,
                [](double sum, const std::complex<double>& z) {
                    return sum + std::norm(z); // somma dei quadrati dei moduli
                }
            )
        );
        if (norm == 0.0 || std::isnan(norm) || std::isinf(norm)) {
            std::cerr << "Warning: invalid norm at step " << i << "\n";
            norm = 1.0;
        }
        for (int j = 0; j < D; ++j) 
        {
            psiPrev[j] /= norm;
        }

        // Save every Nprint
        if (i % Nprint == 0 || i == Nstep) 
        {
            psiOut.push_back(psiPrev);
            envOut.push_back(env[2]);
            tOut.push_back(t[i]);
        }

    }

    std::cout << "Calculation completed...\n";


    std::string outfile = prefix + ".txt"; 

    std::cout << "Writing output file...\n";
    
    //print out result
    FILE* f = std::fopen(outfile.c_str(), "w");
if (!f) { std::cerr<<"fopen failed\n"; }
else {

    for (int i = 0; i < (int)(tOut.size()); ++i) {
        std::fprintf(f, "%g %g ", tOut[i], envOut[i]);
        for (int k = 0; k < D; ++k) {
            std::fprintf(f, "%g+%gj ", std::real(psiOut[i][k]), std::imag(psiOut[i][k]));
        }
        std::fprintf(f, "\n");
    }
    std::fclose(f);
}

    std::cout << "Output file written correctly...\n";

}
