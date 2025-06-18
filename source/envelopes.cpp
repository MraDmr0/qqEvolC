#include "envelopes.h"
#include <iostream>
#include <cmath> 

//potential constant and equal to 0
void off(const json& input, double* tvec, double* env, double* env2)
{
    for (int k = 0; k < 3; k++)
    {
        env[k] = 0.0;
    }
}

//potential constant and equal to F1
void constant(const json& input, double* tvec, double* env, double* env2)
{
    
    double F1 = input["F1"];
    for (int k = 0; k < 3; k++)
    {
        env[k] = F1;
    }
    
}

//square potential between t1 and t2 equal to F1
void impulse(const json& input, double* tvec, double* env, double* env2)
{
    
    double F1 = input["F1"];
    double t1 = input["t1"];
    double t2 = input["t2"];
    
    for (int k = 0; k < 3; k++)
    {
        if (tvec[k] < t1 || tvec[k] > t2)
        {
            env[k] = 0.0;
        }
        else
        {
            env[k] = F1;
        }
    }
}

void double_impulse(const json& input, double* tvec, double* env, double* env2)
{
    
    double F1 = input["F1"];
    double t1 = input["t1"];
    double t2 = input["t2"];

    double F2 = input["F2"];
    double t3 = input["t3"];
    double t4 = input["t4"];
    
    for (int k = 0; k < 3; k++)
    {
        if (tvec[k] < t1 || tvec[k] > t2)
        {
            env[k] = 0.0;
        }
        else
        {
            env[k] = F1;
        }
    }

    for (int k = 0; k < 3; k++)
    {
        if (tvec[k] < t3 || tvec[k] > t4)
        {
            env[k] = 0.0;
        }
        else
        {
            env2[k] = F2;
        }
    }
}

//gussian potential centered in t1, strength F1 and amplitude sigma1
void gauss(const json& input, double* tvec, double* env, double* env2)
{
    
    double F1 = input["F1"];
    double t1 = input["t1"];
    double sigma1 = input["sigma1"];

    for (int k = 0; k < 3; k++)
    {
        env[k] = F1/(std::sqrt(2.0*M_PI)*sigma1)*std::exp(-pow((tvec[k]-t1)*1.0e6/sigma1,2.0)/2.0);
    }

}


void double_gauss(const json& input, double* tvec, double* env, double* env2)
{
    
    double F1 = input["F1"];
    double t1 = input["t1"];
    double sigma1 = input["sigma1"];

    double F2 = input["F2"];
    double t2 = input["t2"];
    double sigma2 = input["sigma2"];


    for (int k = 0; k < 3; k++)
    {
        env[k] = F1/(std::sqrt(2.0*M_PI)*sigma1)*std::exp(-pow((tvec[k]-t1)*1.0e6/sigma1,2.0)/2.0);
        env2[k] = F2/(std::sqrt(2.0*M_PI)*sigma2)*std::exp(-pow((tvec[k]-t2)*1.0e6/sigma2,2.0)/2.0);

    }

}

