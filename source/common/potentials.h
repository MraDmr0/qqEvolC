#ifndef POTENTIALS_H
#define POTENTIALS_H

#include "json.hpp"
#include <functional>
#include <complex>

using json = nlohmann::json;
using EnvelopeFunction = std::function<void(const json&, double*, double*, double*)>;
using PotentialFunction = std::function<void(const json& , int, double, double , double* , std::complex<double>** , EnvelopeFunction , std::complex<double>**, double*, double*)>;

void UpdatePotential(const json& input, int D, double t, double dt, double* wl, std::complex<double>** wr, EnvelopeFunction envelope, std::complex<double>** Vmatrices , double* env, double* env2);
void UpdatePotential2(const json& input, int D, double t, double dt, double* wl, std::complex<double>** wr, EnvelopeFunction envelope, std::complex<double>** Vmatrices , double* env, double* env2);
#endif