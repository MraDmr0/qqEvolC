#ifndef SIMULATION_H
#define SIMULATION_H

#include "json.hpp"
#include <functional>
#include <complex>

using json = nlohmann::json;
using EnvelopeFunction  = std::function<void(const json&, double*, double*, double*)>;
using PotentialFunction = std::function<void(const json& , int, double, double , double* , std::complex<double>** , EnvelopeFunction , std::complex<double>**, double*, double*)>;
using SimulationFunction = std::function<void(const json&, PotentialFunction, EnvelopeFunction)>;

void EvolveRK4(const json& input, PotentialFunction potential, EnvelopeFunction envelope);
#endif
