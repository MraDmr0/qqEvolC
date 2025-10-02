#ifndef SIMULATION_H
#define SIMULATION_H

#include "json.hpp"
#include <functional>
#include <complex>
#include <vector>

using json = nlohmann::json;
using EnvelopeFunction  = std::function<void(const json&, double*, std::vector<double>& , std::vector<double>&)>;
using PotentialFunction = std::function<void(const json& , int, double, double , std::vector<double>& , std::vector<std::complex<double>>& , EnvelopeFunction , std::vector<std::vector<std::complex<double>>>& , std::vector<double>&, std::vector<double>&)>;
using SimulationFunction = std::function<void(const json&, PotentialFunction, EnvelopeFunction)>;

void EvolveRK4(const json& input, PotentialFunction potential, EnvelopeFunction envelope);
#endif
