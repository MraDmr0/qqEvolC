#ifndef POTENTIALS_H
#define POTENTIALS_H

#include "json.hpp"
#include <functional>
#include <complex>
#include <vector>

using json = nlohmann::json;
using EnvelopeFunction = std::function<void(const json&, double*, std::vector<double>&, std::vector<double>&)>;
using PotentialFunction = std::function<void(const json& , int, double, double , std::vector<double>& , std::vector<std::complex<double>>& , EnvelopeFunction , std::vector<std::vector<std::complex<double>>>&, std::vector<double>&, std::vector<double>&)>;

void UpdatePotential(const json& input, int D, double t, double dt, std::vector<double>& wl, std::vector<std::complex<double>>& wr, EnvelopeFunction envelope, std::vector<std::vector<std::complex<double>>>& Vmatrices , std::vector<double>& env, std::vector<double>& env2);
void UpdatePotential2(const json& input, int D, double t, double dt, std::vector<double>& wl, std::vector<std::complex<double>>& wr, EnvelopeFunction envelope, std::vector<std::vector<std::complex<double>>>& Vmatrices , std::vector<double>& env, std::vector<double>& env2);
#endif