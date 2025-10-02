#ifndef ENVELOPES_H
#define ENVELOPES_H

#include "json.hpp"
#include <functional>
#include <vector>

using json = nlohmann::json;
using EnvelopeFunction = std::function<void(const json&, double*, std::vector<double>&, std::vector<double>&)>;

void off(const json& input, double* tvec, std::vector<double>& env, std::vector<double>& env2);
void constant(const json& input, double* tvec, std::vector<double>& env, std::vector<double>& env2);
void impulse(const json& input, double* tvec, std::vector<double>& env, std::vector<double>& env2);
void gauss(const json& input, double* tvec, std::vector<double>& env, std::vector<double>& env2);
void double_impulse(const json& input, double* tvec, std::vector<double>& env, std::vector<double>& env2);
void double_gauss(const json& input, double* tvec, std::vector<double>& env, std::vector<double>& env2);
#endif