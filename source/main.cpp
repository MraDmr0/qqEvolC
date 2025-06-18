//main executable of qqEvol package 
//to simulate time evolution of 
//D-state quantum system 
//under an external potential

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <complex>
#include "json.hpp"
#include "algorithms.h"
#include "potentials.h"
#include "envelopes.h"

using json = nlohmann::json;
//define possible types for input data
enum FieldType {STRING, INT, FLOAT, BOOLEAN, ARRAY, MATRIX};
//structure name and type for input data
struct FieldRequirement {
    std::string name;
    FieldType type;
};

// determines the type of input and returns corresponding string
std::string getTypeName(const json& value) 
{
    if (value.is_string())          return "string";
    if (value.is_number_integer())  return "integer";
    if (value.is_number())          return "number";
    if (value.is_boolean())         return "boolean";
    if (value.is_array())           return "array";
    if (value.is_object())          return "object";
    if (value.is_null())            return "null";
    return "unknown";
}

//return string of expected data type of input
std::string expectedTypeName(FieldType type) 
{
    switch (type) 
    {
        case STRING:  return "string";
        case INT:     return "integer";
        case FLOAT:   return "number";
        case BOOLEAN: return "boolean";
        case ARRAY:   return "array";
        case MATRIX:  return "matrix";
        default:      return "unknown";
    }
}

//check if expected type coincides with actual type of input
bool isTypeValid(const json& value, FieldType type) 
{
    switch (type) 
    {
        case STRING:  return value.is_string();
        case INT:     return value.is_number_integer();
        case FLOAT:   return value.is_number();
        case BOOLEAN: return value.is_boolean();
        case ARRAY:   return !value.empty() && value.is_array() ;  // extra check below
        case MATRIX:  return !value.empty() && value.is_array() && value.size() > 0  && value[0].is_array();  // basic matrix check
        default:      return false;
    }
}

//check if mandatory data are in input file and if the type or dimension is the expected one
bool validateFields(const json& input, const std::vector<FieldRequirement>& fields) 
{
    // Assumes that basic fields like 'Dstates' and 'qbmode' have been validated before calling.

    int D = input["Dstates"];

    //check if D is non negative
    if (D <= 0) 
    {
        std::cerr << "Invalid dimension 'Dstates': must be positive!\n";
        return false;
    }
    //check if qbmode and D are an allowed combination
    if (input["qbmode"] == "on") //&& D != 2) 
    {
        std::cerr << "Only supported value for 'qbmode' is 'off'!\n";

        //std::cerr << "Only allowed dimension for 'qbmode' = 'on' is 2!\n";
        return false;
    }
    //check if mandatory data are given as input
    for (const auto& field : fields) 
    {
        if (!input.contains(field.name)) 
        {
            std::cerr << "Missing mandatory input data '" << field.name
                      << "' of type '" << expectedTypeName(field.type) << "'!\n";
            return false;
        }
        //if field is in input, check if the type is correct
        const auto& val = input[field.name];

        if (!isTypeValid(val, field.type)) 
        {
            std::cerr << "Wrong type '" << getTypeName(val)
                      << "' for input data '" << field.name
                      << "', expected '" << expectedTypeName(field.type) << "'!\n";
            return false;
        }

        //special case for array data
        if (field.type == ARRAY) 
        {
            //check dimensions and type of elements of "psi"
            if (field.name == "psi") 
            {
                if (val.size() != D) 
                {
                    std::cerr << "Wrong size '" << val.size() << "' for array '" << field.name
                              << ", expected '" << D << "'!\n";
                    return false;
                }
                for (const auto& item : val) 
                {
                    if (!item.is_number()) 
                    {
                        std::cerr << "Wrong type '" << getTypeName(item) << "' for elements of array '" << field.name << "', expected 'number'!\n";
                        return false;
                    }
                }
            }
            //check dimensions and type of elements of "wl" depending on "qbmode"
            else if (field.name == "wl" && input["qbmode"] == "on") 
            {
                if (val.size() != 1) 
                {
                    std::cerr << "Wrong size '"<< val.size() <<"' for array '" << field.name
                              << "', expected '1'!\n";
                    return false;
                }
                for (const auto& item : val) 
                {
                    if (!item.is_number()) 
                    {
                        std::cerr << "Wrong type '" << getTypeName(item) << "' for elements of array '" << field.name << "', expected 'number'!\n";
                        return false;
                    }
                }
            }
            else if (input["qbmode"] == "off") 
            {
                if (val.size() != D) 
                {
                    std::cerr << "Wrong size '" << val.size() << "' for array '" << field.name
                              << "', expected '" << D << "'!\n";
                    return false;
                }
                for (const auto& item : val) 
                {
                    if (!item.is_number()) 
                    {
                        std::cerr << "Wrong type '" << getTypeName(item) << "' for elements of array '" << field.name << "', expected 'number'!\n";
                        return false;
                    }
                }
            }
        }

        //special case for matrix data
        if (field.type == MATRIX) 
        {
            if (val.size() != D) 
            {
                std::cerr << "Wrong number of rows '" << val.size() << "' for matrix '" << field.name
                              << "', expected '" << D << "'!\n";
                
                return false;
            }
            for (const auto& row : val) 
            {
                if (!row.is_array() || row.size() != D) 
                {
                    std::cerr << "Wrong number of columns '" << row.size() << "' for matrix '" << field.name
                              << "', expected '" << D << "'!\n";
                    return false;
                }
                for (const auto& item : row) 
                {
                    if (!item.is_number()) 
                    {
                        std::cerr << "Wrong type '" << getTypeName(item) << "' for elements of matrix '" << field.name << "', expected 'number'!\n";
                        return false;
                    }
                }
            }
        }
    }

    return true;  //all fields pass the check
}


//merge base and optional data
std::vector<FieldRequirement> mergeFields(
    const std::vector<FieldRequirement>& a,
    const std::vector<FieldRequirement>& b)
{
    std::vector<FieldRequirement> result = a;
    result.insert(result.end(), b.begin(), b.end());
    return result;
}

int main(int argc, char* argv[]) 
{
    //check if input file exists
    if (argc < 2) 
    {
        std::cerr << "Missing mandatory argument: json input data file!\n";
        return 1;
    }
    //check if input file can be opened
    std::ifstream file(argv[1]);
    if (!file) 
    {
        std::cerr << "Impossible to open the input file '" << argv[1] << "'!\n";
        return 1;
    }
    //json read input file
    json input;
    file >> input;
    //next checks are needed to determine which are the optional input data
    //check if qbmode is specified in input file
    if (!input.contains("qbmode"))
    {
        std::cerr << "Missing mandatory input data 'qbmode' of type 'string'!\n";
        return 1;
    }
    //check if qbmode is of correct type
    if (!input["qbmode"].is_string())
    {
        std::cerr << "Wrong type '" << getTypeName(input["qbmode"])
        << "' for input data 'qbmode', expected 'string'!\n";
        return 1;           
    }
    if (input["qbmode"] == "on")
    {
        std::cerr << "The specified qbmode " << input["qbmode"] << " is not supported\n";
        return 1;

    }

    //check if envelope is specified in input file
    if (!input.contains("envelope"))
    {
        std::cerr << "Missing mandatory input data 'envelope' of type 'string'!\n";
        return 1;
    }
    //check if envelope is of correct type
    if (!input["envelope"].is_string())
    {
        std::cerr << "Wrong type '" << getTypeName(input["envelope"])
        << "' for input data 'envelope', expected 'string'!\n";
        return 1;           
    }
    //check if Dstates is specified in input file
    if (!input.contains("Dstates"))
    {
        std::cerr << "Missing mandatory input data 'Dstates' of type 'string'!\n";
        return 1;
    }
    //check if Dstates is of correct type
    if (!input["Dstates"].is_number_integer())
    {
        std::cerr << "Wrong type '" << getTypeName(input["Dstates"])
        << "' for input data 'envelope', expected 'integer'!\n";
        return 1;           
    }
    //assign qbmode, envelope, Dstates to local variables
    std::string qbmode   = input["qbmode"];
    std::string envelope = input["envelope"];
    int dimension        = input["Dstates"];


    //map of qbmodes algorithms
    std::unordered_map<std::string, SimulationFunction> qbmodes = 
    {
        {"off", EvolveRK4},
        {"on", EvolveRK4}
    };

    //map of envelope functions
    std::unordered_map<std::string, std::pair<EnvelopeFunction, std::vector<FieldRequirement>>> envelopes = 
    {
        {"off",     {off, {}}},
        {"const",   {constant, {{"F1", FLOAT}}}},
        {"impulse", {impulse, {{"F1", FLOAT},{"t1", FLOAT},{"t2", FLOAT}}}},
        {"gauss" ,  {gauss, {{"F1", FLOAT}, {"t1", FLOAT}, {"sigma1", FLOAT} }}},
        {"double_impulse", {double_impulse, {{"F1", FLOAT},{"t1", FLOAT},{"t2", FLOAT},{"w2", FLOAT},{"t3", FLOAT},{"t4", FLOAT},{"F2",FLOAT}}}},
        {"double_gauss", {double_gauss, {{"F1", FLOAT},{"t1", FLOAT},{"w2", FLOAT},{"F2",FLOAT},{"sigma2", FLOAT}}}}

    };

    //map of potential functions (depend on qbmode and envelope)
    std::unordered_map<std::string, PotentialFunction> potentials = 
    {
        {"off:off",   UpdatePotential},
  
        {"const:off",   UpdatePotential},

        {"impulse:off",   UpdatePotential},

        {"gauss:off",   UpdatePotential},

        {"double_impulse:off", UpdatePotential2},
        {"double_gauss:off",   UpdatePotential2}
    };


    //base mandatory input data
    std::vector<FieldRequirement> baseFields = 
    {
        {"prefix", STRING}, {"qbmode", STRING}, {"envelope", STRING},
        {"Dstates", INT}, {"ti", FLOAT}, {"tf", FLOAT},
        {"Nstep", INT}, {"Nprint", INT}, {"psi", ARRAY}, {"wl", ARRAY}, {"wr", MATRIX}, {"w1", FLOAT}
    };

    //allowed dimensions
    std::unordered_map<int,std::pair<std::string, std::string>> dimensions = 
    {
        {1,{"off","off"}},{2,{"off","off"}},{3,{"off","off"}},{4,{"off","off"}}
    };

    //check if specified qbmode is supported
    if (qbmodes.find(qbmode) == qbmodes.end()) 
    {
        std::cerr << "The specified qbmode '" << qbmode << "' is not supported!\n";
        return 1;
    }
    //check if specified envelope is supported
    if (envelopes.find(envelope) == envelopes.end()) 
    {
        std::cerr << "The specified envelope '" << envelope << "' is not supported!\n";
        return 1;
    }
    //check if specified dimension is supported
    if (dimensions.find(dimension) == dimensions.end()) 
    {
        std::cerr << "The specified dimension '" << dimension << "' is not supported!\n";
        return 1;
    }

    std::string potential;
    potential = envelope + ":" + qbmode ;
    //check if enveope and qbmode specified are compatible
    if (potentials.find(potential) == potentials.end())
    {
        std::cerr << "The specified combination of envelope function and qbmode is not supported!\n";
        return 1;
    }

    //list of needed input data (base+optional)
    const auto& potFields = envelopes[envelope].second;
    auto totalFields = mergeFields(baseFields, potFields);
    //check if needed data are in input file
    if (!validateFields(input, totalFields)) 
    {
        return 1;
    }

    //start simulation
    qbmodes[qbmode](input,  potentials[potential], envelopes[envelope].first);
    return 0;
}
