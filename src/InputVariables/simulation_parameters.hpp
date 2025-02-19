#ifndef SIMULATION_PARAMETERS_HPP
#define SIMULATION_PARAMETERS_HPP

#include "Json/json.hpp"
#include "input_schema.hpp"
#include "config.hpp"

void compose_default_json(nlohmann::json const& schema__, nlohmann::json& output__);
void compose_json(nlohmann::json const& schema__, nlohmann::json const& in__, nlohmann::json& inout__);

class Config : public config_t
{
    public: 
        Config();
        void import(nlohmann::json const& in__);
};


/// Set of basic parameters of a simulation.
/** This class provides shortcuts to the mostly used input parameters, for example `verbosity`. */
class Simulation_parameters
{
  private:
    /// All user-provided paramters are stored here.
    Config cfg_;

  public:
    Config&
    cfg()
    {
        return cfg_;
    }

    Config const&
    cfg() const
    {
        return cfg_;
    }

    Simulation_parameters()
    {
    }

    /// Import parameters from a file or a serialized json string.
    void
    import(std::string const& str__);

    /// Import parameters from a json dictionary.
    void
    import(nlohmann::json const& dict__);
};

#endif