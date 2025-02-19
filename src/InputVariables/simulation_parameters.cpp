/* This file is part of SIRIUS electronic structure library.
 *
 * Copyright (c), ETH Zurich.  All rights reserved.
 *
 * Please, refer to the LICENSE file in the root directory.
 * SPDX-License-Identifier: BSD-3-Clause
 */
#include <unordered_set>
#include "GlobalFunctions.hpp"
#include "simulation_parameters.hpp"
#include "Json/json_sirius.hpp"

/// Compose JSON dictionary with default parameters based on input schema.
/** Traverse the JSON schema and add nodes with default parameters to the output dictionary. The nodes without
 *  default parameters are ignored. Still, user has a possibility to add the missing nodes later by providing a
 *  corresponding input JSON dictionary. See compose_json() function. */
void
compose_default_json(nlohmann::json const& schema__, nlohmann::json& output__)
{
    for (auto it : schema__.items()) {
        auto key = it.key();
        /* this is a final node with the description of the data type */
        if (it.value().contains("type") &&  it.value()["type"] != "object") {
            /* check if default parameter is present */
            if (it.value().contains("default")) {
                output__[key] = it.value()["default"];
            }
            else if(it.value().contains("items")) {
                if (!output__.contains(key)) {
                    output__[key] = nlohmann::json{};
                    compose_default_json(it.value()["items"], output__[key]);
                }
            }
        } else { /* otherwise continue to traverse the schema */
            if (!output__.contains(key)) {
                output__[key] = nlohmann::json{};
            }
            if (it.value().contains("properties")) {
                compose_default_json(it.value()["properties"], output__[key]);
            } 
        }
    }
}

/// Append the input dictionary to the existing dictionary.
/** Use JSON schema to traverse the existing dictionary and add on top the values from the input dictionary. In this
 *  way we can add missing nodes which were not defined in the existing dictionary. */
void
compose_json(nlohmann::json const& schema__, nlohmann::json const& in__, nlohmann::json& inout__)
{
    std::unordered_set<std::string> visited;

    for (auto it : in__.items()) {
        visited.insert(it.key());
    }

    for (auto it : schema__.items()) {
        auto key = it.key();

        // Remove visited items.
        auto found = visited.find(key);
        if (found != visited.end()) {
            visited.erase(found);
        }

        /* this is a final node with the description of the data type */
        if (it.value().contains("type") && it.value()["type"] != "object") {
            if (in__.contains(key)) {
                /* copy the new input */
                inout__[key] = in__[key];
            }
        } else { /* otherwise continue to traverse the schema */
            /* not simple data type : a section with parameter, a dictionary map, etc.*/
            if (it.value().contains("properties")) {
                compose_json(it.value()["properties"], in__.contains(key) ? in__[key] : nlohmann::json{}, inout__[key]);
            } else if (in__.contains(key)) {
                inout__[key] = in__[key];
            } else if (!inout__.contains(key)) {
                inout__[key] = nlohmann::json();
            }
        }
    }

    // Emit warnings about keys that were set but unused.
    if (!visited.empty()) {
        std::stringstream ss;
        ss << "The following configuration parameters were not recognized and ignored: ";
        std::copy(visited.begin(), visited.end(), std::ostream_iterator<std::string>(ss, " "));
        output::print(ss.str());
    }
}

Config::Config()
{
    /* initialize JSON dictionary with default parameters */
    compose_default_json(EDUS::input_schema["properties"], this->dict_);
}

void
Config::import(nlohmann::json const& in__)
{
    /* overwrite the parameters by the values from the input dictionary */
    compose_json(EDUS::input_schema["properties"], in__, this->dict_);
}

void
Simulation_parameters::import(nlohmann::json const& dict__)
{
    cfg_.import(dict__);
}

void
Simulation_parameters::import(std::string const& str__)
{
    auto dict = read_json_from_file_or_string(str__);
    this->import(dict);
}

