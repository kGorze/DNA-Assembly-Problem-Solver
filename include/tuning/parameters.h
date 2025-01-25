#pragma once

#include <string>
#include <vector>
#include <map>
#include <memory>

/**
 * Structure representing a tunable parameter
 */
struct Parameter {
    std::string name;
    double minValue;
    double maxValue;
    double currentValue;
    bool isInteger;
    
    Parameter(const std::string& n, double min, double max, double current, bool integer = false)
        : name(n), minValue(min), maxValue(max), currentValue(current), isInteger(integer) {}
};

/**
 * Class managing tunable parameters
 */
class Parameters {
public:
    Parameters() = default;
    
    /**
     * Add a new parameter
     */
    void addParameter(const std::string& name, double min, double max, double initial, bool isInteger = false);
    
    /**
     * Get parameter by name
     */
    Parameter& getParameter(const std::string& name);
    
    /**
     * Get all parameters
     */
    std::vector<Parameter>& getParameters() { return m_parameters; }
    
    /**
     * Set parameter value
     */
    void setParameterValue(const std::string& name, double value);
    
    /**
     * Get parameter value
     */
    double getParameterValue(const std::string& name) const;
    
private:
    std::vector<Parameter> m_parameters;
    std::map<std::string, size_t> m_parameterIndices;
}; 