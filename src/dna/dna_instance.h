#pragma once

#include <string>
#include <vector>
#include <mutex>
#include <stdexcept>

class DNAInstance {
public:
    void setSpectrum(const std::vector<std::string>& value) {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (value.empty()) {
            LOG_ERROR("Cannot set empty spectrum");
            throw std::invalid_argument("Cannot set empty spectrum");
        }
        m_spectrum = value;
    }
    
    void clearSpectrum() {
        std::lock_guard<std::mutex> lock(m_mutex);
        LOG_WARNING("Clearing spectrum - this operation should be used with caution");
        m_spectrum.clear();
    }

private:
    std::vector<std::string> m_spectrum;
    std::mutex m_mutex;
}; 