#include "utils/logging.h"
#include <iostream>

// Initialize static members
std::ofstream Logger::logFile;
std::streambuf* Logger::coutBuffer = nullptr;
std::streambuf* Logger::cerrBuffer = nullptr; 