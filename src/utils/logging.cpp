#include "../../include/utils/logging.h"
#include <iostream>
#include <filesystem>

// Define static members
std::ofstream Logger::logFile;
std::streambuf* Logger::coutBuffer = nullptr;
std::streambuf* Logger::cerrBuffer = nullptr;
bool Logger::initialized = false;

// Initialize static members
void Logger::initialize() {
    if (!initialized) {
        logFile.open("log.txt");
        coutBuffer = std::cout.rdbuf();
        cerrBuffer = std::cerr.rdbuf();
        std::cout.rdbuf(logFile.rdbuf());
        std::cerr.rdbuf(logFile.rdbuf());
        initialized = true;
    }
}

void Logger::uninitialize() {
    if (initialized) {
        std::cout.rdbuf(coutBuffer);
        std::cerr.rdbuf(cerrBuffer);
        logFile.close();
        initialized = false;
    }
}

void Logger::log(const std::string& message) {
    if (initialized) {
        logFile << message << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category) {
    if (initialized) {
        logFile << "[" << category << "] " << message << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << " and " << extra2 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", and " << extra3 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", and " << extra4 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", and " << extra5 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", and " << extra6 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", and " << extra7 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", and " << extra8 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", and " << extra9 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", " << extra9 << ", and " << extra10 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10, const std::string& extra11) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", " << extra9 << ", " << extra10 << ", and " << extra11 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10, const std::string& extra11, const std::string& extra12) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", " << extra9 << ", " << extra10 << ", " << extra11 << ", and " << extra12 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10, const std::string& extra11, const std::string& extra12, const std::string& extra13) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", " << extra9 << ", " << extra10 << ", " << extra11 << ", " << extra12 << ", and " << extra13 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10, const std::string& extra11, const std::string& extra12, const std::string& extra13, const std::string& extra14) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", " << extra9 << ", " << extra10 << ", " << extra11 << ", " << extra12 << ", " << extra13 << ", and " << extra14 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10, const std::string& extra11, const std::string& extra12, const std::string& extra13, const std::string& extra14, const std::string& extra15) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", " << extra9 << ", " << extra10 << ", " << extra11 << ", " << extra12 << ", " << extra13 << ", " << extra14 << ", and " << extra15 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10, const std::string& extra11, const std::string& extra12, const std::string& extra13, const std::string& extra14, const std::string& extra15, const std::string& extra16) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", " << extra9 << ", " << extra10 << ", " << extra11 << ", " << extra12 << ", " << extra13 << ", " << extra14 << ", " << extra15 << ", and " << extra16 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10, const std::string& extra11, const std::string& extra12, const std::string& extra13, const std::string& extra14, const std::string& extra15, const std::string& extra16, const std::string& extra17) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", " << extra9 << ", " << extra10 << ", " << extra11 << ", " << extra12 << ", " << extra13 << ", " << extra14 << ", " << extra15 << ", " << extra16 << ", and " << extra17 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10, const std::string& extra11, const std::string& extra12, const std::string& extra13, const std::string& extra14, const std::string& extra15, const std::string& extra16, const std::string& extra17, const std::string& extra18) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", " << extra9 << ", " << extra10 << ", " << extra11 << ", " << extra12 << ", " << extra13 << ", " << extra14 << ", " << extra15 << ", " << extra16 << ", " << extra17 << ", and " << extra18 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10, const std::string& extra11, const std::string& extra12, const std::string& extra13, const std::string& extra14, const std::string& extra15, const std::string& extra16, const std::string& extra17, const std::string& extra18, const std::string& extra19) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", " << extra9 << ", " << extra10 << ", " << extra11 << ", " << extra12 << ", " << extra13 << ", " << extra14 << ", " << extra15 << ", " << extra16 << ", " << extra17 << ", " << extra18 << ", and " << extra19 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10, const std::string& extra11, const std::string& extra12, const std::string& extra13, const std::string& extra14, const std::string& extra15, const std::string& extra16, const std::string& extra17, const std::string& extra18, const std::string& extra19, const std::string& extra20) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", " << extra9 << ", " << extra10 << ", " << extra11 << ", " << extra12 << ", " << extra13 << ", " << extra14 << ", " << extra15 << ", " << extra16 << ", " << extra17 << ", " << extra18 << ", " << extra19 << ", and " << extra20 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10, const std::string& extra11, const std::string& extra12, const std::string& extra13, const std::string& extra14, const std::string& extra15, const std::string& extra16, const std::string& extra17, const std::string& extra18, const std::string& extra19, const std::string& extra20, const std::string& extra21) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", " << extra9 << ", " << extra10 << ", " << extra11 << ", " << extra12 << ", " << extra13 << ", " << extra14 << ", " << extra15 << ", " << extra16 << ", " << extra17 << ", " << extra18 << ", " << extra19 << ", " << extra20 << ", " << extra21 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10, const std::string& extra11, const std::string& extra12, const std::string& extra13, const std::string& extra14, const std::string& extra15, const std::string& extra16, const std::string& extra17, const std::string& extra18, const std::string& extra19, const std::string& extra20, const std::string& extra21, const std::string& extra22) {
    if (initialized) {
        logFile << "[" << category << "] " << message << " (" << file << ":" << line << " in " << function << " with " << extra << ", " << extra2 << ", " << extra3 << ", " << extra4 << ", " << extra5 << ", " << extra6 << ", " << extra7 << ", " << extra8 << ", " << extra9 << ", " << extra10 << ", " << extra11 << ", " << extra12 << ", " << extra13 << ", " << extra14 << ", " << extra15 << ", " << extra16 << ", " << extra17 << ", " << extra18 << ", " << extra19 << ", " << extra20 << ", " << extra21 << ", " << extra22 << ")" << std::endl;
    }
}

void Logger::log(const std::string& message, const std::string& category, const std::string& file, int line, const std::string& function, const std::string& extra, const std::string& extra2, const std::string& extra3, const std::string& extra4, const std::string& extra5, const std::string& extra6, const std::string& extra7, const std::string& extra8, const std::string& extra9, const std::string& extra10, const std::string& extra11, const std::string& extra12, const std::string& extra13, const std::string& extra14, const std::string& extra15, const std::string& extra16, const std::string& extra17, const std::string& extra18, const std::string& extra19, const std::string& extra20, const std::string& extra21, const std::string& extra22, const std::string& extra23) {
std::streambuf* Logger::cerrBuffer = nullptr; 