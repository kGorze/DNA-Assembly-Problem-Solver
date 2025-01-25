# CMake configuration file for Valgrind build

# Set build type to Debug for better Valgrind analysis
set(CMAKE_BUILD_TYPE Debug)

# Add debug flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

# Add flags for better memory checking
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")

# Enable address sanitizer if available
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
endif()

# Disable optimizations that could interfere with Valgrind
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-inline")

# Add debug symbols
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gdwarf-4") 