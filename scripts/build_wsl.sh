#!/bin/bash

# Exit on error
set -e

# Create build directory if it doesn't exist
mkdir -p build_wsl

# Navigate to build directory
cd build_wsl

# Configure with CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
cmake --build . -j$(nproc)

# Run tests
ctest --output-on-failure

echo "Build and tests completed successfully!" 