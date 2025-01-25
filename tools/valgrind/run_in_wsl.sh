#!/bin/bash

# Install required packages if not already installed
echo "Installing required packages..."
sudo apt-get update
sudo apt-get install -y cmake g++ valgrind

# Create build directory
mkdir -p build_valgrind
cd build_valgrind

# Configure with CMake
echo "Configuring CMake for Valgrind build..."
cmake -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_FLAGS="-g -O0 -fno-inline -fno-omit-frame-pointer" \
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
      ..

# Build the project
echo "Building project..."
cmake --build . --config Debug

# Create logs directory
mkdir -p ../tools/valgrind/valgrind_logs

# Run Valgrind
echo "Running Valgrind analysis..."
valgrind --tool=memcheck \
         --leak-check=full \
         --show-leak-kinds=all \
         --track-origins=yes \
         --verbose \
         --log-file=../tools/valgrind/valgrind_logs/valgrind_test.log \
         ./optymalizacja_kombinatoryczna --test_instance

echo "Valgrind analysis completed. Results saved to tools/valgrind/valgrind_logs/valgrind_test.log" 