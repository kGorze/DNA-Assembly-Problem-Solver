#!/bin/bash

# Exit on any error
set -e

# Create build directory if it doesn't exist
mkdir -p build_wsl

# Clean build directory
rm -rf build_wsl/*

# Go to build directory
cd build_wsl

# Configure with debug symbols and no optimization for better valgrind output
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-g -O0"

# Build
cmake --build . -- -j$(nproc)

# Run tests with valgrind
echo "Running valgrind..."
valgrind --leak-check=full \
         --show-leak-kinds=all \
         --track-origins=yes \
         --verbose \
         --log-file=valgrind-out.txt \
         ./optymalizacja_kombinatoryczna

# Print valgrind output
echo "Valgrind output:"
cat valgrind-out.txt 