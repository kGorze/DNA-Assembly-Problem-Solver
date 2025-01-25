# PowerShell script to configure and build for Valgrind

# Create build directory if it doesn't exist
if (-not (Test-Path "build_valgrind")) {
    New-Item -ItemType Directory -Path "build_valgrind"
}

# Navigate to build directory
Set-Location -Path "build_valgrind"

# Configure with CMake
Write-Host "Configuring CMake for Valgrind build..."
& "C:\Program Files\JetBrains\CLion 2023.3.4\bin\cmake\win\x64\bin\cmake.exe" -DCMAKE_BUILD_TYPE=Debug `
      -DCMAKE_CXX_FLAGS="-g -O0 -fno-inline -fno-omit-frame-pointer" `
      -DCMAKE_EXPORT_COMPILE_COMMANDS=ON `
      ..

# Build the project
Write-Host "Building project..."
& "C:\Program Files\JetBrains\CLion 2023.3.4\bin\cmake\win\x64\bin\cmake.exe" --build . --config Debug

# Return to original directory
Set-Location -Path ".."

Write-Host "Build completed. You can now run the Valgrind analysis using run_valgrind.sh" 