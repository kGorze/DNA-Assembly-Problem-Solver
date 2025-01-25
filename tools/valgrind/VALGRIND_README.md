# Valgrind Analysis Setup

This directory contains scripts and configuration files for running Valgrind memory analysis on the project.

## Prerequisites

1. Install Valgrind in WSL (Windows Subsystem for Linux):
```bash
sudo apt-get update
sudo apt-get install valgrind
```

2. Make sure you have CMake installed in WSL:
```bash
sudo apt-get install cmake
```

## Files

- `run_valgrind.sh` - Main script for running Valgrind analysis on different modes of the program
- `build_valgrind.cmake` - CMake configuration file with debug and memory checking options
- `configure_valgrind_build.ps1` - PowerShell script to configure and build the project for Valgrind

## Usage

1. First, configure and build the project for Valgrind analysis:
```powershell
.\configure_valgrind_build.ps1
```

2. Switch to WSL and navigate to the project directory:
```bash
cd /mnt/c/Users/konrad_guest/CLionProjects/optymalizacja_kombinatoryczna
```

3. Run the Valgrind analysis:
```bash
./run_valgrind.sh
```

## Output

The Valgrind analysis results will be stored in the `valgrind_logs` directory:
- `valgrind_generate_instance.log` - Results for generate_instance mode
- `valgrind_test_instance.log` - Results for test_instance mode
- `valgrind_debug.log` - Results for debug mode
- `valgrind_tuning.log` - Results for tuning mode
- `valgrind_tuning_hybrid.log` - Results for tuning_hybrid mode
- `combined_valgrind.log` - Combined results from all modes

## Modes

The Valgrind analysis will be run on all available modes:

1. Generate Instance:
```
generate_instance -n 400 -k 8 -dk 1 -ln 10 -lp 10 -o generated_instance.txt
```

2. Test Instance:
```
test_instance -i generated_instance.txt -o test_results.txt -pid 12345 -diff Medium
```

3. Debug Mode:
```
debug --debug
```

4. Tuning:
```
tuning -cfg config.cfg -out tuning_results.csv
```

5. Tuning Hybrid:
```
tuning_hybrid -cfg config.cfg -out hybrid_tuning_results.csv
```

## Customization

You can modify the parameters for each mode in the `run_valgrind.sh` script. The script includes a `run_valgrind` function that takes three parameters:
1. Mode name
2. Log file path
3. Command line arguments

## Troubleshooting

1. If you get permission errors when running the scripts:
   ```bash
   chmod +x run_valgrind.sh
   ```

2. If Valgrind reports "Command not found":
   - Make sure Valgrind is installed in WSL
   - Make sure you're running the script from WSL, not Windows PowerShell

3. If the build fails:
   - Check CMake configuration
   - Make sure all dependencies are installed
   - Check compiler errors in the build output 