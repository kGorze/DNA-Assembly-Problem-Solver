# PowerShell script to build and run Valgrind analysis through WSL

# First, configure and build the project
Write-Host "Configuring and building project for Valgrind..."
& "$PSScriptRoot\configure_valgrind_build.ps1"

if ($LASTEXITCODE -ne 0) {
    Write-Host "Error: Build failed"
    exit 1
}

# Get the current directory path and convert it to WSL format
$currentPath = (Get-Location).Path
$wslPath = $currentPath -replace '^([A-Za-z]):', '/mnt/$($1.ToLower())' -replace '\\', '/'

Write-Host "Running Valgrind analysis through WSL..."
Write-Host "Project path in WSL: $wslPath"

# Create the valgrind_logs directory if it doesn't exist
$logsDir = "$PSScriptRoot/valgrind_logs"
if (-not (Test-Path $logsDir)) {
    New-Item -ItemType Directory -Path $logsDir | Out-Null
}

# Run Valgrind through WSL
$valgrindCmd = "cd $wslPath && valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes ./build_valgrind/optymalizacja_kombinatoryczna --test_instance 2>&1 | tee $logsDir/valgrind_test.log"
Write-Host "Executing: $valgrindCmd"

wsl bash -c "$valgrindCmd"

if ($LASTEXITCODE -eq 0) {
    Write-Host "Valgrind analysis completed successfully"
    Write-Host "Results saved to: $logsDir/valgrind_test.log"
} else {
    Write-Host "Error: Valgrind analysis failed"
    exit 1
} 