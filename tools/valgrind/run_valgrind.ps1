# PowerShell script to run Valgrind analysis in WSL

# Get the current directory path and convert it to WSL format
$currentPath = (Get-Location).Path
$driveLetter = $currentPath[0].ToString().ToLower()
$wslPath = "/mnt/$driveLetter" + $currentPath.Substring(2).Replace('\', '/')

Write-Host "Running Valgrind analysis in WSL..."
Write-Host "Project path in WSL: $wslPath"

# Make sure the shell script is executable
$wslCmd = "chmod +x $wslPath/tools/valgrind/run_in_wsl.sh"
Write-Host "Making script executable: $wslCmd"
wsl bash -c $wslCmd

# Run the shell script in WSL
$wslCmd = "cd $wslPath && ./tools/valgrind/run_in_wsl.sh"
Write-Host "Running Valgrind: $wslCmd"
wsl bash -c $wslCmd

if ($LASTEXITCODE -eq 0) {
    Write-Host "Valgrind analysis completed successfully"
    Write-Host "Results saved to: tools/valgrind/valgrind_logs/valgrind_test.log"
} else {
    Write-Host "Error: Valgrind analysis failed"
    exit 1
} 