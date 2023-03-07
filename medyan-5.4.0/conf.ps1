# System requirements for this script:
#   - git >= 2.7.0
#   - cmake >= 3.12.0
#   - Latest C++ compiler

$medyan_root_dir = $PSScriptRoot
$medyan_conf_current_dir = $(Get-Location)

try {
    Write-Host "Bootstrapping..."
    & "$medyan_root_dir\scripts\bootstrap.ps1"
}
finally {
    Set-Location $medyan_conf_current_dir
}
