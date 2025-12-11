param(
    [string]$BuildType = "Debug",
    [string]$Compiler = "g++",
    [switch]$Tests
)

Write-Host "🔨 Starting CMake build ($BuildType)..."

$current_dir=$pwd
$project_dir=$pwd
if ($PSScriptRoot -eq $pwd){
    $project_dir = Split-Path -parent $PSScriptRoot
}

if (-not (Test-Path "$project_dir/build")) {
    Write-Host "📁 Creating build directory"
    New-Item -ItemType Directory -Path "$project_dir/build" | Out-Null
}

Set-Location "$project_dir/build"

# Determine C++ compiler (use environment detection)
$cxx_compiler = switch ($Compiler) {
    "clang" { "clang++" }
    default { "g++" }
}

# Check if the chosen compiler exists in PATH
if (-not (Get-Command $cxx_compiler -ErrorAction SilentlyContinue)) {
    Write-Host "❌ Compiler '$cxx_compiler' not found in PATH." -ForegroundColor Red
    exit 1
}


$cmakeArgs = @("..", "-DCMAKE_BUILD_TYPE=$BuildType", "-DCMAKE_CXX_COMPILER=$cxx_compiler")

# Enable tests if requested
$cmakeArgs += "-DPSO-CPP_BUILD_TESTS=$((($Tests.IsPresent -replace $true,'ON') -replace $false,'OFF'))"
# Configure
Write-Host "⚙️  Configuring project..."
cmake @cmakeArgs

# Build
Write-Host "🚧 Building..."
cmake --build . --config $BuildType

Set-Location $current_dir
Write-Host "✅ Build finished!"
