# FastMDXplora PowerShell dynamic banner wrapper
# Version v8: ASCII-safe, generic, no mini box, compact reporting box at end
# Usage preview only:
#   powershell -ExecutionPolicy Bypass -File .\fastmdxplora_run_with_banner.ps1
# Usage run with banner:
#   powershell -ExecutionPolicy Bypass -File .\fastmdxplora_run_with_banner.ps1 --system 1L2Y --output test_run ...

$FastMDXArgs = @($args)

# -----------------------------
# Defaults shown in banner
# -----------------------------
$FMDX_SYSTEM = "ANY_PDB_ID"
$FMDX_OUTPUT = "output_folder"
$FMDX_INCLUDE = "setup -> simulation -> analysis -> report"
$FMDX_FORCEFIELD = "charmm36"
$FMDX_PH = "7.0"
$FMDX_ION_CONC = "0.15"
$FMDX_NVT_STEPS = "default"
$FMDX_NPT_STEPS = "default"
$FMDX_PROD_STEPS = "default"
$FMDX_TIMESTEP = "2.0"
$FMDX_TEMP = "300"
$FMDX_PLATFORM = "auto"
$FMDX_PRECISION = "mixed"
$FMDX_TRAJ_INTERVAL = "default"
$FMDX_ANALYSES = "rmsd -> rg"
$FMDX_REPORT_TITLE = "FastMDXplora Run"

# -----------------------------
# Helpers
# -----------------------------
function Join-UntilNextFlag {
    param(
        [string[]]$Items,
        [int]$StartIndex,
        [string]$Joiner
    )
    $vals = New-Object System.Collections.Generic.List[string]
    $i = $StartIndex
    while ($i -lt $Items.Count -and -not $Items[$i].StartsWith("--")) {
        $vals.Add($Items[$i])
        $i++
    }
    return @($vals -join $Joiner, $i)
}

function Shorten-Text {
    param([string]$Text, [int]$Width)
    if ($null -eq $Text) { $Text = "" }
    if ($Text.Length -le $Width) { return $Text }
    if ($Width -le 3) { return $Text.Substring(0, $Width) }
    return $Text.Substring(0, $Width - 3) + "..."
}

function Pad-RightSafe {
    param([string]$Text, [int]$Width)
    $Text = Shorten-Text $Text $Width
    return $Text + (" " * [Math]::Max(0, $Width - $Text.Length))
}

function Write-BoxTop {
    param([int]$Width)
    Write-Host ("+" + ("-" * $Width) + "+") -ForegroundColor Cyan
}

function Write-BoxLine {
    param([string]$Text, [int]$Width = 95)
    Write-Host -NoNewline "|" -ForegroundColor Cyan
    Write-Host -NoNewline (Pad-RightSafe $Text $Width) -ForegroundColor White
    Write-Host "|" -ForegroundColor Cyan
}

function Write-BoxBlank {
    param([int]$Width = 95)
    Write-BoxLine "" $Width
}

function Write-ReportingBox {
    $reportWidth = 52
    Write-Host ("+" + ("-" * $reportWidth) + "+") -ForegroundColor Blue
    Write-Host -NoNewline "|" -ForegroundColor Blue
    Write-Host -NoNewline (Pad-RightSafe " REPORTING & OUTPUTS" $reportWidth) -ForegroundColor White
    Write-Host "|" -ForegroundColor Blue
    Write-Host -NoNewline "|" -ForegroundColor Blue
    Write-Host -NoNewline (Pad-RightSafe " Markdown reports  |  HTML summaries  |  PDF figures" $reportWidth) -ForegroundColor White
    Write-Host "|" -ForegroundColor Blue
    Write-Host -NoNewline "|" -ForegroundColor Blue
    Write-Host -NoNewline (Pad-RightSafe " PowerPoint slides |  PNG/SVG plots   |  ZIP result bundles" $reportWidth) -ForegroundColor White
    Write-Host "|" -ForegroundColor Blue
    Write-Host ("+" + ("-" * $reportWidth) + "+") -ForegroundColor Blue
    Write-Host ""
    Write-Host -NoNewline "[OK] " -ForegroundColor Green
    Write-Host -NoNewline "Reproducible  " -ForegroundColor White
    Write-Host -NoNewline "[OK] " -ForegroundColor Green
    Write-Host -NoNewline "Modular  " -ForegroundColor White
    Write-Host -NoNewline "[OK] " -ForegroundColor Green
    Write-Host -NoNewline "Publication-ready  " -ForegroundColor White
    Write-Host -NoNewline "[OK] " -ForegroundColor Green
    Write-Host "Open and extensible" -ForegroundColor White
}

# -----------------------------
# Parse fastmdx explore arguments
# -----------------------------
$i = 0
while ($i -lt $FastMDXArgs.Count) {
    switch ($FastMDXArgs[$i]) {
        "--system" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_SYSTEM = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "-s" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_SYSTEM = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--output" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_OUTPUT = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--setup-forcefield" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_FORCEFIELD = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--setup-ph" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_PH = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--setup-ion-concentration-M" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_ION_CONC = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--simulate-nvt-steps" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_NVT_STEPS = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--simulate-npt-steps" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_NPT_STEPS = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--simulate-production-steps" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_PROD_STEPS = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--simulate-timestep-fs" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_TIMESTEP = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--simulate-temperature-K" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_TEMP = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--simulate-platform" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_PLATFORM = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--simulate-precision" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_PRECISION = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--simulate-trajectory-interval-steps" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_TRAJ_INTERVAL = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--report-title" { if ($i + 1 -lt $FastMDXArgs.Count) { $FMDX_REPORT_TITLE = $FastMDXArgs[$i + 1] }; $i += 2; continue }
        "--include" {
            $joined = Join-UntilNextFlag $FastMDXArgs ($i + 1) " -> "
            if ($joined[0]) { $FMDX_INCLUDE = $joined[0] }
            $i = [int]$joined[1]
            continue
        }
        "--analyze-analyses" {
            $joined = Join-UntilNextFlag $FastMDXArgs ($i + 1) " -> "
            if ($joined[0]) { $FMDX_ANALYSES = $joined[0] }
            $i = [int]$joined[1]
            continue
        }
        default { $i++ }
    }
}

# -----------------------------
# Print banner
# -----------------------------
$width = 95
Write-Host ""
Write-Host "FastMDXplora" -ForegroundColor White
Write-Host "Fully Automated SysTem for Molecular Dynamics eXploration" -ForegroundColor Cyan
Write-Host ""

Write-BoxTop $width
Write-BoxLine " FASTMDXPLORA CORE" $width
Write-BoxBlank $width
Write-BoxLine (" System:      {0}" -f $FMDX_SYSTEM) $width
Write-BoxLine (" Output:      {0}" -f $FMDX_OUTPUT) $width
Write-BoxLine (" Pipeline:    {0}" -f $FMDX_INCLUDE) $width
Write-BoxLine (" Platform:    {0}    Precision: {1}" -f $FMDX_PLATFORM, $FMDX_PRECISION) $width
Write-BoxBlank $width
Write-BoxLine (" Setup:       pH {0}    Ion concentration: {1} M    Force field: {2}" -f $FMDX_PH, $FMDX_ION_CONC, $FMDX_FORCEFIELD) $width
Write-BoxLine (" Simulation:  NVT {0} steps    NPT {1} steps    Production {2} steps" -f $FMDX_NVT_STEPS, $FMDX_NPT_STEPS, $FMDX_PROD_STEPS) $width
Write-BoxLine (" Runtime:     timestep {0} fs    temperature {1} K    trajectory every {2} steps" -f $FMDX_TIMESTEP, $FMDX_TEMP, $FMDX_TRAJ_INTERVAL) $width
Write-BoxLine (" Analysis:    {0}" -f $FMDX_ANALYSES) $width
Write-BoxLine (" Report:      {0}" -f $FMDX_REPORT_TITLE) $width
Write-BoxTop $width
Write-Host ""

Write-ReportingBox
Write-Host ""

# -----------------------------
# Run FastMDXplora only when arguments are provided
# -----------------------------
if ($FastMDXArgs.Count -gt 0) {
    Write-Host "Starting FastMDXplora..." -ForegroundColor Green
    & fastmdx explore @FastMDXArgs
    exit $LASTEXITCODE
} else {
    Write-Host "Preview only. Add fastmdx explore flags after this script to run a study." -ForegroundColor DarkGray
}
