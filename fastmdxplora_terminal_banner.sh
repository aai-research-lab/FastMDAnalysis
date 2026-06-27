#!/usr/bin/env bash
# FastMDXplora terminal welcome banner - core-only version
# No ASCII logo diagram; shows core simulation startup details and outputs.

BOLD='\033[1m'
RESET='\033[0m'
CYAN='\033[38;5;45m'
GREEN='\033[38;5;40m'
BLUE='\033[38;5;33m'
ORANGE='\033[38;5;208m'
WHITE='\033[38;5;255m'
GRAY='\033[38;5;245m'

printf "\n${WHITE}${BOLD}FastMDXplora${RESET}\n"
printf "${CYAN}${BOLD}Fully Automated SysTem for Molecular Dynamics eXploration${RESET}\n\n"

printf "${CYAN}┌──────────────────────────────────────────────────────────────────────────────────────────────┐${RESET}\n"
printf "${CYAN}│${WHITE}${BOLD} FASTMDXPLORA CORE                                                                           ${CYAN}│${RESET}\n"
printf "${CYAN}│${RESET} Single-command molecular dynamics workflow for setup, simulation, analysis, and reporting.        ${CYAN}│${RESET}\n"
printf "${CYAN}│${RESET}                                                                                              ${CYAN}│${RESET}\n"
printf "${CYAN}│${RESET} ${CYAN}${BOLD}Trigger:${RESET}    Runs when a simulation phase is requested from CLI or config.                ${CYAN}│${RESET}\n"
printf "${CYAN}│${RESET} ${CYAN}${BOLD}Input:${RESET}      PDB ID or local PDB/mmCIF; optional ligand/topology; YAML/TOML config.       ${CYAN}│${RESET}\n"
printf "${CYAN}│${RESET} ${CYAN}${BOLD}Setup:${RESET}      pH, ion concentration, force field, solvation, reproducible parameters.     ${CYAN}│${RESET}\n"
printf "${CYAN}│${RESET} ${CYAN}${BOLD}Simulation:${RESET} OpenMM CPU/GPU, precision, NVT/NPT/production, timestep, checkpoints.     ${CYAN}│${RESET}\n"
printf "${CYAN}│${RESET} ${CYAN}${BOLD}Analysis:${RESET}   RMSD, Rg, RMSF, H-bonds, contacts, PCA/clustering, ligand metrics.        ${CYAN}│${RESET}\n"
printf "${CYAN}│${RESET} ${CYAN}${BOLD}Energy:${RESET}     energy.csv -> energy trace, simulation health plot, summary table.       ${CYAN}│${RESET}\n"
printf "${CYAN}│${RESET} ${CYAN}${BOLD}Reports:${RESET}    report.md, slides.pptx, captions, manifest.json, project bundle.          ${CYAN}│${RESET}\n"
printf "${CYAN}└──────────────────────────────────────────────────────────────────────────────────────────────┘${RESET}\n\n"

printf "${ORANGE}┌──────────────────────────────────────────────────────────────────────────────────────────────┐${RESET}\n"
printf "${ORANGE}│${WHITE}${BOLD} SIMULATION COMMAND PATTERN                                                                 ${ORANGE}│${RESET}\n"
printf "${ORANGE}│${GREEN}$${RESET} fastmdx explore --system 1L2Y --output trpcage_100steps                                ${ORANGE}│${RESET}\n"
printf "${ORANGE}│${RESET}   --include setup simulation analysis report --setup-ph 7.4 --setup-ion-concentration-M 0.15     ${ORANGE}│${RESET}\n"
printf "${ORANGE}│${RESET}   --simulate-nvt-steps 100 --simulate-npt-steps 0 --simulate-production-steps 100              ${ORANGE}│${RESET}\n"
printf "${ORANGE}│${RESET}   --simulate-timestep-fs 0.5 --simulate-temperature-K 100 --simulate-platform CPU               ${ORANGE}│${RESET}\n"
printf "${ORANGE}│${RESET}   --simulate-precision double --simulate-trajectory-interval-steps 10                           ${ORANGE}│${RESET}\n"
printf "${ORANGE}│${RESET}   --simulate-checkpoint-interval-steps 0 --analyze-analyses rmsd rg                              ${ORANGE}│${RESET}\n"
printf "${ORANGE}│${RESET}   --report-title \"Trp-cage 100 Step Test\"                                                   ${ORANGE}│${RESET}\n"
printf "${ORANGE}└──────────────────────────────────────────────────────────────────────────────────────────────┘${RESET}\n\n"

printf "${BLUE}┌──────────────────────────────────────────────────────────────────────────────────────────────┐${RESET}\n"
printf "${BLUE}│${WHITE}${BOLD} REPORTING & OUTPUTS                                                                        ${BLUE}│${RESET}\n"
printf "${BLUE}│${RESET} Reports:      Markdown report | PowerPoint slides | project_bundle.zip                         ${BLUE}│${RESET}\n"
printf "${BLUE}│${RESET} Figures:      RMSD/Rg plots | energy_trace.png | simulation_health.png                        ${BLUE}│${RESET}\n"
printf "${BLUE}│${RESET} Data:         energy_summary.csv | trajectory + logs | analysis_manifest.json                   ${BLUE}│${RESET}\n"
printf "${BLUE}│${RESET} Review-ready: figure captions | reproducibility context | README/GitHub graphic                  ${BLUE}│${RESET}\n"
printf "${BLUE}└──────────────────────────────────────────────────────────────────────────────────────────────┘${RESET}\n\n"

printf "${GREEN}${BOLD}✓${RESET} Reproducible  ${GREEN}${BOLD}✓${RESET} Config-driven  ${GREEN}${BOLD}✓${RESET} OpenMM CPU/GPU  ${GREEN}${BOLD}✓${RESET} Energy-aware  ${GREEN}${BOLD}✓${RESET} Publication-ready\n"
