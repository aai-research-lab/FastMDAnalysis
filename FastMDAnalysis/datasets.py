"""
Datasets Module for FastMDAnalysis

This module provides dataset classes for example MD systems.
For each dataset (e.g., Ubiquitin, Trp-cage) it specifies the trajectory and topology file paths,
as well as additional attributes describing simulation conditions such as:
  - time_step: simulation time step (in picoseconds)
  - force_field: the force field used
  - integrator: the integration algorithm used
  - temperature: the simulation temperature (in Kelvin)
  - pressure: the simulation pressure (in atm or bar as defined)
  - rtc: a string representing the run-time, creation date, or other relevant tag

The data directory is assumed to be located in the parent of the FastMDAnalysis package directory.
"""

from pathlib import Path

# Set the data directory to be the grandparent of this file joined with "data".
DATA_DIR = Path(__file__).resolve().parent.parent / "data"

class Ubiquitin:
    """
    Ubiquitin Dataset

    Attributes
    ----------
    traj : str
        Absolute path to the ubiquitin trajectory file.
    top : str
        Absolute path to the ubiquitin topology file.
    time_step : float
        Time step used in the simulation (picoseconds).
    force_field : str
        The force field used.
    integrator : str
        The integration method used.
    temperature : float
        Simulation temperature in Kelvin.
    pressure : float
        Simulation pressure (e.g., in atm or bar).
    md_engine : str
        molecular dynamics simulation engine.
    """
    traj = str((DATA_DIR / "ubiquitin.dcd").resolve())
    top = str((DATA_DIR / "ubiquitin.pdb").resolve())
    time_step = 0.002
    force_field = "CHARMM36m"
    integrator = "Legenvin"
    temperature = 300
    pressure = 1.0
    md_engine = "Gromacs" 

class TrpCage:
    """
    Trp-cage Dataset

    Attributes
    ----------
    traj : str
        Absolute path to the trp-cage trajectory file.
    top : str
        Absolute path to the trp-cage topology file.
    time_step : float
        Time step used in the simulation (picoseconds).
    force_field : str
        The force field used.
    integrator : str
        The integration algorithm used.
    temperature : float
        Simulation temperature in Kelvin.
    pressure : float
        Simulation pressure.
    md_engine : str
        molecular dynamics simulation engine.
    """
    traj = str((DATA_DIR / "trp_cage.dcd").resolve())
    top = str((DATA_DIR / "trp_cage.pdb").resolve())
    time_step = 0.002
    force_field = "CHARMM36m"
    integrator = "LegenvinMiddleIntegrator"
    temperature = 300
    pressure = 1.0
    md_engine = "OpenMM 8.2"  

# Convenience shortcuts for easy import.
ubiquitin = Ubiquitin
trp_cage = TrpCage

