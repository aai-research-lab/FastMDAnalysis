These audit records and PDB snapshots were written by PreparationRecorder
from this repository, with OpenMM 8.6.1 PDBFile and Topology. No chemical
preparation, minimization, integration or protonation prediction ran.

The four-atom GLY topology from stage_states was recorded, then one H atom
was explicitly appended to the Topology at (0.07, 0.2, 0.3) nm and bonded
to N. Recording this deliberate change exercises the reader's atom
correspondence; it is not evidence that a chemistry backend chose a state.
The recorded pH is a requested setting, not an independently verified
chemical protonation state.
