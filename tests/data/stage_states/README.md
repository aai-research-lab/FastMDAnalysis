These two files were written by OpenMM 8.6.1, using its Reference platform.
A Context was created but no integration or minimization was run. The
four-atom GLY Topology was written with app.PDBFile.writeFile(keepIds=True).
The State was written with XmlSerializer.serialize after setPositions and
getState(getPositions=True, enforcePeriodicBox=True).

Topology positions in nm: (0.1,0.2,0.3), (0.24,0.2,0.3), (0.37,0.2,0.3),
(0.47,0.25,0.3). State positions in nm: (0.3,0.4,0.5), (0.44,0.4,0.5),
(0.57,0.4,0.5), (0.67,0.45,0.5). Both use periodic box vectors in nm:
(3,0,0), (0.4,2.8,0), (0.2,0.3,3.2). Different coordinates ensure that
reading the topology instead of the State cannot pass the coordinate test.
