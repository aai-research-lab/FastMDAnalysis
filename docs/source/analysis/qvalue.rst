.. currentmodule:: FastMDAnalysis.analysis.qvalue

Fraction of Native Contacts (Q-Value)
=====================================

The :class:`QAnalysis` module computes the fraction of native contacts (Q-value)
according to the Best-Hummer-Eaton definition. This metric measures how many
native contacts (heavy atom pairs present in the reference structure) are
maintained in each frame of the trajectory, providing a scalar descriptor of
protein folding state and native structure retention.

Theory
------

The fraction of native contacts is defined as:

.. math::

   Q(X) = \frac{1}{|S|} \sum_{(i,j) \in S} \frac{1}{1 + \exp\left[\beta(r_{ij}(X) - \lambda r^0_{ij})\right]}

where:

- :math:`X` is a conformation (trajectory frame)
- :math:`r_{ij}(X)` is the distance between heavy atoms :math:`i` and :math:`j` in frame :math:`X`
- :math:`r^0_{ij}` is the reference distance (native state)
- :math:`S` is the set of native contact pairs: heavy atom pairs >3 residues apart with :math:`r^0_{ij} < 0.45` nm
- :math:`\beta = 50` nm⁻¹ (steepness parameter)
- :math:`\lambda = 1.8` (distance scaling factor)

Q-values range from 0 (native contacts lost) to 1 (all native contacts maintained).

Reference
---------

Best, R. B., Hummer, G., & Eaton, W. A. (2013). Native contacts determine 
protein folding mechanisms in atomistic simulations. *Proceedings of the 
National Academy of Sciences*, 110(44), 17874-17879. 
https://doi.org/10.1073/pnas.1311599110

Implementation Highlights
-------------------------

- Automatically identifies native contacts from a reference frame (default: frame 0)
- Native contacts are heavy atom pairs >3 residues apart within 0.45 nm in reference
- Applies the sigmoid formula for each native contact pair
- Reports native contact count and stores all parameters in metadata
- Includes metadata annotation on plot for reproducibility
- Fully configurable parameters with sensible defaults

Usage
-----

**API**

.. code-block:: python

   from fastmdanalysis import FastMDAnalysis

   fastmda = FastMDAnalysis("traj.dcd", "top.pdb")
   q = fastmda.qvalue(
       reference_frame=0,
       beta_const=50.0,      # nm^-1
       lambda_const=1.8,
       native_cutoff=0.45    # nm
   )
   q_values = q.results["qvalue"]
   q.plot()

**CLI**

.. code-block:: bash

   fastmda qvalue -traj traj.dcd -top top.pdb \
       --reference-frame 0 \
       --beta 50.0 \
       --lambda 1.8 \
       --cutoff 0.45 \
       -o q_output

**Multi-Analysis**

.. code-block:: python

   # Include Q-value in multi-analysis workflow
   results = fastmda.analyze(
       include=["rmsd", "rg", "qvalue"],
       output="analysis_output"
   )

Parameters
----------

``reference_frame`` : int, default=0
    Frame index to use as the native/reference state.
    Aliases: ``ref``, ``reference``

``beta_const`` : float, default=50.0
    Beta constant in nm⁻¹. Controls the steepness of the sigmoid function.
    Alias: ``--beta``

``lambda_const`` : float, default=1.8
    Lambda constant (dimensionless). Scaling factor for native distances.
    Alias: ``--lambda``

``native_cutoff`` : float, default=0.45
    Cutoff distance in nm for identifying native contacts.
    Only heavy atom pairs within this distance in the reference frame are considered.
    Alias: ``--cutoff``

``atoms`` : str, optional
    MDTraj atom selection string for restricting analysis to a subset of atoms.
    If not provided, heavy atoms are used automatically.
    Aliases: ``atom_indices``, ``selection``

Outputs
-------

- ``qvalue.dat`` — Q-values (one per frame, range [0, 1])
- ``qvalue_stats.dat`` — mean and standard deviation when ``compute_stat=True``
- ``qvalue.png`` — Line plot of Q-value vs. frame with metadata annotation
- ``qvalue_metadata.json`` — Complete metadata including:
  
  - ``native_contacts_count`` — Number of native contacts identified
  - ``reference_frame`` — Reference frame index used
  - ``beta_const_nm_inv`` — Beta parameter value
  - ``lambda_const`` — Lambda parameter value
  - ``native_cutoff_nm`` — Cutoff distance used
  - ``n_frames`` — Total frames analyzed
  - ``n_atoms`` — Total atoms in trajectory

Interpretation
--------------

- **Q ≈ 1.0**: Native structure well-preserved; strong similarity to reference
- **Q ≈ 0.5**: Intermediate state; half of native contacts maintained
- **Q ≈ 0.0**: Native structure significantly disrupted; few contacts retained

The trajectory slice showing how Q evolves provides insights into:

- Protein stability and folding kinetics
- Unfolding or refolding pathways
- Conformational dynamics relative to native state
