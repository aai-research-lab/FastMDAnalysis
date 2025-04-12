from FastMDAnalysis import RMSDAnalysis, load_trajectory

traj = load_trajectory("traj.dcd", "top.pdb")
analyzer = RMSDAnalysis(traj, ref_frame=0)
results = analyzer.run()
analyzer.save_results()
analyzer.plot()
