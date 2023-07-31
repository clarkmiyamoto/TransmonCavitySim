import pandas as pd

from qiskit_metal.analyses.quantization import LOManalysis

class AnsysQ3D:

  def __init__(self, design):
    self.design = design
    self._define_setup(design)

  def _define_setup(self):
    """
    Adjust Ansys -> Setup settings here.
    """ 
    c1 = LOManalysis(self.design, "q3d")
    c1.sim.setup.reuse_selected_design = False
    c1.sim.setup.reuse_setup = False
    c1.sim.setup.max_passes = 30
    c1.sim.setup.min_converged_passes = 1
    c1.sim.setup.percent_error = 0.1

    self.c1 = c1

  def extract_CapMatrix(self, **kwargs) -> pd.DataFrame:
    """
    Runs Q3D simulation for capacitance matrix.

    Args:
      kwargs (dict): Goes into Qiskit Metal's LOManalysis.sim.run(...).
    
    Returns:
      CapMatrix (pd.DataFrame): Capacitance Matrix
    """
    self.c1.sim.run(**kwargs)
    CapMatrix = c1.sim.capacitance_matrix

    return CapMatrix
