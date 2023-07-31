import pandas as pd
from datetime import datetime

from qiskit_metal.analyses.quantization import LOManalysis
from qiskit_metal.analyses.quantization import EPRanalysis
import pyEPR as epr


### Global variables in `ansys.py`
current_date = datetime.now()
current_date_only = current_date.strftime("%Y%m%d")

project_name = f'CavityQubit_proj_{current_date_only}'
eigenmodal_design_name = f'CavityQubit_eigenmodal_{current_date_only}'
q3d_design_name = f'CavityQubit_q3d_{current_date_only}'

### Ansys Renderers
# (1/2) Eigenmodal Renderer
class AnsysEigenmodal:
  """Ansys Eigenmodal Renderer in Qiskit Metal"""

  def __init__(self, design: "MultiPlanar"):
    self.design = design
    self.renderer = self._define_renderer(design)

  def _define_renderer(self) -> EPRanalysis:
    renderer = EPRanalysis(design, "hfss")
    return renderer

  def _define_setup(self, qubit_name) -> None:
    """
    Updates `self.renderer`
    """
    # For Ansys -> Setup
    renderer = self.renderer
    renderer.sim.setup.name = 'Setup'
    renderer.sim.setup.min_freq_ghz = 1
    renderer.sim.setup.n_modes = 2
    renderer.sim.setup.max_passes = 23
    renderer.sim.setup.max_delta_f = 0.1
    renderer.sim.setup.min_converged = 1

    # For pyEPR
    Lj = self.design.components[qubit_name].options.hfss_inductance
    Cj = self.design.components[qubit_name].options.hfss_capacitance    
    if (type(Lj) != str) or (type(Cj) != str):
      raise ValueError(f"In `self.design.components[{qubit_name}].options`, `hfss_inductance` and `hfss_capacitance` must be strings w/ units."
    
    renderer.sim.setup.vars = Dict({'Lj': Lj, 'Cj': Cj})
    renderer.setup.junctions.jj.rect = f'JJ_rect_Lj_{qubit_name}_rect_jj'
    renderer.setup.junctions.jj.line = 'JJ_Lj_{qubit_name}_rect_jj_'

    self.renderer = renderer

  def run(self, 
          qubit_name: str,
          cpw_name: str,
          connector_pad_name: str,
          open_pin: tuple[str, str],
          pyEPR_options: dict = None):
    """Main functionality, run eigenmode simulation & extract pyEPR results."""
    self._define_setup(qubit_name)
    self._create_design_and_preseed_mesh(qubit_name, cpw_name, connector_pad_name, open_pin)

    # Run ANSYS Simulation
    self.renderer.sim.renderer.pinfo.setup.analyze()

    # Run EPR analysis
    self.renderer.run_epr(**pyEPR_options)

    # Release Ansys session
    self.renderer.sim.close()

  def _create_design_and_preseed_mesh(self, 
                                      qubit_name: str, 
                                      cpw_name: str,
                                      connector_pad_name: str,
                                      open_pin: tuple[str, str]):
    """Automatically render design to Ansys and preseed Mesh"""
    hfss = self.renderer.sim.renderer

    ### Renderer design
    hfss.start()
    hfss.activate_ansys_design(eigenmodal_design_name, 'eigenmode')                                    
    hfss.clean_active_design()

    hfss._render(
      selection=[qubit_name, cpw_name],
      solution_type="eigenmode",
      vars_to_initialize=self.renderer.setup.vars,
      open_pins=[open_pin]
    )

    ### Pre-Seeding Mesh
    mesh_sigfigs = 3
    divide_characteristic_length = 3

    ## CPW Mesh
    # Parsing `cpw_width`
    cpw_width = self.design.components[cpw_name].options.cpw_width
    if (type(cpw_width) == str):
      cpw_width_value = float(re.search(r'\d+\.?\d*', cpw_width).group())
      cpw_width_unit = re.search(r'[a-zA-Z]+', cpw_width).group()
    elif (type(cpw_width) == float) or (type(cpw_width) == int):
      cpw_width_value = cpw_width
      cpw_width_unit = "mm"
    else:
      raise ValueError(f"Couldn't find a `cpw_width` from inputted cpw_name: `{cpw_name}`")

    # Adding CPW mesh
    cpw_MaxLength = str(round(cpw_width_value / divide_characteristic_length, mesh_sigfigs)) + cpw_width_unit
    hfss.modeler.mesh_length(
                    'cpw_mesh',
                    [f'trace_{cpw_name}', f'{connector_pad_name}_connector_arm_{qubit_name}'],
                    MaxLength=cpw_MaxLength)

    ## Qubit Mesh
    # Parsing `cross_width`
    cross_width = self.design.components[qubit_name].options.cross_width
    if (type(cpw_width) == str):
      cross_width_value = float(re.search(r'\d+\.?\d*', cross_width).group())
      cross_width_unit = re.search(r'[a-zA-Z]+', cross_width).group()
    elif (type(cpw_width) == float) or (type(cpw_width) == int):
      cross_width_value = cross_width
      cross_width_unit = "mm"
    else:
      raise ValueError(f"Couldn't find a `cpw_width` from inputted cpw_name: `{cpw_name}`")

    # Adding Qubit mesh
    qubit_MaxLength = str(round(cross_width_value / divide_characteristic_length, mesh_sigfigs)) + cross_width_unit
    hfss.modeler.mesh_length(
                    'qubit_mesh',
                    [f'cross_{qubit_name}'],
                    MaxLength=qubit_MaxLength)


# (2/2) Q3D Renderer
class AnsysQ3D:
  """Ansys Q3D Renderer in Qiskit Metal"""

  def __init__(self, 
               design: "MultiPlanar", 
               qubit_name: str, 
               connector_arm_name: str, 
               cavity_name: str):
    self.design = design
    self.renderer = self._define_renderer(design)
    
  def _define_renderer(self) -> LOManalysis:
    """Defines Qiskit Metal renderer"""
    renderer = LOManalysis(self.design, "q3d")
    renderer = self._define_setup(renderer)

    return renderer

  def _define_setup(self):
    """Adjust Ansys -> Setup settings here."""
    renderer = self.renderer
    
    renderer.sim.setup.reuse_selected_design = False
    renderer.sim.setup.reuse_setup = False
    renderer.sim.setup.max_passes = 30
    renderer.sim.setup.min_converged_passes = 1
    renderer.sim.setup.percent_error = 0.1

    self.renderer

  def run(self, 
          qubit_name: str,
          cpw_name: str,
          connector_pad_name: str,
          open_pin: tuple[str, str],) -> pd.DataFrame:
    """
    Runs Q3D simulation for capacitance matrix.
    
    Returns:
      CapMatrix (pd.DataFrame): Capacitance Matrix
    """
    self.renderer.run()
    CapMatrix = self.renderer.sim.capacitance_matrix

    return CapMatrix




