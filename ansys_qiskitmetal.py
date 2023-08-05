from sim_template import AbstractSim

from qiskit_metal import Dict
from qiskit_metal.analyses.quantization import EPRanalysis
from qiskit_metal.analyses.quantization import LOManalysis
import pyEPR as epr

from pyaedt import Hfss, Q3d

from datetime import datetime


class AnsysQiskitMetal(AbstractSim):

    def __init__(self, 
                 design: "Planar Design", 
                 qubit_name: str,
                 connection_pad_name: str,
                 cpws_names: list[str],
                 feedline_name: str,
                 other_names: list[str] = [],
                 open_pins: list[tuple[str, str]] = None):
        super().__init__(design=design, selection=[qubit_name] + cpws_names + [feedline_name] + other_names, open_pins=open_pins)

        ### Renderers
        # Qiskit Metal
        self.hfss_renderer = EPRanalysis(design, "hfss")
        self.q3d_renderer  = LOManalysis(design, "q3d")

        # EPR analysis modules
        self.pinfo = None
        self.eprd = None
        self.epra = None

        ### Naming
        # Names QComponents
        self.qubit_name = qubit_name
        self.connection_pad_name = connection_pad_name
        self.cpws_names = cpws_names
        self.feedline_name = feedline_name
        self.other_names = other_names

        # Name of renderers
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M")

        self.hfss_renderer.default_options.project_name = "QubitCavity" + self.timestamp
        self.q3d_renderer.default_options.project_name = "QubitCavity" + self.timestamp
        self.eigenmode_design_name = "QubitCavity_eigenmode" + self.timestamp
        self.q3d_design_name = "QubitCavity_q3d" + self.timestamp

        ### Simulation status
        self.got_EigenModes = False
        self.got_EPR        = False
        self.got_CapMatrix  = False

        ### Results
        # Raw results
        self.undressed_freqs = None  # Classical frequencies from Eigenmodal sim; in `run_EigenModes`.        
        self.other_data_EPR = None   # Instance of `pyEPR.QuantumAnalysis.data` after running epr analysis; in `run_EPR`.
        self.cap_matrix = None       # Instance of `LOManalysis.sim.capacitance_matrix`; in `run_CapMatrix`.

        # Parsed results from simulaton
        self.qubit_freq        = None # Linear MHZ
        self.cavity_freq       = None # Linear MHZ
        self.anharmonicity     = None # Linear MHZ
        self.dispersive_shift  = None # Linear MHZ
        self.coupling_strength = None # Linear MHZ


    def run_all_simulations(self, 
                            eigenmode_options: dict = None,
                            epr_options: dict = None,
                            capmatrix_options: dict = None):
        
        self.get_EigenModes(**eigenmode_options)
        self.run_EPR(**epr_options)
        self.get_CapMatirx(**capmatrix_options)

        self._parse_all_results()

    def run_EigenModes(self,
                       design_name: str = None,
                       setup_name: str = "Setup",
                       min_freq_ghz: int = 2,
                       n_modes: int = 2,
                       max_delta_f: float = 0.1,
                       max_passes: int = 30,
                       min_passes: int = 1,
                       min_converged: int = 1,
                       pct_refinement: int = 30,
                       basis_order: int = -1,
                       qubit_mesh_MaxLength: str = '6um',
                       cavity_mesh_MaxLength: str = '3um',
                       other_mesh_MaxLength: dict = None) -> dict:
        """
        Runs ANSYS HFSS Eigenmode Simulation and returns mode frequencies. These are classical / undressed frequencies.

        Args:
            design_name (str, optional): The name of the HFSS design to be used for the simulation. Defaults to self.eigenmode_design_name
            setup_name (str, optional): The name of the simulation setup to be used.
            min_freq_ghz (int, optional): The minimum frequency of interest for the simulation, in GHz.
            n_modes (int, optional): The number of electromagnetic modes to compute.
            max_delta_f (float, optional): The maximum frequency error allowed for each mode, in GHz.
            max_passes (int, optional): The maximum number of passes.
            min_passes (int, optional): The minimum number of passes.
            min_converged (int, optional): Number of passes to run after reaching `max_delta_f`.
            pct_refinement (int, optional): Maximum allowed refinement for the meshing per pass. 
                                            For example: 30 means 30% of the entire mesh can change per pass.
            basis_order (int, optional): The basis order for the simulation.
            qubit_mesh_MaxLength (str, optional): The maximum mesh length for the qubit region. Include units in str.
            cavity_mesh_MaxLength (str, optional): The maximum mesh length for the cavity region. Include units in str.
            other_mesh_MaxLength (dict, optional): A dictionary of additional regions and their maximum mesh lengths
                in micrometers. The keys should be the region names, and values should be the maximum mesh length.

        """
        ### Naming of design & setup
        if (design_name == None):
            design_name = self.eigenmode_design_name

        self.eigenmode_setup_name = setup_name

        if (n_modes != 2):
            print("WARNING: This simulation is designed for Qubit + Cavity. It is normal to set `n_modes = 2`.")
            print("PROCED WITH CAUTION.")

        ### Notation
        hfss = self.hfss_renderer.sim.renderer

        ### Start ANSYS, Active Design
        hfss.start()
        hfss.activate_ansys_design(design_name, 'eigenmode')

        ### Render design to ANSYS
        hfss.clean_active_design()
        hfss.render_design(selection=self.selection)

        ### Add mesh
        # Qubit Mesh
        hfss.modeler.mesh_length('qubit',
                                 [f'cross_{self.qubit_name}'],
                                 MaxLength=qubit_mesh_MaxLength)

        # CPW Mesh
        trace_names = [f'trace_{cpw_name}' for cpw_name in self.cpws_names]
        claw_name = [f'{self.connection_pad_name}_connector_arm_{self.qubit_name}']
        hfss.modeler.mesh_length('cpw',
                                trace_names + claw_name,
                                MaxLength=cavity_mesh_MaxLength)

        # Other Meshes
        if type(other_mesh_MaxLength) == dict:
            for name, MaxLength in other_mesh_MaxLength.items():
                hfss.modeler.mesh_length(name,
                                        [name],
                                        MaxLength=MaxLength)

        ### Change silicon
        self._pyAEDT_functionality(solutiontype='Eigenmode')

        ### Add Setup
        hfss.add_eigenmode_setup(name=setup_name,
                                 min_freq_ghz=min_freq_ghz,
                                 n_modes=n_modes,
                                 max_delta_f=max_delta_f,
                                 max_passes=max_passes,
                                 min_passes=min_passes,
                                 min_converged=min_converged,
                                 pct_refinement=pct_refinement,
                                 basis_order=basis_order)

        ### Analyze Setup
        hfss.analyze_setup(setup_name)
        
        ### Release ANSYS Session
        self.hfss_renderer.sim.renderer = hfss
        self.hfss_renderer.sim.close()

        ### Log succsessful HFSS Eigenmode Simulation
        self.got_EigenModes = True


    def run_EPR(self, 
                cos_trunc: int = 8, 
                fock_trunc: int = 15,
                print_result: bool = True) -> epr.QuantumAnalysis:
        """Connect to HFSS Eigenmode, run `pyEPR` analysis on it.
        
        Args:
            cos_trunc (int, optional): Truncate Taylor expansion of cosine terms in Transmon Hamiltonian. Defaults to 8.
            fock_trunc (int, optional): Truncate fock space Hamiltonian for numerical diagonalization. Defaults to 15.
             (bool, optional): Display results? Defaults to True.

        Returns
            self.epra (epr.QuantumAnalysis): Obj containing characteristics of Jaynes Cummings Hamiltonian.
            
        """
        if (self.got_EigenModes != True):
            raise RuntimeError("Must run `run_EigenModes` before calling `run_EPR`.")

        ### Connect EPR to ANSYS
        # Start ANSYS w/ Qiskit Metal, Active Design
        hfss.start()
        hfss.activate_ansys_design(design_name, 'eigenmode')
        self.hfss_renderer.sim.renderer.activate_ansys_design(self.eigenmode_design_name, 'eigenmode')
        hfss.close()

        # Launch pyEPR
        self.pinfo = epr.ProjectInfo()

        ### Tells pyEPR where Johsephson Junctions are located in ANSYS
        self.pinfo.junctions = Dict()
        self.pinfo.junctions[f'jj'] = Dict(rect=f'JJ_rect_Lj_{self.qubit_name}_rect_jj', 
                                              line=f'JJ_Lj_{self.qubit_name}_rect_jj_',
                                              Lj_variable=f'Lj', 
                                              Cj_variable=f'Cj')
        self.pinfo.validate_junction_info()

        # Tells pyEPR which components have dissipative elements
        self.pinfo.dissipative['dielectrics_bulk'] = ['main'] 

        ### Extract Energies
        self.eprd = epr.DistributedAnalysis(self.pinfo)
        
        ℰ_elec = self.eprd.calc_energy_electric()
        ℰ_elec_substrate = self.eprd.calc_energy_electric(None, 'main')
        ℰ_mag = eprd.calc_energy_magnetic()

        if print_result:
            print(f"""
            ℰ_elec_all       = {ℰ_elec}
            ℰ_elec_substrate = {ℰ_elec_substrate}
            EPR of substrate = {ℰ_elec_substrate / ℰ_elec * 100 :.1f}%

            ℰ_mag    = {ℰ_mag}
            """)

        ### Run EPR analysis
        self.eprd.do_EPR_analysis()
        self.epra = epr.QuantumAnalysis(self.eprd.data_filename)
        self.epra.analyze_all_variations(cos_trunc=cos_trunc, 
                                         fock_trunc=fock_trunc,
                                         print_result=print_result)

        ###Print results?
        if print_result:
            self.epra.report_results(swp_variable='Lj', numeric=True)

        ### Release ANSYS from pyEPR script
        self.pinfo.disconnect()

        ### Log succsessful HFSS Eigenmode Simulation
        self.got_EPR = True

        return self.epra
            
    def run_CapMatirx(self,
                      design_name: str = None,
                      setup_name: str = "Setup",
                      freq_ghz: int = 2,
                      max_delta_f: float = 0.1,
                      max_passes: int = 30,
                      min_passes: int = 1,
                      min_converged: int = 1,
                      pct_refinement: int = 30,
                      auto_increase_solution_order = None,
                      solution_order = None,
                      solver_type = None,
                      qubit_mesh_MaxLength: str = '6um',
                      cavity_mesh_MaxLength: str = '3um',
                      other_mesh_MaxLength: dict = None) -> "pd.DataFrame":
        """
        Runs ANSYS Q3D Simulation and returns the capacitance matrix.

        Args:
            design_name (str, optional): The name of the ANSYS Q3D design.
            setup_name (str, optional): The name of the simulation setup.
            freq_ghz (int, optional): The frequency at which the simulation will be performed, in GHz.
                This parameter should not affect the result of the capacitance matrix.
            max_delta_f (float, optional): The maximum error allowed for the simulation, in GHz.
            max_passes (int, optional): The maximum number of passes.
            min_passes (int, optional): The minimum number of passes.
            min_converged (int, optional): Number of passes to run after reaching `max_delta_f`.
            pct_refinement (int, optional): Maximum allowed refinement for the meshing per pass. 
                                            For example: 30 means 30% of the entire mesh can change per pass.
            auto_increase_solution_order (bool, optional): If True, allows the solver to automatically
                increase the solution order to improve accuracy
            solution_order (int, optional): The solution order to be used for the simulation.
                If auto_increase_solution_order is set to True, this value will be ignored.
            solver_type (str, optional): The solver type to be used for the simulation.
            qubit_mesh_MaxLength (str, optional): The maximum mesh length for the qubit region, in micrometers.
            cavity_mesh_MaxLength (str, optional): The maximum mesh length for the cavity region, in micrometers.
            other_mesh_MaxLength (dict, optional): A dictionary of additional regions and their maximum mesh lengths
                in micrometers. The keys should be the region names, and values should be the maximum mesh length.

        Returns:
            self.cap_matrix (pd.DataFrame): Capacitance matrix.
        """
        
        ### Naming of design & setup
        if (design_name == None):
            design_name = self.q3d_design_name
        
        self.q3d_setup_name = setup_name

        ### Notation
        q3d = self.q3d_renderer.sim.renderer

        ### Start ANSYS
        q3d.start()
        q3d.activate_ansys_design(design_name, 'capacitive')

        ### Render design to ANSYS
        q3d.clean_active_design()
        q3d.render_design(self.selection, self.open_pins)

        ### Add mesh
        # Qubit Mesh
        hfss.modeler.mesh_length('qubit',
                                 [f'cross_{self.qubit_name}'],
                                 MaxLength=qubit_mesh_MaxLength)

        # CPW Mesh
        trace_names = [f'trace_{cpw_name}' for cpw_name in self.cpws_names]
        claw_name = [f'{self.connection_pad_name}_connector_arm_{self.qubit_name}']
        hfss.modeler.mesh_length('cpw',
                                trace_names + claw_name,
                                MaxLength=cavity_mesh_MaxLength)

        # Other Meshes
        if type(other_mesh_MaxLength) == dict:
            for name, MaxLength in other_mesh_MaxLength.items():
                hfss.modeler.mesh_length(name,
                                        [name],
                                        MaxLength=MaxLength)
       
        ### Use pyAEDT to do custom functionality
        self._pyAEDT_functionality(solutiontype='Q3d')

        ### Add Setup
        q3d.add_q3d_setup(
                      name=setup_name,
                      freq_ghz=freq_ghz,
                      save_fields=True,
                      enabled=True,
                      max_passes=max_passes,
                      min_passes=min_passes,
                      min_converged_passes=min_converged,
                      percent_error=max_delta_f,
                      percent_refinement=pct_refinement,
                      auto_increase_solution_order=auto_increase_solution_order,
                      solution_order=solution_order,
                      solver_type=solver_type)

        ### Analyze Setup
        q3d.analyze_setup(setup_name)

        ### Get Capacitance Matrix
        self.q3d_renderer.sim.capacitance_matrix, self.q3d_renderer.sim.units = q3d.get_capacitance_matrix()
        self.cap_matrix = self.q3d_renderer.sim.capacitance_matrix

        ### Release ANSYS
        self.q3d_renderer.sim.renderer = q3d
        self.q3d_renderer.sim.close()

        ### Log succsessful HFSS Eigenmode Simulation
        self.got_CapMatrix = True

        return self.cap_matrix

    def _pyAEDT_functionality(self, solutiontype):
        """Interfaces w/ ANSYS via pyEPR for more custom automation.
        1. Connect to ANSYS
        2. Change Silicon permitivity to 11.45; represents ultra cold silicon.
        3. Checks for prexisting Setups, deletes them...
        """

        if solutiontype == 'Eigenmode':
            projectname = self.hfss_renderer.sim.renderer.pinfo.project_name
            designname = self.hfss_renderer.sim.renderer.get_active_design_name()
            aedt = Hfss(projectname=projectname, 
                        designname=designname, 
                        solution_type=solutiontype,
                        new_desktop_session=False, 
                        close_on_exit=False)
        elif solutiontype == 'Q3d':
            projectname = self.q3d_renderer.sim.renderer.pinfo.project_name
            designname = self.q3d_renderer.sim.renderer.get_active_design_name()
            aedt = Q3d(projectname=projectname,
                       designname=designname,
                       new_desktop_session=False, 
                       close_on_exit=False)
        else:
            raise NotImplementedError('`solutiontype` not implemented yet. Only supports ["Eigenmode", "Q3d"].')
        
        self._ultra_cold_silicon(aedt)
        self._delete_old_setups(aedt)

        aedt.release_desktop(close_projects=False, close_desktop=False)

    def _ultra_cold_silicon(self, aedt):
        """Change silicon properties to ultra cold silicon
        
        Args:
            aedt (pyAEDT Desktop obj)
        """
        materials = aedt.materials
        silicon = materials.checkifmaterialexists('silicon')
        silicon.permittivity = 11.45
        silicon.dielectric_loss_tangent = 1E-7

    def _delete_old_setups(self, aedt):
        """Delete old setups
        
        Args:
            aedt (pyAEDT Desktop obj)
        """
        # Clear setups
        if len(aedt.setups) != 0:
            aedt.setups[0].delete()
        
        

    def _parse_all_results(self, print_result=True):
        raise NotImplementedError('Need to set self.qubit_freq, self.cavity_freq, etc...')

        if (self.got_EPR == False) or (self.got_CapMatrix == False) or (self.got_EigenModes == False):
            raise RuntimeError("Must run `run_EigenMode`, `run_EPR`, and `run_CapMatrix` before calling `_parse_all_results`.")
            
        self._parse_EPR()
        self._parse_CapMatrix()
        self._calc_CouplingStrength()

        if print_result:
            print('________________')
            print(f'Qubit Frequency (f_q) = {self.qubit_freq} Linear MHz')
            print(f'Cavity Frequency (f_cav) = {self.cavity_freq} Linear MHz')
            print(f'Qubit Anharmonicity (f_q) = {self.anharmonicity} Linear MHz')
            print(f'Coupling Strength (g) = {self.coupling_strength} Linear MHz')
            print('________________')

    def _parse_EPR(self) -> dict:
        if (self.got_EPR == False):
            raise RuntimeError("Must run `run_EPR` before calling `parse_EPR`.")
        # Extraction of variables
        omegas = self.epra.get_frequencies()
        chis = self.epra.get_chis()
        other_data = self.epra.data

        
        self.qubit_freq = omegas['0'][0] # Linear MHz
        self.cavity_freq = omegas['1'][0] # Linear MHz
        self.anharmonicity = chis[0][0] # Linear MHz
        self.other_data_EPR = str(other_data)

        package = {'qubit_freq_MHz': self.qubit_freq,
                   'cavity_freq_MHz': self.cavity_freq,
                   'qubit_anharmonicity_MHz': self.anharmonicity,
                   'other_data': self.other_data_EPR}

        return package

    
    def _parse_CapMatrix(self, wavelength: str) -> dict:

        if (wavelength == 'quarter'): # Shorted to ground
            raise NotImplementedError()
            # package = {'cross_to_ground_pF': ,
            #            'cross_to_cpw_pF': ,
            #            'cpw_to_ground_pF': ,
            #            'full_matirx': self.cap_matrix}
        elif (wavelength == "half"): # Open to ground
            raise NotImplementedError()
            # package = {'cross_to_ground_pF': ,
            #        'cross_to_cpw_pF': ,
            #        'full_matirx': self.cap_matrix}
        else:
            raise ValueError("Supported cavity wavelengths are ['quarter', 'half']. \
                             'quarter' refers to a cavity which is shorted to ground. \
                             'half' refers to a cavity which is open to ground.")

        return package
        
    def _calc_CouplingStrength(self) -> dict:
        raise NotImplementedError()

