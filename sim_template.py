class AbstractSim:

    def __init__(self, 
                 design,
                 selection: list[str],
                 open_pins: list[tuple[str, str]]
                 ):
        self.design = design

        if self.valid_selection(selection):
            self.selection = selection
            self.open_pins = open_pins

    def valid_selection(self, selection):
        for component_name in self.design.components.keys():
            if component_name in selection:
                raise ValueError(f'`{component_name}` not in `design.components`')
        return True
        
    def get_EigenModes(self):
        """Use HFSS Eigenmode to get frequencies"""
        raise NotImplementedError('Must implement `get_EigenModes` method.')

    def get_EPR(self):
        """Connect to HFSS Eigenmode, run `pyEPR` analysis on it."""
        raise NotImplementedError('Must implement `run_EPR` method.')
    
    def get_CapMatirx(self):
        """Use Q3D to get capacitance matricies"""
        raise NotImplementedError('Must implement `get_CapMatirx` method.')