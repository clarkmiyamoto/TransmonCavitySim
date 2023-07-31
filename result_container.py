import pandas as pd

class Characteristics:

  def __init__(self, 
               freq_qubit: float, 
               freq_cavity: float, 
               anharmonicity: float, 
               coupling_strength: float,
               units = "Linear GHz"):

    self.freq_qubit = freq_qubit
    self.freq_cavity = freq_cavity
    self.anharmonicity = anharmonicity
    self.coupling_strength = coupling_strength
    self.units = units

  def to_df(self) -> pd.DataFrame:
    """
    Exports object as pd.DataFrame
    """
    d = {
      "freq_qubit"        : self.freq_qubit,
      "freq_cavity"       : self.freq_cavity,
      "anharmonicity"     : self.anharmonicity,
      "coupling_strength" : self.coupling_strength,
      "units"             : self.units
    }

    return pd.DataFrame(d, ignore_index=True)
    
    
