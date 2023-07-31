import numpy as np
import scqubits as scq

e = 1.602e-19  # elementary charge in C
hbar = 1.054e-34  # reduced Planck constant in Js
Z_0 = 50 # in Ohms

def calc_g(C_crossTOcavity, 
           C_crossTOground, 
           f_cavity, 
           f_qubit, 
           anharmonicity,
           freq_units: str = "GHz"):
  """
  Calculate coupling strength $g$. Returns units used for `f_r` and `f_qubit`:

  Args:
    C_crossTOcavity (float): Capacitance between cross & cavity. Units pF.
    C_crossTOground (float): Capacitance between cross & ground. Units pF.
    f_cavity (float): Resonant frequency of cavity, obtained from pyEPR. Units in linear `freq_units`.
    f_qubit  (float): Resonant frequency of qubit, obtained from pyEPR. Units in linear `freq_units`.
    anharmonicity (float): Qubit anharmonicity, obtained from diagonal of kerr matrix in pyEPR. Units in linear `freq_units`.
    freq_units (str): Choose from ["GHz", "MHz", "kHz", "Hz"]

  Returns:
    g (float): Coupling strength $g$. Units are linear `freq_units`.
    
  """
  C_Sigma = C_crossTOcavity + C_crossTOground + 1.5e-15

  scq.set_units(freq_units)
  if (anharmonicity > 0):
    anharmonicity = -1 * anharmonicity
    print("WARNING: You forgot the negative sign on `anharmonicity`, I've corrected it for you :)")
  EJ, EC = scq.Transmon.find_EJ_EC(E01=f_qubit, anharmonicity=anharmonicity)

  g = (C_crossTOcavity / C_Sigma) * f_cavity * np.sqrt(Z_0 * e**2 / hbar) * np.sqrt(1/2) * (EJ/(8*EC))**(1/4)

  return g_linear 
