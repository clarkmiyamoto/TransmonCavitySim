# TransmonCavitySim todo

# Todo
- [ ] fix issue of `run_Cap` asap
- [ ] implement `_calc_CouplingStrength`
- [ ] implement `_parse_CapMatrix`
- [ ] test `parser` functionality
- [ ] screenshots of fields after HFSS sim
- [ ] closing ansys sessions + qiskit-metal and freeing memory once analysis is done via method
- [ ] screenshots of `mesh` should be implemented
- [ ] get reproducible `conda env` and `requirements.txt`

# Finished

- [x] ensure best sim hyper parameters are set to default
    - Yes they default to best! `~Clark`
- [x] issue of c1
    - Yes `c1` was suppposed to be `self.q3d_renderer`. `~Clark`
- [x] make all ANSYS project names end with `unique_date_string` so that there is no issue in execution
    - Implemented. `~Clark`
