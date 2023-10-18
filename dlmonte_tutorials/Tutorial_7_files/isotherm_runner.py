import dlmontepython.simtask.dlmonteinterface as interface
import dlmontepython.simtask.measurement as measurement
import dlmontepython.simtask.analysis as analysis
import dlmontepython.simtask.task as task

interface = interface.DLMonteInterface("DLMONTE-SRL.X")

energy_obs = task.Observable( ("energy",) )
nmol_obs = task.Observable( ("nmol",1) )
observables = [ energy_obs, nmol_obs ]

precisions = { nmol_obs : 0.2 }

m_template = measurement.Measurement(interface, observables, maxsims=20,
                                     precisions=precisions, outputdir='./fixedprecision')

molchempots = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3]

sweep = measurement.MeasurementSweep(param="molchempot", paramvalues=molchempots,
                                     measurement_template=m_template, outputdir="isotherm_simulation")

sweep.run()