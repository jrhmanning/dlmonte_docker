# See 'fixedtime.py' before considering this script
#
# This script performs multiple DL_MONTE simulations to calculate the mean and uncertainty in the energy and 
# number of molecules in the system reflected in the CONTROL, CONFIG and FIELD files in the current directory 
# - which correspond to a Lennard-Jones fluid - over a range of temperatures to obtain the energy vs. temperature
# and density vs. temperature for this system. The simulations used to calculate the mean and uncertainty at each 
# temperature are constrained for convenience to take 30 seconds. The 10 temperatures from 8000 to 17000 
# (inclusive) in gaps of 1000 are considered; thus this script should take about 300s = 5 mins to complete.

import logging

import dlmontepython.simtask.dlmonteinterface as interface
import dlmontepython.simtask.converge as converge
import dlmontepython.simtask.analysis as analysis
import dlmontepython.htk.sources.dlcontrol as control
import dlmontepython.htk.sources.dlfeddat as feddat
import dlmontepython.htk.sources.dlfedorder as fedorder
import dlmontepython.htk.sources.tmatrix as tmatrix
import dlmontepython.simtask.task as task
import os
import errno
import numpy as np

# Set up the logger, which determines the nature of information output by the machinery in the 'task' package.
# The code below results in logging information being output to stdout

handler = logging.StreamHandler()
handler2 = logging.FileHandler('test.log')

# Replacing 'logging.INFO' with 'logging.DEBUGGING' below results in more information being output by the
# logger. Using 'logging.WARNING' results in less information being output: only 'warnings'

converge.logger.setLevel(logging.INFO)
converge.logger.addHandler(handler)
converge.logger.addHandler(handler2)

# Set up the relevant TaskInterface object: which tells the low-level machinery in the 'task' package
# which code will be used to perform the simulations, and how to perform various tasks specific to that
# code, e.g. extracting the energy from output files created by the code.
# In this case we use DL_MONTE to perform our simulations; thus the TaskInterface object we will use
# is in fact a DLMonteInterface object (DLMonteInterface is a subclass of TaskInterface). The line
# below sets up a DL_MONTE-specific interface. Note that the interface must know the location of the 
# DL_MONTE executable - which is specified as the argument to the DLMonteInterface constructor.

interface = interface.DLMonteInterface("DLMONTE-SRL.X")

# Set up a list of 'observables' to track and analyse. Observables must be Observable objects, and the nature
# of Observable objects may vary between simulation codes. For DL_MONTE only observables corresponding to variables
# output periodically in YAMLDATA are currently supported. For a variable 'foo' specified in the YAMLDATA file
# the corresponding Observable object is returned by the command 'task.Observable( ("foo",) )'. Note the essential
# comma after "foo"! For a variable in YAMLDATA which is an array (e.g., 'nmol'), the observable corresponding to 
# the nth element in the array is returned by the command 'task.Observable( ("foo",n-1) )'. See below: 'energy_obs'
# corresponds to the 'energy' variable in YAMLDATA, and 'nmol_obs' corresponds to the 1st element in the 'nmol'
# array in YAMLDATA (which in fact is the number of molecules belonging to the 1st molecular species)

bias_obs = task.Observable( ("fedbias",) )
tmatrix_obs = task.Observable( ("tmatrix",) )
observables = [ bias_obs, tmatrix_obs ]
precisions = {bias_obs: 1}

# Set up a Measurement object which will determine the nature of the simulations and data analysis at each
# temperature. We specify that no more than 20 simulations will be performed (via the 'axsims' argument), and 
# that the maximum time we will allow over all simulations at a given temperature is 30s (via the 'maxtime' 
# argument).


converge_template = converge.Converge(
    interface,
    observables,
    maxsims=5,
    maxtime=600,
    precisions=precisions
)

# Set up the list of fed values to consider


min_val = 0
max_val = 260
stride = 20  # Minimum size of windows

windows = ['fed order param nmols {0} {1} {2} 1 win {3} {4}'.format(
    (max_val - min_val) + 1,
    min_val - 0.5,
    max_val + 0.5,
    win_min - 0.5,
    win_min + stride + 0.5) for win_min in range(min_val, max_val, stride)]

# Set up a MeasurementSweep object which - which will actually perform the simulations and data analysis.
# Note that all the simulations and output files pertaining to analysis will be created in the directory
# 'fixedtimesweep' (via the 'outputdir' argument). Note also that the control parameter 'param' is set to 
# "temperature" while the list of control parameters to explore 'paramvalues' is set to the the temperatures
# mentioned above. The DLMonteInterface, linked to the MeasurementSweep object via 'interface' in 
# 'measurement_template', thus knows to treat the temperature in the CONTROL file as a control parameter, 
# and explore the temperatures in 'temperatures' accordingly

sweep = converge.ConvergeSweep(
    param='orderparam',
    paramvalues=windows,
    converge_template=converge_template,
    outputdir='fedsweep'
)

# Run the task

sweep.run()
print(sweep.report)

import json
with open('converge_test.json', 'w') as f:
    json.dump(sweep.report ,fp=f, indent=4)


# Once the task is complete in the 'fixedtimesweep' directory there will be directories 'param_8000', 'param_9000',
# etc. containing the data pertaining to each temperature. The contents of each of these directories is the
# same as described in 'fixedtime.py'. However the salient results are contained in the files 'energy_sweep.dat'
# and 'nmol_0_sweep.dat', which contain plots of, respectively, the energy and its uncertainty vs. temperature,
# and the number of molecules in the system and its uncertainty (which is proportional to the density of the system)
# vs. temperature. Note that the 'task' module automatically detects whether the system has equilibrated, and
# also automatically determines a block size to use in block averaging. Thus it is possible that a simulation at one
# of the considered temperatures does not equilibrate within the allowed simulation time, or that it does equilibrate
# and that there is not enough post-equilibration data to obtain the mean or uncertainty. If the mean, say, energy
# cannot be calculated at a given temperature then that temperature is omitted from 'energy_sweep.dat'. If the mean
# can be calculated but the uncertainty cannot (note that the uncertainty requires more post-equilibration to
# calculate) then the uncertainty is quoted in 'energy_sweep.dat' as 'NaN'

for x in sweep.report['Results']['tmatrix'].values():
    try: 
        combined_tmatrix
    except (NameError, UnboundLocalError):
        combined_tmatrix = np.zeros_like(np.array(x['Final fed values']))
    
    combined_tmatrix += np.array(x['Final fed values'])

# Deposit it into a new directory, and set up a simulaiton there

combined_tmatrix = tmatrix.TMATRX(data= combined_tmatrix)

with open('./TMATRX.000', 'w') as f:
    f.write(combined_tmatrix.__str__())

final_simulation = converge.Converge(
    interface,
    observables = observables,
    maxsims=100,
    maxtime=600,
    precisions = precisions,
    outputdir='./TMMC_final'
)

final_simulation.run()
with open('tmmc_final.json', 'w') as f:
    json.dump(final_simulation.report, f, indent=4)