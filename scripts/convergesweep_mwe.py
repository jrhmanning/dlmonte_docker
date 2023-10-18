import logging
import os
import dlmontepython.simtask.dlmonteinterface as interface
import dlmontepython.simtask.measurement as measurement
import dlmontepython.simtask.converge as converge
import dlmontepython.simtask.analysis as analysis
import dlmontepython.simtask.task as task
import dlmontepython.htk.sources.dlfeddat as fed
import dlmontepython.htk.sources.dlfedorder as fedorder
import dlmontepython.htk.sources.dlcontrol as control
import dlmontepython.htk.sources.dlfedmethod as fedmethod
import fedsweep_control_generator as fedsweep
from glob import glob
import pandas as pd
import errno

import numpy as np



interface = interface.DLMonteInterface("/usr/local/bin/DLMONTE-SRL.X")

# Set up a list of 'observables' to track and analyse. Observables must be Observable objects, and the nature
# of Observable objects may vary between simulation codes. For DL_MONTE only observables corresponding to variables
# output periodically in YAMLDATA are currently supported. For a variable 'foo' specified in the YAMLDATA file
# the corresponding Observable object is returned by the command 'task.Observable( ("foo",) )'. Note the essential
# comma after "foo"! For a variable in YAMLDATA which is an array (e.g., 'nmol'), the observable corresponding to 
# the nth element in the array is returned by the command 'task.Observable( ("foo",n-1) )'. See below: 'energy_obs'
# corresponds to the 'energy' variable in YAMLDATA, and 'nmol_obs' corresponds to the 1st element in the 'nmol'
# array in YAMLDATA (which in fact is the number of molecules belonging to the 1st molecular species)

bias_obs = task.Observable(("fedbias",))
observables = [bias_obs]

# Set up a dictionary containing the threshold precisions to which to determine the observables. Here we only
# specify a threshold precision for the energy observable: the simulations will terminate when the energy is
# determined to an uncertainty less than 0.2 (energy units).

precisions = {bias_obs: 1}

# Set up the relevant Measurement object - which will actually perform the simulations and data analysis.
# Note that all the simulations and output files pertaining to analysis will be created in the directory
# 'fixedprecision' (via the 'outputdir' argument), and that we will impose a threshold precision on the
# energy as implied by the dictionary 'precisions' described above. Furthermore, we specify that no more 
# than 20 simulations will ever be performed (via the 'maxsims' argument)

c_template = converge.Converge(interface, observables, maxsims=3, precisions=precisions, outputdir='TMMC_test')

# set up the parameters for the sweep

min_val = 0
max_val = 480
stride = int(np.ceil((max_val-min_val)/50))  # Minimum size of windows

windows = ['fed order param nmols {0} {1} {2} 1 win {3} {4}'.format(
    (max_val - min_val) + 1,
    min_val - 0.5,
    max_val + 0.5,
    win_min - 0.5,
    win_min + stride + 0.5) for win_min in range(min_val, max_val, stride)]

sweep = converge.ConvergeSweep(
    param="orderparam",
    paramvalues=windows,
    converge_template=c_template,
    outputdir="TMMC_test"
)
# Run the task

sweep.run()

# Now combine the TMATRX files together

output_df = pd.DataFrame({
    0: [0. for _ in range(min_val, max_val + 1)],
    1: [0. for _ in range(min_val, max_val + 1)],
    2: [0. for _ in range(min_val, max_val + 1)],
    3: [0. for _ in range(min_val, max_val + 1)]
})

logging.info(f'TMATRX data length: {len(output_df)}')

for x in windows:
    logging.info(f'Window {x} of {len(windows)}')
    final_sim = sorted([int(y.rsplit('_', 1)[-1]) for y in glob(f'TMMC_test/param_{x}/sim*')], reverse=True)[0]
    logging.info(x, '\tfinal sim = ', final_sim)
    target_file = f'TMMC_test/param_{x}/sim_{final_sim}/TMATRX.000'
    with open(target_file, 'r') as f:
        f.readline()
        f.readline()
        working = [x.strip().split() for x in f.readlines()]
    working_df = pd.DataFrame(working, dtype=float)
    output_df = output_df.add(working_df)

output_str = '# Python generated TMATRX file\n# Tridiagonal format\n'
for _, y in output_df.iterrows():
    output_str += '    '.join([f'{y2:23.17E}' for y2 in y])
    output_str += '\n'

# Deposit it into a new directory, and set up a simulaiton there

try:
    os.makedirs('./TMMC_final')
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

with open('./TMMC_final/TMATRX', 'w') as f:
    f.write(output_str)

# run the simulation, avoiding the Measurement.run() method because it won't handle extra input files

interface.copy_input_files('.', './TMMC_final')
final_control = control.from_file('./TMMC_final/CONTROL')
final_control.use_block.fed_block.method.mode = "res"
final_control.use_block.fed_block.orderparam = fedorder.from_string(
    'fed order param nmols {0} {1} {2} 1'.format(
        (max_val - min_val) + 1,
        min_val - 0.5,
        max_val + 0.5
    )
)
with open('./TMMC_final/CONTROL', 'w') as f:
    f.write(str(final_control))

interface.run_sim('./TMMC_final')

# Extract the final feddat file from the directory

final_fed_data = fed.load('./TMMC_final')

with open(f'{args.OutputFolder}/final_feddat.dat', 'w') as f:
    f.write(str(final_fed_data))
