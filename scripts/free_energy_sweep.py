import logging
import os
import dlmontepython.simtask.dlmonteinterface as interface
import dlmontepython.simtask.measurement as measurement
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
import cif2config as c2c
from glob import glob
import argparse
import pathlib
import json
import sorbates


def Pa_to_katm(pressure: str) -> float:
    try:
        float(pressure)
    except TypeError:
        raise TypeError('Impossible to parse pressure values as floats!')
    return float(pressure) * 9.86923e-9


def str_to_floats(input_string: str) -> list:
    return [float(x) for x in input_string.split(',')]


def pressure_preprocess(input_string: str) -> list:
    return [Pa_to_katm(x) for x in str_to_floats(input_string)]


# Set up the logger, which determines the nature of information output by the machinery in the 'task' package.
# The code below results in logging information being output to stdout

handler = logging.StreamHandler()

# Replacing 'logging.INFO' with 'logging.DEBUGGING' below results in more information being output by the
# logger. Using 'logging.WARNING' results in less information being output: only 'warnings'

measurement.logger.setLevel(logging.INFO)
measurement.logger.addHandler(handler)

# Define external arguments for file input/output locations

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--InputFolder',
                    type=str,
                    action='store',
                    required=False,
                    metavar='INPUT_FOLDER',
                    default='.',
                    help='Location of the framework CIF files.')

parser.add_argument('-o', '--OutputFolder',
                    type=str,
                    action='store',
                    required=False,
                    metavar='OUTPUT_FOLDER',
                    default='.',
                    help='Intended location for output files.')

parser.add_argument('-f', '--FrameworkName',
                    type=str,
                    action='store',
                    required=True,
                    metavar='FRAMEWORK_NAME',
                    help='Name of the framework used (and associated .cif file).')

parser.add_argument('-c', '--GasComposition',
                    action='store',
                    type=eval,
                    required=False,
                    default='{"CO2": 1.0}',
                    metavar='GAS_COMPOSITION',
                    help="Gas composition, as a string json object (e.g. {\"CO2\": 1.0}).")

parser.add_argument('-t', '--Temperature',
                    type=float,
                    action='store',
                    required=False,
                    metavar='TEMPERATURE',
                    default=298.0,
                    help='Specified temperature (in K).')

parser.add_argument('-p', '--Pressure',
                    action='store',
                    required=False,
                    metavar='PRESSURE',
                    type=Pa_to_katm,
                    default=4500,
                    help='Specified pressures in Pa (e.g. 4500).')

parser.add_argument('-q', '--Charges',
                    action='store',
                    required=False,
                    metavar='CHARGES',
                    type=bool,
                    default=True,
                    help='Ewald summation charges - set to False to turn off charges')
parser.add_argument('-nmin',
                    action='store',
                    required=False,
                    metavar='MINIMUM_SORBATES',
                    type=int,
                    default=0,
                    help='Minimum number of sorbates to simulate in free energy simulations. Defaults to 0')
parser.add_argument('-nmax',
                    action='store',
                    required=False,
                    metavar='MAXIMUM_SORBATES',
                    type=int,
                    default=200,
                    help='Maximum number of sorbates to simulate in free energy simulations. Defaults to 200')

args = parser.parse_args()

# Now let's set up the paths to the input and output directories, and check they exist
input_file = pathlib.Path(args.InputFolder, args.FrameworkName).with_suffix('.cif')
assert input_file.exists(), '''Cannot find input file from specified location: {0}
current directory: {1}
{2}
input directory: {3}
{4}'''.format(input_file, os.getcwd(), os.listdir(), args.InputFolder, os.listdir(args.InputFolder))

output_folder = pathlib.Path(args.OutputFolder)
output_folder.mkdir(parents=True, exist_ok=True)  # Makes the directory, if it didn't already exist

logging.debug(args)
logging.info(f"""-------------------
Beginning Automated free energy curve simulation
-------------------
Meta-variables:
-------------------
Input folder: {args.InputFolder}
Output folder: {args.OutputFolder}
Simulation framework: {args.FrameworkName}
Simulation temperature: {args.Temperature}
Simulation composition: {args.GasComposition} (N.B. gas pressures will only change for the first species)
Simulation pressure (in Pa): {args.Pressure}
Minimum number of sorbate molecules to consider: {args.nmin}
Maximum number of sorbate molecules to consider: {args.nmax}
-------------------
""")

# Set up the CONTROL input file from the example in fedsweep_control_generator.py.
control_location = pathlib.Path('/run/CONTROL')

control_obj = fedsweep.TMMCExample
control_obj.use_block.use_statements.pop('ortho')
if not args.Charges:
    control_obj.main_block.statements['noewald'] = 'all'

control_obj.main_block.statements['temperature'] = args.Temperature

control_obj.main_block.moves = fedsweep.define_molecule_movers(list(args.GasComposition.keys())[0],
                                                               molpot=args.Pressure)
with open(control_location, 'w') as f:
    f.write(str(control_obj))

config_field_location = pathlib.Path('/run/')
c2c.create_config_field(input_file=input_file,
                        output_directory=config_field_location,
                        use_cif_hack=True,
                        sorbate_molecules=[sorbates.lookup[list(args.GasComposition.keys())[0]]])

# DEBUG: print out the locations of the input files

logging.info(f'Control file at {control_location} exists? {control_location.exists()}')
config_location = config_field_location / 'CONFIG'
logging.info(f'Config file at {config_location} exists? {config_location.exists()}')
field_location = config_field_location / 'FIELD'
logging.info(f'Field file at {field_location} exists? {field_location.exists()}')

# Set up the relevant TaskInterface object: which tells the low-level machinery in the 'task' package
# which code will be used to perform the simulations, and how to perform various tasks specific to that
# code, e.g. extracting the energy from output files created by the code.
# In this case we use DL_MONTE to perform our simulations; thus the TaskInterface object we will use
# is in fact a DLMonteInterface object (DLMonteInterface is a subclass of TaskInterface). The line
# below sets up a DL_MONTE-specific interface. Note that the interface must know the location of the
# DL_MONTE executable - which is specified as the argument to the DLMonteInterface constructor.

# Set up the relevant TaskInterface object: which tells the low-level machinery in the 'task' package
# which code will be used to perform the simulations, and how to perform various tasks specific to that
# code, e.g. extracting the energy from output files created by the code.
# In this case we use DL_MONTE to perform our simulations; thus the TaskInterface object we will use
# is in fact a DLMonteInterface object (DLMonteInterface is a subclass of TaskInterface). The line
# below sets up a DL_MONTE-specific interface. Note that the interface must know the location of the 
# DL_MONTE executable - which is specified as the argument to the DLMonteInterface constructor.

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

m_template = measurement.Measurement(interface, observables, maxsims=2, precisions=precisions, outputdir='TMMC_test')

# set up the parameters for the sweep

min_val = args.nmin
max_val = args.nmax
stride = 3  # Minimum size of windows

windows = ['fed order param nmols {0} {1} {2} 1 win {3} {4}'.format(
    (max_val - min_val) + 1,
    min_val - 0.5,
    max_val + 0.5,
    win_min - 0.5,
    win_min + stride + 0.5) for win_min in range(min_val, max_val, stride)]

sweep = measurement.MeasurementSweep(
    param="orderparam",
    paramvalues=windows,
    measurement_template=m_template,
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
