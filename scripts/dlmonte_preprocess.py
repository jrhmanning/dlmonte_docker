"""
This script aims to take in a series of high-level simulation arguments for a DL_MONTE molecular simulaiton, producing a set of template input files and an instruction file of arguments required for a DL_MONTE simulation. 

"""

#TODO: there's a bug when -i and -o point to different folders. GHotta fix that!

import logging
import os
import dlmontepython.htk.sources.dlconfig as config
import dlmontepython.htk.sources.dlfield as field
import dlmolecule as dlm
import isotherm_control_generator as isotherm
from copy_files import copy_cif_file
import cif2config as c2c
from glob import glob
import argparse
import pathlib
import json
import sorbates
from typing import Union
import tarfile

def composition_preprocess(input_string:str) -> dict:
    return json.loads(eval(input_string))

# Set up the logger, which determines the nature of information output by the machinery in the 'task' package.
# The code below results in logging information being output to stdout

handler = logging.StreamHandler()

# Replacing 'logging.INFO' with 'logging.DEBUGGING' below results in more information being output by the
# logger. Using 'logging.WARNING' results in less information being output: only 'warnings'

# Define external arguments for file input/output locations

parser = argparse.ArgumentParser()
parser.add_argument('-i','--InputFolder',
                    type=str,
                    action='store',
                    required=False,
                    metavar='INPUT_FOLDER',
                    default='.',
                    help='Location of the framework CIF files.')

parser.add_argument('-o','--OutputFolder',
                    type=str,
                    action='store',
                    required=False,
                    metavar='OUTPUT_FOLDER',
                    default='.',
                    help='Intended location for output files.')

parser.add_argument('-f','--FrameworkName',
                    type=str,
                    action='store',
                    required=True,
                    metavar='FRAMEWORK_NAME',
                    help='Name of the framework used (and associated .cif file). Can also handle float values (for an empty cube of that side length), or "CONFIG" (for predefined CONFIG and FIELD files in hte input folder).')

parser.add_argument('--FrameworkSource',
                    type=str,
                    required=True,
                    action='store',
                    metavar='FRAMEWORK_SOURCE',
                    choices=['local',
                             'CSD',
                             'hMOF',
                             'BWDB',
                             'BW20K',
                             'ARABG',
                             'CoRE2019',
                             'CoRE_DDEC',
                             'CURATED-COF',
                             'baburin_2008',
                             'simperler_2005',
                             'database_zeolite_structures'],
                    help='Source of the CIF file describing the nanoporous material structure.')

parser.add_argument('-c', '--GasComposition',
                    action='store',
                    type=eval,
                    required=False,
                    default='{"CO2": 1.0}',
                    metavar='GAS_COMPOSITION',
                    help="Gas composition, as a string json object (e.g. '{\"CO2\": 1.0}').")

parser.add_argument('-t','--Temperature',
                    type=float,
                    action='store',
                    required=False,
                    metavar='TEMPERATURE',
                    default=298.0,
                    help='Specified temperature (in K).')

parser.add_argument('-q', '--Charges',
                    action='store',
                    required=False,
                    metavar='CHARGES',
                    type=bool,
                    default=True,
                    help='Ewald summation charges - set to False to turn off charges')

parser.add_argument('--Cutoff',
                    action='store',
                    required=False,
                    metavar='CUTOFF_LENGTH',
                    type=float,
                    default=21,
                    help='simulation cutoff length - defaults to 12 Angstrom. Automatically creates supercells based on the minimum image convention for cif file inputs (but not empty boxes).')
args = parser.parse_args()

logging.debug(args)

# Now let's set up the paths to the input and output directories, and check they exist

# Copy files to output directory
copy_cif_file(args.InputFolder, args.FrameworkSource, args.FrameworkName, args.OutputFolder)


input_file = pathlib.Path(args.InputFolder,  args.FrameworkName).with_suffix('.cif')
assert input_file.exists(), '''Cannot find input file from specified location: {0}
current directory: {1}
{2}
input directory: {3}
{4}'''.format(input_file, os.getcwd(), os.listdir(), args.InputFolder, os.listdir(args.InputFolder))

output_folder = pathlib.Path(args.OutputFolder)
output_folder.mkdir(parents=True, exist_ok=True) # Makes the directory, if it didn't already exist

logging.info(f"""-------------------
Generating simulation input files for DL_MONTe
-------------------
Meta-variables:
-------------------
Input folder: {args.InputFolder}
Output folder: {args.OutputFolder}
Simulation framework: {args.FrameworkName}
Simulation temperature: {args.Temperature}
Simulation composition: {args.GasComposition} (N.B. gas pressures will only change for the first species)
-------------------
""")


def create_control_file(outdir, control_obj, charges = False, gascomp = args.GasComposition):

    # Set up the CONTROL input file from the example in isotherm.py.
    control_location = pathlib.Path(f'{outdir}/CONTROL')

    try:
        control_obj.use_block.use_statements.pop('ortho')
    except KeyError:
        pass
    if not charges:
        control_obj.main_block.statements['noewald'] = 'all'

    control_obj.main_block.statements['temperature'] = args.Temperature

    control_obj.main_block.moves = isotherm.define_molecule_movers(list(gascomp.keys())[0])
    with open(control_location, 'w') as f:
        f.write(str(control_obj))

    logging.debug(f'Control file at {control_location} exists? {control_location.exists()}')

    return control_location


control_location = create_control_file(outdir = args.OutputFolder, 
                    control_obj = isotherm.AdsorptionExample, 
                    charges = args.Charges, 
                    gascomp = args.GasComposition
                    )

# Set up the FIELD and CONFIG files from the generator in cif2config.
#TODO: add multiple sorbate functionality
config_field_location = pathlib.Path(args.OutputFolder)

def parse_input_file(
        input_name: str,
        indir = args.InputFolder,  
        outdir: Union[str, pathlib.Path] = config_field_location,
        gascomp = args.GasComposition,
        cutoff = args.Cutoff) -> Union[str,config.CONFIG]:
    # flexibly parses input file information to handle empty cells

    ##Empty box
    try:
        boxsize = float(input_name)
        empty_box = config.CONFIG(
            title="Empty cubic cell",
            level=1,
            dlformat=0,
            vcell=[boxsize,boxsize,boxsize],
            nummol=1000,
            molecules=0
        )

        with open(pathlib.Path(outdir) / 'CONFIG') as f:
            f.write(str(empty_box))
        
        field_file = dlm.make_field(
            framework_molecule=None, 
            sorbate_molecules=[sorbates.lookup[list(gascomp.keys())[0]]],
            cutoff=cutoff,
            sim_title= 'Bulk simulation'   
        )

        with open(pathlib.Path(outdir) / 'FIELD') as f:
            f.write(str(field_file))

        return None

    except (TypeError, ValueError):

        ## preformatted CONFIG and FIELD
        if input_file == "CONFIG":
            try:
                config_file = config.CONFIG.from_file(pathlib.path(indir) / 'CONFIG')
                field_file = field.from_file(pathlib.path(indir) / 'FIELD')

                with open(pathlib.Path(outdir) / 'CONFIG') as f:
                    f.write(str(config_file))

                with open(pathlib.Path(outdir) / 'FIELD') as f:
                    f.write(str(field_file))
            except: 
                raise FileNotFoundError("Can't import both your CONFIG and FIELD file in the input folder - please check your input arguments and try again")
                
            # raise NotImplementedError("Can't handle a preformmated config file without a preformatted FIELD file")

            return None
        else:
            c2c.create_config_field(
                input_file=pathlib.Path(indir, input_name).with_suffix('.cif'),
                output_directory=outdir,
                use_cif_hack=True,
                sorbate_molecules=[sorbates.lookup[list(gascomp.keys())[0]]],
                cutoff=cutoff
            )
        return None
    except:
        raise TypeError(f"Can't interpret the file name {input_name} as a float, CONFIG, or cif file. Please check your spelling and try again")

parse_input_file(
    input_name=args.FrameworkName,
    indir = args.InputFolder,
    outdir = config_field_location, 
    gascomp = args.GasComposition,
    cutoff = args.Cutoff)
# DEBUG: print out the locations of the input files

config_location = config_field_location / 'CONFIG'
logging.debug(f'Config file at {config_location} exists? {config_location.exists()}')
field_location = config_field_location / 'FIELD'
logging.debug(f'Field file at {field_location} exists? {field_location.exists()}')


with tarfile.open(f'{args.OutputFolder}/dlm_preprocess.tgz', 'w:gz') as tar:
    tar.add(config_location, arcname = 'CONFIG')
    tar.add(field_location, arcname = 'FIELD')
    tar.add(control_location, arcname = 'CONTROL')


logging.info('Simulation preprocessing complete!')