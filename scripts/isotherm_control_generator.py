#!/usr/bin/env python
# coding: utf-8

from dlmontepython.htk.sources import dlcontrol, dlmove
from collections import OrderedDict
from random import randint

# A worked example of a control file (as a python object)
AdsorptionExample = dlcontrol.CONTROL(
    title='Adsorption worked example',

    use_block=dlcontrol.UseBlock(  # overarching simulation parameters e.g. free energy schemes
        use_statements=OrderedDict(
            {
                'gaspressure': None,  # Use fugacity, not chemical potential
                'ortho': None  # Use various speedups of an orthogonal simulation cell
            }
        )
    ),

    main_block=dlcontrol.MainBlock(  # General simulation parameters
        statements=OrderedDict(
            {
                # Defining random seeds
                'seeds': OrderedDict(
                    {
                        'seed0': randint(0, 99),
                        'seed1': randint(0, 99),
                        'seed2': randint(0, 99),
                        'seed3': randint(0, 99)
                    }
                ),

                # statements on simulation length and output frequency
                'archiveformat': 'dlmonte',  # Trajectory output file format
                'equilibration': 0,  # Simulation burn in steps
                'steps': int(1e6),  # Total simulation steps
                'check': int(1e3),  # Frequency of checking internal energy calculations
                'stack': int(1e3),  # Stack size for block averaging
                'yamldata': int(1e3),  # Frequency of yamldata file output
                'print': int(1e4),  # Frequency of configuration printing output

                # Statements on simulation physical parameters
                'temperature': 77.,  # In Kelvin
                'ewald precision': 1e-6,

                # Statements on Monte Carlo move updates
                'acceptatmmoveupdate': 200,
                'acceptmolrotupdate': 200
            }
        ),
        # MC move definition
        moves=[
            dlmove.InsertMoleculeMove(
                pfreq=50,  # proportional frequency of move type
                rmin=0.7,  # auto-reject distance for insert attemts, in angstrom
                movers=[
                    {
                        'id': 'Nitrogen',  # Molecule name, cross-references against FIELD and CONFIG files
                        'molpot': 5e-7  # chemical potential as defined in the use block
                    }
                ]
            ),
            dlmove.MoleculeMove(
                pfreq=25,
                movers=[
                    {
                        'id': 'Nitrogen'
                    }
                ]
            ),
            dlmove.RotateMoleculeMove(
                pfreq=25,
                movers=[
                    {
                        'id': 'Nitrogen'
                    }
                ]
            )
        ],

        # Writing information to other files
        samples=OrderedDict(
            {
                'coords': OrderedDict(  # Atom coordinates
                    {
                        'nfreq': int(1e6)
                    }
                )
            }
        )
    )

)

if __name__ == '__main__':
    AdsorptionExample.use_block.use_statements.pop('ortho')
    AdsorptionExample.main_block.statements['noewald'] = 'all'
    with open('./input_files/CONTROL', 'w') as f:
        f.write(str(AdsorptionExample))


