#!/usr/bin/env python
# coding: utf-8

from dlmontepython.htk.sources import dlcontrol, dlmove
from dlmontepython.htk.sources import dlfedorder, dlfedmethod, dlfedflavour
from collections import OrderedDict
from random import randint

def define_molecule_movers(molname: str, molpot: float = 5e-7) -> list:
    insert = dlmove.InsertMoleculeMove(
        pfreq=50,  # proportional frequency of move type
        rmin=0.7,  # auto-reject distance for insert attempts, in angstrom
        movers=[{'id': molname, # Molecule name, cross-references against FIELD and CONFIG files
                 'molpot': molpot}]
    )
    move = dlmove.MoleculeMove(
        pfreq=25,
        movers=[{'id': molname}]
    )
    rot = dlmove.RotateMoleculeMove(
        pfreq=25,
        movers=[{'id': molname}]
    )
    output = [insert, move, rot]
    return output


# A worked example of a control file (as a python object)
TMMCExample = dlcontrol.CONTROL(
    title='Transition matrix Monte Carlo worked example',

    use_block=dlcontrol.UseBlock(  # overarching simulation parameters e.g. free energy schemes
        use_statements=OrderedDict(
            {
                'gaspressure': None,  # Use fugacity, not chemical potential
                'ortho': None  # Use various speedups of an orthogonal simulation cell
            }
        ),
        fed_block = dlcontrol.FEDBlock(
            flavour= dlfedflavour.Generic(),
            method= dlfedmethod.TransitionMatrix(nout=int(1e6),
                                                 n_upd=int(1e3),
                                                 mode='new',
                                                 tri=True),
            orderparam= dlfedorder.FEDOrderParameter(name='nmols',
                                                     ngrid=201,
                                                     xmin=-0.5,
                                                     xmax=200.5,
                                                     npow=1)
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
        moves= define_molecule_movers('Nitrogen'),

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
