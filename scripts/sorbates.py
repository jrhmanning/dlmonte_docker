import dlmolecule as dlm
from ase import Atoms

Nitrogen = dlm.DLMolecule(
    name='Nitrogen',
    molecule=Atoms(
        'NXN',
        positions=[
            (0, 0, 0),
            (0, 0, 0.55),
            (0, 0, 1.1)
        ],
        tags=[0, 1, 0],
        charges=[-0.482, 0.964, -0.482]),
    tags={0: 'N', 1: 'COM'},
    potentials={0: [36, 3.31, -0.482],
                1: [0, 0, 0.964]
                }  # eps, sigma, q
)  # Nitrogen from TRaPPE

THF_twisted = dlm.DLMolecule(
    name='THF',
    molecule=Atoms(
        'OCCCC',
        positions=[
            (0.0000000000, 0.0000000000, 0.0000000000),
            (1.1689000000, 0.7885000000, 0.0000000000),
            (-1.1689000000, 0.7885000000, 0.0000000000),
            (0.7601760000, 2.2682000000, 0.1226000000),
            (-0.7601760000, 2.2682000000, -0.1226000000)
        ],
        tags=[0, 1, 1, 2, 2],
        charges=[-0.41, 0.16, 0.16, 0.045, 0.045]),
    tags={0: 'O', 1: 'C_0', 2: 'C_C'},
    potentials={0: [190, 2.2, -0.41],
                1: [56.3, 3.88, 0.16],
                2: [56.3, 3.88, 0.045]
                }
)  # Twisted pentagon THF configuration from TRAPPE

THF_envelope = dlm.DLMolecule(
    name='THF',
    molecule=Atoms(
        'OCCCC',
        positions=[
            (0.7700000000, 2.2290000000, 0.1925000000),
            (-0.4115000000, 1.4840000000, 0.0000000000),
            (1.9520000000, 1.4840000000, 0.0000000000),
            (1.5400000000, 0.0000000000, 0.0000000000),
            (0.0000000000, 0.0000000000, 0.0000000000)
        ],
        tags=[0, 1, 1, 2, 2],
        charges=[-0.41, 0.16, 0.16, 0.045, 0.045]),
    tags={0: 'O', 1: 'C_0', 2: 'C_C'},
    potentials={0: [190, 2.2, -0.41],
                1: [56.3, 3.88, 0.16],
                2: [56.3, 3.88, 0.045]
                }
)  # Envelope THF configuration from TRAPPE

MeOH = dlm.DLMolecule(
    name='MeOH',
    molecule=Atoms(
        'HOX',
        positions=[(-1.430, 0.000, 0.000), (0, 0, 0), (0.300, 0.896, 0.000)],
        tags=[0, 1, 2]
    ),
    tags={0: 'H', 1: 'O', 2: 'Me'},
    potentials={
        0: [0, 0, 0.435],
        1: [93, 3.02, -0.7],
        2: [98, 3.75, 0.265]
    }  # eps, sigma, q
)  # TRAPPE

DMF = dlm.DLMolecule(
    name='DMF',
    molecule=Atoms(
        'HCONXX',
        positions=[
            (1.594, 1.019, 0.000),  # H
            (1.130, 0.000, 0.000),  # C
            (1.732, -0.945, 0.000),  # O
            (0, 0, 0),  # N
            (-0.720, 1.247, 0.000),  # Me
            (-0.698, -1.259, 0.000)  # Me
        ],
        tags=[0, 1, 2, 3, 4, 4]),
    tags={0: 'H', 1: 'C', 2: 'O', 3: 'N', 4: 'Me'},
    potentials={
        0: [7.18, 2.2, 0.06],
        1: [47.3, 3.7, 0.45],
        2: [226, 2.96, -0.5],
        3: [144, 3.2, -0.57],
        4: [69, 3.8, 0.28]
    }
)  # TRAPPE

Ace = dlm.DLMolecule(
    name='Acetone',
    molecule=Atoms(
        'XCOX',
        positions=[
            (-0.792, 1.297, 0.000),
            (0.000, 0.000, 0.000),
            (1.229, 0.000, 0.000),
            (-0.792, -1.297, 0.000)],
        tags=[0, 1, 2, 0]),
    tags={0: 'Me', 1: 'C', 2: 'O'},
    potentials={
        0: [98, 3.75, 0],
        1: [40, 3.82, 0.424],
        2: [79, 3.05, -0.424]}
)  # TRAPPE

EtOH = dlm.DLMolecule(
    name='EtOH',
    molecule=Atoms(
        'HOCC',
        positions=[
            (0.300, 0.896, 0.000),
            (0, 0, 0),
            (-1.430, 0.000, 0.000),
            (-1.944, -1.452, 0.000)],
        tags=[0, 1, 2, 3]),
    tags={0: 'H', 1: 'O', 2: 'CH2', 3: 'Me'},
    potentials={
        0: [0, 0, 0.435],
        1: [93, 3.02, -0.7],
        2: [46, 3.95, 0.265],
        3: [98, 3.75, 0]}  # eps, sigma, q
)  # TRAPPE

Dioxane = dlm.DLMolecule(
    name='Dioxane',
    molecule=Atoms(
        'OCCOCC',
        positions=[
            (1.320, 0.000, 0.596),
            (0.780, 1.194, 0.000),
            (-0.780, 1.194, 0.000),
            (0.780, -1.194, 0.000),
            (-0.780, -1.194, 0.000),
            (-1.320, 0.000, -0.596)],
        tags=[0, 1, 1, 0, 1, 1]),
    tags={0: 'O', 1: 'CH2'},
    potentials={
        0: [155, 2.39, -0.38],
        1: [52.5, 3.91, 0.19]}
)  # TRAPPE

ACN = dlm.DLMolecule(
    name='Acetonitrile',
    molecule=Atoms(
        'NCX',
        positions=[
            (1.157, 0, 0),
            (0, 0, 0),
            (-1.54, 0, 0)],
        tags=[0, 1, 2]),
    tags={0: 'N', 1: 'C', 2: 'Me'},
    potentials={
        0: [98, 3.75, 0.269],
        1: [60, 3.55, 0.129],
        2: [60, 2.95, -0.398]}
)  # TRAPPE

Chloroform = dlm.DLMolecule(
    name='Chloroform_Kamath',
    molecule=Atoms(
        'HCCl3',
        positions=[
            (-0.39735, -0.26191, 1.13681),
            (-0.03496, 0.00302, 0.12229),
            (1.74297, -0.00071, 0.10829),
            (-0.63937, -1.18691, -1.05244),
            (-0.62768, 1.62127, -0.31495)],
        tags=[0, 1, 2, 2, 2]),
    tags={0: 'H', 1: 'C', 2: 'Cl'},
    potentials={
        0: [10.06, 2.81, 0.355],
        1: [138.58, 3.41, -0.235],
        2: [68.94, 3.45, -0.04]}
)  # after 10.1021/jp0535238

Chloroform_2 = dlm.DLMolecule(
    name='Chloroform_CDP',
    molecule=Atoms(
        'HCCl3',
        positions=[
            (-0.39735, -0.26191, 1.13681),
            (-0.03496, 0.00302, 0.12229),
            (1.74297, -0.00071, 0.10829),
            (-0.63937, -1.18691, -1.05244),
            (-0.62768, 1.62127, -0.31495)],
        tags=[0, 1, 2, 2, 2]),
    tags={0: 'H', 1: 'C', 2: 'Cl'},
    potentials={
        0: [10.06, 2.81, -0.0551],
        1: [138.58, 3.41, 0.5609],
        2: [68.94, 3.45, -0.1686]}
)  # CDP forcefield, after 10.1021/jp9638550

Chloroform_3 = dlm.DLMolecule(
    name='Chloroform_OPLS',
    molecule=Atoms(
        'XCl3',
        positions=[
            (-0.03496, 0.00302, 0.12229),
            (1.74297, -0.00071, 0.10829),
            (-0.63937, -1.18691, -1.05244),
            (-0.62768, 1.62127, -0.31495)],
        tags=[0, 1, 1, 1]),
    tags={0: 'CH', 1: 'Cl'},
    potentials={
        0: [40.26, 3.8, 0.42],
        1: [150.98, 3.47, -0.14]}
)  # OPLS-UA

Chloroform_4 = dlm.DLMolecule(
    name='Chloroform_gupta',
    molecule=Atoms(
        'HCCl3',
        positions=[
            (-0.39735, -0.26191, 1.13681),
            (-0.03496, 0.00302, 0.12229),
            (1.74297, -0.00071, 0.10829),
            (-0.63937, -1.18691, -1.05244),
            (-0.62768, 1.62127, -0.31495)],
        tags=[0, 1, 2, 2, 2]),
    tags={0: 'H', 1: 'C', 2: 'Cl'},
    potentials={
        0: [0, 0, 0.185],
        1: [37.7417, 3.8, 0. - 0.05],
        2: [150.943, 3.47, -0.045]}
)  # after 10.1016/j.chemphys.2011.03.029

CCl4 = dlm.DLMolecule(
    name='CCl4',
    molecule=Atoms(
        'CCl4',
        positions=[
            (0, 0, 0),
            (1.18, 1.18, 1.18),
            (-1.18, -1.18, 1.18),
            (-1.18, 1.18, -1.18),
            (1.18, -1.18, -1.18)],
        tags=[0, 1, 1, 1, 1]),
    tags={0: 'C', 1: 'Cl'},
    potentials={
        0: [12.37, 2.81, -0.362],
        1: [212.6, 3.25, -0.235]}
)  # after 10.1063/1.4943395

CH2Cl2 = dlm.DLMolecule(
    name='Dichloromethane',
    molecule=Atoms(
        'XCl2',
        positions=[
            (0.000000, 0.762012, 0.000016),
            (-1.474470, -0.215523, 0.000013),
            (1.474468, -0.215525, 0.000014)],
        tags=[0, 1, 1]),
    tags={0: 'CH2', 1: 'Cl'},
    potentials={
        0: [123.34, 3.60, 0.4044],
        1: [123.34, 3.42, -0.2022]}
)  # intermolecular interactions after 10.1021/ct500853q, positions by Mat Tolladay

cyclohexane = dlm.DLMolecule(
    name='cyclohexane',
    molecule=Atoms(
        'CCCCCC',
        positions=[
            (0.0, 0.0, 0.0),
            (1.212, 0.7, 0.0),
            (-1.212, 0.7, 0.0),
            (1.212, 2.1, 0.0),
            (-1.212, 2.1, 0.0),
            (0.0, 2.8, 0.0)
        ],
        tags=[0, 0, 0, 0, 0, 0]
    ),
    tags={0: 'C'},
    potentials={0: [52.5, 3.91, 0]}
) # TRAPPE

CO2 = dlm.DLMolecule(
    name='CO2',
    molecule=Atoms(
        'OCO',
        positions=[
            (-1.16, 0, 0),
            (0, 0, 0),
            (1.16, 0, 0)
        ],
        tags=[0, 1, 0]
    ),
    tags={0: 'O', 1: 'C'},
    potentials={0: [79, 2.8, -0.35],
                1: [27, 3.05, 0.7]}
)

lookup = {
    'Nitrogen': Nitrogen,
    'CO2': CO2,
    'THF': THF_twisted,
    'MeOH': MeOH,
    'EtOH': EtOH,
    'Acetone': Ace,
    'Acetonitrile': ACN,
    'Chloroform': Chloroform_3,
    'CCl4': CCl4,
    'Dichloromethane': CH2Cl2,
    'Cyclohexane': cyclohexane

}
