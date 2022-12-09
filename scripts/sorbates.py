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
) # TRAPPE

DMSO = dlm.DLMolecule(
    name='DMSO',
    molecule=Atoms(
        'CCSO',
        positions=[
            (-0.517, 1.352, 1.0606),
            (-0.517, -1.352, 1.0606),
            (0,0,0),
            (1.53, 0 0)
        ],
        tags = [0,0,1,2]
    ),
    tags = {0: 'C', 1: 'S', 2: 'O'},
    potentials = {
        0: [125, 3.81, 0.16],
        1: [214, 3.47, 0.139],
        2: [176, 2.83, -0.459]
    }
) # after 10.1039/C4CP05961A

methanol_cgenff = dlm.DLMolecule(
    name = 'MeOH',
    molecules = Atoms(
        'COHHHH',
        positions=[
            (-4.025, 1.427, 0.000),
            (-2.955, 1.427, 0.000),
            (-4.381, 1.041, -0.932),
            (-4.381, 2.427, 0.131),
            (-4.381, 0.814, 0.801),
            (-2.631, 1.983, -0.726)
        ],
        tags = [0,1,2,2,2,3]
    ),
    tags = {0: 'C', 1: 'O', 2: 'H', 3: 'H_O'},
    potentials = {
        0: [39.25, 4.1,   -0.04],
        1: [96.67, 3.53,  -0.650],
        2: [12.08, 2.68,  0.09],
        3: [23.15, 0.449, 0.42]
    }
)

ethanol_cgenff = dlm.DLMolecule(
    name = 'EtOH',
    molecules = Atoms(
        'COHHCHHHH',
        positions=[
            (-4.024,   1.543,  -0.151),
            (-2.624,   1.552,  -0.163),
            (-4.409,   1.128,  -1.109),
            (-4.409,   2.576,   0.000),
            (-4.514,   0.667,   0.992),
            (-2.361,   2.134,  -0.923),
            (-4.145,  -0.373,   0.861),
            (-5.624,   0.655,   1.007),
            (-4.145,   1.064,   1.962)
        ],
        tags = [0,1,2,2,3,4,5,5,5]
    ),
    tags = {0: 'C1', 1: 'O', 2: 'H_C1', 3: 'C2', 4: 'H_O', 5: 'H_C2'},
    potentials = {
        0: [28.18, 4.02, 0.053], 
        1: [96.67, 3.53, -0.65],
        2: [17.61, 2.68, 0.09],
        3: [39.25, 4.1, -0.272], 
        4: [23.15, 0.449, 0.419],
        5: [12.08, 2.68, 0.09]
    }
)

IPA_cgenff = dlm.DLMolecule(
    name = 'IPA',
    molecules = Atoms(
        'COHCCHHHHHHH',
        positions=[
            (-3.968,   1.440,  -0.137),
            (-2.565,   1.470,  -0.114),
            (-4.325,   1.001,  -1.097),
            (-4.521,   2.859,   0.015),
            (-4.464,   0.542,   0.995),
            (-2.270,   1.774,  -1.012),
            (-4.061,  -0.485,   0.864),
            (-5.573,   0.489,   0.987),
            (-4.127,   0.933,   1.979),
            (-4.165,   3.494,  -0.824),
            (-4.182,   3.308,   0.972),
            (-5.632,   2.843,  -0.005)
        ],
        tags = [0,1,2,3,3,4,5,5,5,5,5,5]
    ),
    tags = {0: 'C_O', 1: 'O', 2: 'H_C1', 3: 'C', 4: 'H_O', 5: 'H'},
    potentials = {
        0: [16.1, 4, 0.139],
        1: [96.67, 5.53, -0.641],
        2: [22.65, 2.68, 0.09],
        3: [39.25, 4.1, -0.268],
        4: [23.15, 0.449, 0.408],
        5: [12.08, 2.68, 0.09]
    }
)

DMF_cgenff = dlm.DLMolecule(
    name = 'DMF',
    molecules = Atoms(
        'CONHCCHHHHHH',
        positions=[
            (-1.532,   1.434,  -0.289),
            (-0.306,   1.468,  -0.331),
            (-2.207,   0.699,   0.642),
            (-2.088,   2.032,  -0.999),
            (-1.484,  -0.093,   1.639),
            (-3.654,   0.687,   0.667),
            (-4.066,   1.337,  -0.128),
            (-4.045,   1.072,   1.627),
            (-4.030,  -0.352,   0.568),
            (-0.375,  -0.050,   1.509),
            (-1.806,  -1.157,   1.562),
            (-1.715,   0.268,   2.679)
        ],
        tags = [0,1,2,3,4,4,5,5,5,5,5,5]
    ),
    tags = {0: 'C_O', 1: 'O', 2: 'N', 3: 'H_CO', 4: 'C', 5: 'H_CN'},
    potentials = {
        0: [55.36, 4, 0.423],
        1: [60.39, 3.4, -0.523],
        2: [100.65, 3.7, -0.333],
        3: [23.15, 1.8, 0.079],
        4: [39.25, 4.1, -0.093],
        5: [12.08, 2.68, 0.09]
    }
)

TIP4P = dlm.DLMolecule(
    name = 'TIP4P',
    molecules = Atoms(
        'HOHHe',
        positions=[
            (0.585882,0.756950,0),
            (0,0,0),
            (0.585882,-0.756950,0),
            (0.1546,0,0)
        ],
        tags = [0,1,0,2]
    ),
    tags = {0: 'H', 1: 'O', 2: 'He'},
    potentials = {
        0: [0,0,0.5564],
        1: [93.192,3.1589,0],
        2: [0,0,-1.1128]
    }
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
    'Cyclohexane': cyclohexane,
    'water': TIP4P,
    'cgenff_meoh': methanol_cgenff,
    'cgenff_etoh': ethanol_cgenff,
    'cgenff_IPA': IPA_cgenff,
    'cgenff_DMF': DMF_cgenff,
    'DMSO': DMSO
}
