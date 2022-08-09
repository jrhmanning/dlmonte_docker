from ase import Atoms
from ase.io import read
import dlmolecule as dlm
import pathlib
import sorbates
import cif_hack

# TODO: add argparsing funcionality
# TODO: hack ASE so it imports charges. Done?

MOF_LJ = {
    'Zn_S': [62.4, 2.46],
    'O_S': [30.19, 3.12],
    'C_S': [47.86, 3.47],
    'H_S': [7.65, 2.85],
    'N_S': [38.95, 3.26],
    'Br_S': [186.19, 3.52],
    'Se_S': [146.3, 3.75],
    'Al_S': [253.9, 4.01],
    'Cu_S': [0, 0],

}

UFF_LJ = {
    'C_S': [52.800000, 3.431000],
    'O_S': [30.200000, 3.118000],
    'H_S': [22.140000, 2.571000],
    'N_S': [34.700000, 3.261000],
    'F_S': [25.140000, 2.997000],
    'Na_S': [15.090000, 2.658000],
    'Mg_S': [55.820000, 2.691000],
    'Al_S': [253.940000, 4.008000],
    'Si_S': [202.150000, 3.826000],
    'P_S': [153.370000, 3.695000],
    'S_S': [137.780000, 3.595000],
    'Cl_S': [114.150000, 3.516000],
    'K_S': [17.600000, 3.396000],
    'Ca_S': [119.680000, 3.028000],
    'Sc_S': [9.550000, 2.936000],
    'Ti_S': [8.550000, 2.829000],
    'V_S': [8.050000, 2.801000],
    'Cr_S': [7.540000, 2.693000],
    'Mn_S': [6.540000, 2.638000],
    'Fe_S': [6.540000, 2.594000],
    'Co_S': [7.040000, 2.559000],
    'Ni_S': [7.540000, 2.525000],
    'Cu_S': [2.510000, 3.114000],
    'Zr_S': [34.700000, 2.783000],
    'Mo_S': [28.160000, 2.719000],
    'Be_S': [42.740000, 2.446000],
    'B_S': [90.510000, 3.638000],
    'Zn_S': [62.350000, 2.462000],
    'Ga_S': [208.690000, 3.905000],
    'Ge_S': [190.580000, 3.813000],
    'As_S': [155.380000, 3.769000],
    'Se_S': [146.330000, 3.746000],
    'Br_S': [126.220000, 3.732000],
    'Rb_S': [20.110000, 3.665000],
    'Sr_S': [118.170000, 3.244000],
    'Y_S': [36.210000, 2.980000],
    'Nb_S': [29.670000, 2.820000],
    'Tc_S': [24.140000, 2.671000],
    'Ru_S': [28.160000, 2.640000],
    'Rh_S': [26.650000, 2.609000],
    'Pd_S': [24.140000, 2.583000],
    'Ag_S': [18.100000, 2.805000],
    'Cd_S': [114.650000, 2.537000],
    'In_S': [301.210000, 3.976000],
    'Sn_S': [285.120000, 3.913000],
    'Sb_S': [225.780000, 3.938000],
    'Te_S': [200.140000, 3.982000],
    'Cs_S': [22.630000, 4.024000],
    'Ba_S': [183.040000, 3.299000],
    'La_S': [8.550000, 3.138000],
    'Ce_S': [6.540000, 3.168000],
    'Pr_S': [5.030000, 3.213000],
    'Nd_S': [5.030000, 3.185000],
    'Pm_S': [4.530000, 3.160000],
    'Sm_S': [4.020000, 3.136000],
    'Eu_S': [4.020000, 3.112000],
    'Gd_S': [4.530000, 3.001000],
    'Tb_S': [3.520000, 3.074000],
    'Dy_S': [3.520000, 3.054000],
    'Ho_S': [3.520000, 3.037000],
    'Er_S': [3.520000, 3.021000],
    'Tm_S': [3.020000, 3.006000],
    'Yb_S': [114.650000, 2.989000],
    'Lu_S': [20.620000, 3.243000],
    'Hf_S': [36.210000, 2.798000],
    'Ta_S': [40.730000, 2.824000],
    'W_S': [33.690000, 2.734000],
    'Re_S': [33.190000, 2.632000],
    'Os_S': [18.610000, 2.780000],
    'Ir_S': [36.710000, 2.530000],
    'Pt_S': [40.230000, 2.454000],
    'Au_S': [19.610000, 2.934000],
    'Hg_S': [193.600000, 2.410000],
    'Tl_S': [341.940000, 3.873000],
    'Pb_S': [333.390000, 3.828000],
    'Bi_S': [260.480000, 3.893000],
    'Po_S': [163.430000, 4.195000],
    'At_S': [142.810000, 4.232000],
    'Rn_S': [124.710000, 4.245000],
    'Ra_S': [203.150000, 3.276000],
    'Ac_S': [16.590000, 3.099000],
    'Th_S': [13.070000, 3.025000],
    'Pa_S': [11.060000, 3.050000],
    'U_S': [11.060000, 3.025000],
    'Np_S': [9.550000, 3.050000],
    'Pu_S': [8.050000, 3.050000],
    'Am_S': [7.040000, 3.012000],
    'Cm_S': [6.540000, 2.963000],
    'Bk_S': [6.540000, 2.975000],
    'Cf_S': [6.540000, 2.952000],
    'Es_S': [6.030000, 2.939000],
    'Fm_S': [6.030000, 2.927000],
    'Md_S': [5.530000, 2.917000],
    'No_S': [5.530000, 2.894000],
    'Lw_S': [5.530000, 2.883000]
}


def create_config_field(input_file, output_directory=pathlib.Path('/run/'), sorbate_molecules=[sorbates.Nitrogen],
                        use_cif_hack = False):
    if use_cif_hack:
        try:
            placeholder = cif_hack.parse_cif_ase('Cu_BTC.cif')
            framework = next(placeholder).get_atoms()
        except:
            raise
            # framework = read(input_file, store_tags=True)
    else:
        framework = read(input_file, store_tags=True)
    sim_title = str(input_file.stem)
    config_location = output_directory / 'CONFIG'
    field_location = output_directory / 'FIELD'

    dl_framework = dlm.from_ase(framework, sim_title, UFF_LJ)
    config = dl_framework.make_config_empty_framework()
    print(sim_title)
    with open(config_location, 'w') as f:
        f.write(str(config))

    with open(field_location, 'w') as f:
        f.write(str(dlm.make_field(dl_framework, sorbate_molecules,
                                   sim_title=f'{sim_title} + {[x.name for x in sorbate_molecules]}')))


if __name__ == "__main__":
    infile = 'Cu_BTC'
    create_config_field(f'./{infile}', output_directory='./', sorbate_molecules=[sorbates.THF_envelope])
