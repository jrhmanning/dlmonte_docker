"""Interface objects between ASE Atoms objects and DL_MONTE MOLECULE objects

Defines a general molecule class with enough information to produce a CONFIG and FIELD file.
Contains methods to instantiate DLMolecule objects from ASE Atoms objects,

"""

from dlmontepython.htk.sources.dlconfig import CONFIG
import numpy as np
import dlmontepython.htk.sources.dlfield as FIELD
import dlmontepython.htk.sources.dlfieldspecies as FIELDSPECIES
import dlmontepython.htk.sources.dlinteraction as INT
from itertools import combinations, combinations_with_replacement


class DLMolecule:
    # TODO: add repr and string funcitons
    # TODO: build summarizer of the ase object?
    # TODO: decide how to handle charges - none, by atomtype, by atom explicitly?
    # TODO: talk to tom about explicit charge handling

    def __init__(self, name, molecule, tags, potentials):
        self.name = name
        self.molecule = molecule
        self.tags = tags
        self.potentials = potentials

    def _config_atoms_stringify(self):
        '''
        This short function makes an appropriate molecule string for a CONFIG.CONFIG object.
        For each atom in the molecule, the atom name and DL_MONTE tag 'core' is added on line 1
        On line 2, the cartesian coordinates of the atom is then printed.

        :param ase_molecule: (ASE.Atoms) an ASE.Atoms object for one of your simulation molecules
        :param tag_dict: (dict) a dictionary which connects ASE.Atoms tags of struct to a string for DL_Monte
        :return output: (str) a string of the molecule positions for a DL_monte CONFIG.CONFIG object
        '''
        output = []
        for i in self.molecule:
            counter = i.tag
            output.append(
                '{0} core\n {1} {2} {3} 0'.format(self.tags[counter], i.position[0], i.position[1], i.position[2]))
        return output

    def _config_molecule_dict_maker(self, ase_molecule, tag_dict, molecule_name):
        '''
        This function creates a dictionary of all the information you need to write a molecule to a CONFIG.CONFIG object

        :param ase_molecule: (ASE.Atoms) an ASE.Atoms object for one of your simulation molecules
        :param tag_dict: (dict) a dictionary which connects ASE.Atoms tags of struct to a string for DL_Monte
        :param molecule_name: (str) the name of your molecule in the DL_Monte Field file
        :return output: (dict) a dictionary of all of the key information about your molecule to make a CONFIG.CONFIG object
        '''
        output = {}
        output['name'] = self.name
        output['natom'] = len(self.molecule)
        output['atoms'] = self._config_atoms_stringify()
        return output

    def make_config_empty_framework(self, box_name=None, max_molecules=[1000, 1000]):
        '''
        This function makes a CONFIG.CONFIG object for an empty framework, ready for adsorption simulation.
        NB this could be fairly easily extended to handle multiple starting molecules, but that's a job for another day

        :param max_molecules: (list) max num. of each molecule type in your simulation. Default for single adsorption sims
        :return empty_box: (CONFIG.CONFIG) a DL_monte CONFIG.CONFIG object ready to print to a file
        '''
        if not box_name:
            box_name = self.name
        molecules_list = [self._config_molecule_dict_maker(self.molecule, self.tags, self.name)]
        molecule_numbers = [len(molecules_list), *max_molecules]
        empty_box = CONFIG(box_name,
                           0,
                           1,
                           self.molecule.cell,
                           molecule_numbers,
                           molecules_list)
        return empty_box

    def get_maxatom_moltype(self):
        '''
        This function takes an ASE.Atoms object along with a name and makes a FIELDSPECIES.Moltype object for it
        The object suits a MAXATOM type molecule i.e. a rigid framework-style molecule with no self-interactions.

        :param name: (str) the name of your DL_Monte molecule, as written in the FIELD file
        :param ase_molecule: (ASE.Atoms) an ASE.Atoms object for one of your simulation molecules
        :return: (FIELDSPECIES.MolType) a MAXATOM molecule object ready to be put into a DL_Monte FIELD.FIELD object
        '''
        return FIELDSPECIES.MolType(self.name, len(self.molecule))

    def get_field_rigid_molecule(self):
        '''
        This function takes in a molecule as an ASE.Atoms object and returns it as a rgid molecule for a FIELD.FIELD object

        :param name: (str) the name of your DL_Monte molecule, as written in the FIELD file
        :param ase_molecule: (ASE.Atoms) an ASE.Atoms object for one of your simulation molecules
        :param tag_dict: (dict) a dictionary which connects ASE.Atoms tags of struct to a string for DL_Monte
        :return output: (FIELDSPECIES.MolType) a molecule object ready to be put into a DL_Monte FIELD.FIELD object
        '''
        output = FIELDSPECIES.MolType(self.name, len(self.molecule))
        output.rigid = True
        output.exc_coul_ints = True
        output.atoms = self._get_field_atoms()
        return output

    def _get_field_atomtype(self, tag, atom_name):
        '''
        This function creates an Atomtype entry for an individual atom type from your simulation molecule ASE.Atoms objects
        For a given tag, it identifies from the ASE.Atoms object the atom mass and charge, then invokes an Atomtype object

        :param name: (str) the name of your DL_Monte atom, as written in tag_dict
        :param ase_molecule: (ASE.Atoms) an ASE.Atoms object for one of your simulation molecules
        :param tag_dict: (dict) a dictionary which connects ASE.Atoms tags of struct to a string for DL_Monte
        :return: (FIELDspecies.Atomtype) An atomtype object for puttin ginto the Atomtypes section of a FIELD.FIELD object
        '''
        tag_mask = self.molecule.get_tags() == tag
        if np.any(tag_mask):
            test_atom = self.molecule[tag_mask][0]
        else:
            print('Something weird has happened and you have a tag with no assigned atoms!')
            print('Better check your python script!')
            raise
        return FIELDSPECIES.AtomType(atom_name,
                                     'core',
                                     test_atom.mass,
                                     round(test_atom.charge, 5))

    def get_field_atomtypes(self):
        '''
        This function creates a list of FIELDSPECIES.Atomtype objects for a whole molecule, covering each unique atom type

        :param ase_molecule: (ASE.Atoms) an ASE.Atoms object for one of your simulation molecules
        :param tag_dict: (dict) a dictionary which connects ASE.Atoms tags of struct to a string for DL_Monte
        :return output: (list) a list of FIELDSPECIES.Atomtype objects corresponding to the whole molecule
        '''
        output = []
        for tag, name in self.tags.items():
            output.append(self._get_field_atomtype(tag, name))
        return output

    def _get_field_atoms(self):
        '''
        This function takes in the positions of each atom within one of your simulation molecules.
        It then returns a list of FIELDSPECIES.Atom objects, useful for making a FIELDSPECIES.Molecule object

        :param ase_molecule: (ASE.Atoms) an ASE.Atoms object for one of your simulation molecules
        :param tag_dict: (dict) a dictionary which connects ASE.Atoms tags of struct to a string for DL_Monte
        :return output: (list) a list of FIELDSPECIES.Atom objects which go into a FIELDSPECIES.Molecule object
        '''
        output = []
        site = None
        for i in self.molecule:
            output.append(FIELDSPECIES.Atom(
                self.tags[i.tag],
                'core',
                i.position[0],
                i.position[1],
                i.position[2],
                i.mass,
                i.charge,
                site
            ))
        return output

    def get_field_molecule(self):
        '''
        This function takes in an ASE.Atoms molecule object, its name, and a dictionary of the atom tags for DL_Monte
        It then creates a FIELDSPECIES.Molecule object defining the initial atom positions for your FIELD.FIELD object

        :param name: (str) the name of your DL_Monte molecule, as written in the FIELD file
        :param ase_molecule: (ASE.Atoms) an ASE.Atoms object for one of your simulation molecules
        :param tag_dict: (dict) a dictionary which connects ASE.Atoms tags of struct to a string for DL_Monte
        :return output: (FIELDSPECIES.Molecule) a FIELDSPECIES.Molecule object defining your molecule for the FIELD.FIELD object
        '''
        output = FIELDSPECIES.Molecule(self.name)
        atoms_list = self._get_field_atoms()
        for i in atoms_list:
            output.add_atom(i)
        return output

    @staticmethod
    def from_json(json_object):
        import json
        from ase import Atoms
        placeholder = json.loads(json_object)
        thing = DLMolecule(name=placeholder['name'],
                           molecule=Atoms().fromdict(placeholder['molecule']),
                           tags=placeholder['tags'],
                           potentials=placeholder['potentials'])
        return thing

    def to_json(self):
        import json
        from ase import Atoms
        output = {}
        output['name'] = self.name
        output['tags'] = self.tags
        output['potentuaks'] = self.potentials
        output['molecule'] = self.molecule.todict()
        return output

def get_vdw_interactions(interacting_molecules=[], self_excluding_molecules=[]):
    '''
    This function creates a set of Van der Waals interactions for a FIELD.FIELD object.
    First, it takes in a moldict consisting of the following information:
        {'name': (str) the molecule name,
        'molecule' (ASE.Atoms) the Atoms object for the molecule,
        'tags': (dict) the tag_dict for the molecule,
        'LJ_by_tag': (dict) a dict of {tag: [eps, sigma, q]...}}

    It can optionally take in a moldict of a second species to get interaciton sof both species with one another
    Finally, it can take in a moldict for a molecule whose self LJ-interactions you want to exclude i.e. your framework

    Based on these moldicts, it creates a list of all atom pairs minus the self_exclusion atompairs
    It then calculates the Lorentz-Berthelot mixing interactions between these atom pairs.
    For each atompair, it creates an INT.InteractionLJPlus object followed by a FIELD.VDW object from this
    It finally returns a list of FIELD.VDW objects used to make a FIELD.FIELD file


    TODO: Split this up into smaller constituent functions, for better future flexibility
    TODO: TEST instead of explicitly moldict_1 and moldict_2, make it a list to allow flexible numbers of molecules
    :param moldict_1: (dict) a dictionary containing all the key information about your simulaiton molecule
    :param moldict_2: (dict) as moldict_1, if you want multiples
    :param self_exclusion_moldict: (dict) as moldict_1, if you want to exclude self-interactions for a molecules
    :return output: (list) a list of FIELD.VDW objects used for making a FIELD.FIELD file
    '''

    output = []
    interactions_by_name = {}
    exclusion_interactions = {}

    # atompairs = set(list(combinations_with_replacement([x.molecule.atoms for x in interacting_molecules])))\
    #            - [set(list(combinations_with_replacement([x.molecule.atoms for x in self_excluding_molecules])))]

    for molecule in interacting_molecules:
        # self interaction terms
        for tag, name in molecule.tags.items():
            interactions_by_name[name] = (
                (molecule._get_field_atomtype(tag, name), molecule.potentials[tag]))

    for molecule in self_excluding_molecules:
        for tag, name in molecule.tags.items():
            exclusion_interactions[name] = (
                (molecule._get_field_atomtype(tag, name), molecule.potentials[tag]))

    atompairs = set([i for i in combinations_with_replacement(interactions_by_name.keys(), 2)])

    exclusion_atompairs = set([i for i in combinations_with_replacement(exclusion_interactions.keys(), 2)])
    tested_atompairs = list(atompairs.difference(exclusion_atompairs))
    print(tested_atompairs)

    for i in tested_atompairs:
        name1 = interactions_by_name[i[0]][0]
        name2 = interactions_by_name[i[1]][0]
        sigma = round(0.5 * (interactions_by_name[i[0]][1][1] + interactions_by_name[i[1]][1][1]), 3)
        epsilon = round(np.sqrt(interactions_by_name[i[0]][1][0] * interactions_by_name[i[1]][1][0]), 3)
        if epsilon > 0:
            output.append(FIELD.VDW(name1, name2, INT.InteractionLJLRC(epsilon, sigma)))

    return output


def make_field(framework_molecule, sorbate_molecules=[], cutoff=12, sim_title='Test',
               sorbent_self_interactions=False):
    '''
    This function makes an entire FIELD.FIELD object for you.
    It pulls together a lot of functions to make a complete FIELD file for you.
    To do this, it needs a moelcule dictionary for the framework you're using as well as for each sorbate you're using.
    The molecule dictionary takes the following form:
        {'name': (str) the molecule name,
        'molecule' (ASE.Atoms) the Atoms object for the molecule,
        'tags': (dict) the tag_dict for the molecule,
        'LJ_by_tag': (dict) a dict of {tag: [eps, sigma, q]...}}

    First it creates a FIELD.FIELD object and fills in some universal parameters
    Then it creates atomtypes and moltypes for the framework molecule
    After this, it does the same for sorbate molecules
    Foially, it creates a list of VDW parameters for all of the molecules you're simulating

    :param framework_molecule_dict: (dict) a molecule dictionary of your framework
    :param sorbate_molecule_list: (list) a list of molecule dictionaries for all of your sorbates
    :param cutoff: (float) your simulation cutoff, in A
    :param sim_title: (str) the name of your simulation
    :param sorbent_self_interactions: (bool) True if you want to consider framework self-LJ interactions (untested)
    :return output: (FIELD.FIELD) a FIELD.FIELD object of your simulation parameters
    '''
    # Preamble
    output = FIELD.FIELD()
    output.description = sim_title
    output.cutoff = cutoff
    output.units = 'K'

    # Framework atomtypes and moltype
    for atomtype in framework_molecule.get_field_atomtypes():
        output.atomtypes.append(atomtype)
    output.moltypes.append(framework_molecule.get_maxatom_moltype())

    # Sorbate atomtypes and moltypes
    for sorbate in sorbate_molecules:
        for atomtype in sorbate.get_field_atomtypes():
            output.atomtypes.append(atomtype)

        output.moltypes.append(sorbate.get_field_rigid_molecule())

    # VDW interactions
    if len(sorbate_molecules) > 0:
        for interaction in get_vdw_interactions([framework_molecule, *sorbate_molecules], [framework_molecule]):
            output.vdw.append(interaction)

    elif sorbent_self_interactions:
        for interaction in get_vdw_interactions(framework_molecule):
            output.vdw.append(interaction)

    print('number of interactions: ', len(output.vdw), 'max:',
          int(0.5 * (len(output.atomtypes) * (len(output.atomtypes) + 1))))
    return output


# region from ase functions:

def calculate_supercell(ase_object, cutoff=12):
    supercell = []
    for i in range(3):
        dim = np.linalg.norm(ase_object.cell[i])
        if dim < cutoff * 2:
            print('oh no! my x value is too small! ({0})'.format(dim))
            supercell.append(np.ceil((cutoff * 2) / dim))
        else:
            supercell.append(1)
    print(supercell)
    assert len(supercell) == 3
    superstructure = ase_object * (int(supercell[0]), int(supercell[1]), int(supercell[2]))
    return superstructure


def make_framework_indices(struct):
    '''
    This function creates a dictionary of indices in your ASE.Atoms framework.
    For each element type in the framework, a key is made corresponding to the element number.
    This has values 'symbol' (the element symbol as a string), and 'idx_mask'
    idx_mask is an ndarray which is true when the Atoms index == your element, otherwise false
    :param struct: (ASE.Atoms) your framework as an ASE Atoms object
    :return structure_dictionary: (dict) a dictionary of element types, symbols, and their corresponding indices
    '''
    element_types = np.unique(struct.numbers)
    print('Element types detected:', element_types)
    structure_dictionary = {}
    for i in element_types:
        structure_dictionary[i] = {}
        idx_mask = struct.numbers == i
        assert len(idx_mask) == len(struct)
        structure_dictionary[i]['symbol'] = struct[np.where(idx_mask)[0][0]].symbol
        structure_dictionary[i]['idx_mask'] = idx_mask
    return structure_dictionary


def framework_charges(struct, charges=[], tol=None):
    '''
    This simple function takes in your structure as an ASE.Atoms object and list of charges.
    It optionally rounds all the items in the list of charges according to your tolerances
    Then, assuming the charges and structure have the same length, it applies the charges to the structure

    :param struct: (ASE.Atoms) your framework structure as an ASE.Atoms object
    :param charges: (list) a list of charges for each atom in your structure, in theh same order as struct
    :param tol: (int) the number of decimal places to round your charges to, defaults to non-rounded
    :return struct: (ASE.Atoms) the same structure as was input, but with atom charges according to the charges list
    '''
    if tol is not None:
        charges = [round(i, tol) for i in charges]
        assert len(struct) == len(charges)
    struct.set_initial_charges(charges)
    return struct


def assign_all_framework_tags(struct, structure_dictionary, element_label_addition='S'):
    '''
    This parent function takes in an ASE.Atoms framework structure and assigns unique tags to each atom type
    Atom types are subdivided by Lennard Jones rules ('moieties') and coulombic point charges, and so are tags
    It takes in a structure, the structure dictionary (cf. make_framework_indices), and a string to signify it's a sorbent

    Based on the element and option moiety subdivision in the structure dictionary, tags are assigned to unique charges
    A list of unique tags is then made (tag_list) which can be applied directly on the to structure ASE.Atoms object
    This corresponds to keys in a dictionary (tag_dict), which can be read when writing a DL_MONTE FIELD and CONFIG file

    :param struct: (ASE.Atoms) your framework structure as an ASE.Atoms object
    :param structure_dictionary: (dict) a dictionary of element types, symbols, and their corresponding indices
    :return tag_dict: (dict) a dictionary which connects ASE.Atoms tags of struct to a string for DL_Monte
    :return tag_list: (ndarray) a list of the same length as struct, with the assigned tag values in it
    '''
    tag_list = np.zeros_like(struct)
    tag_dict = {}
    tag_counter = 1
    for element, element_dict in structure_dictionary.items():
        print('Making tags for element', element)
        if 'moieties' in element_dict.keys():
            print('Subdividing element {0} into {1} moieties'.format(element, len(element_dict['moieties'])))
            for moiety, indices in element_dict['moieties'].items():
                print('Making tags for element', moiety)
                print('Assigning tag', tag_counter)
                subtag_dict, subtag_list, tag_counter = assign_sub_tags_by_charge(struct, indices, element, moiety,
                                                                                  tag_counter)
                tag_dict = {**tag_dict, **subtag_dict}
                tag_list += subtag_list
        else:
            print('Assigning tag', tag_counter)
            indices = element_dict['idx_mask']
            moiety = '{0}_{1}'.format(element_dict['symbol'], element_label_addition)
            subtag_dict, subtag_list, tag_counter = assign_sub_tags_by_charge(struct, indices, element, moiety,
                                                                              tag_counter)
            tag_dict = {**tag_dict, **subtag_dict}
            tag_list += subtag_list
    print('Final tag assignments:', tag_dict)
    return tag_dict, tag_list


def assign_sub_tags_by_charge(struct, indices, element, name, tag_counter):
    '''
    This function assigns tags to your ase object, according to the number of unique charge states present.
    It takes the framework structure (struct), and first finds each unique charge in the index list.
    Then, for each unique charge state, it adds a new tag to the Atoms object.
    This tag then corresponds to a dictionary key, whose value is the string written to DL_MONTE as the atom type.

    Realistically, it's a very overcomplicated way of automatically splitting all atoms down into the following way:
    [Element]_[index of LJ parameters]_[index of charge state]

    :param struct: (ASE.Atoms) your framework structure as an ASE.Atoms object
    :param indices: (ndarray) ndarray which is true when the Atoms index == your element/moeity, otherwise false
    :param element: (int) the atomic number of your element/moiety of choice
    :param name: (str) the name of your element/moiety of choice e.g. O_COOH
    :param tag_counter: (int) a counter for the number of tags assigned already
    :return sub_tag_dict: (dict) a dictionary which connects ASE.Atoms tags of struct to a string for DL_Monte
    :return sub_tag_list: (ndarray) A list of the same length as struct, with the assigned tag values in it
    :return tag_counter: (int) a counter to make sure there's no duplication of ASE.Atoms.tag values between objects
    '''
    sub_tag_list = np.zeros_like(struct)
    sub_tag_dict = {}
    assert len(indices) == len(struct)
    charge_states = np.unique(struct[indices].get_initial_charges())
    for count, k in enumerate(charge_states):
        k_indices = (struct.get_initial_charges() == k) & (struct.numbers == element) & (indices == True)
        # NB this doesn't protect against 2 atom types having exactly the same charge. How do I do that?
        # But also, is it a problem if the same element has the same charge type? isn't that what I'm trying to do here?
        if len(charge_states) > 1:
            sub_tag_dict[tag_counter] = '{0}_{1}'.format(name, count)
        else:
            sub_tag_dict[tag_counter] = '{0}'.format(name)
        sub_tag_list[k_indices] = tag_counter
        tag_counter += 1
    return sub_tag_dict, sub_tag_list, tag_counter


# endregion


def from_ase(ase_object, name, interactions_dict, cutoff=12, heterogeneous_vdw=False):
    if heterogeneous_vdw:
        raise NotImplementedError
    superstructure = calculate_supercell(ase_object, cutoff)
    atom_masks = make_framework_indices(superstructure)
    charges = [0 for _ in range(len(superstructure))]
    charged_structure = framework_charges(superstructure, charges)

    tag_dict, tag_mask = assign_all_framework_tags(charged_structure, atom_masks)
    charged_structure.set_tags(tag_mask)

    molecule_LJ_by_tag = {}
    for key, value in tag_dict.items():
        split_char = '_'
        LJ_key = split_char.join(value.split(split_char)[:2])
        # print("DEBUGGING: ", key, value, LJ_key)
        molecule_LJ_by_tag[key] = interactions_dict[LJ_key]
    print(molecule_LJ_by_tag)

    output = DLMolecule(
        name=name,
        molecule=charged_structure,
        tags=tag_dict,
        potentials=molecule_LJ_by_tag
    )
    return output


if __name__ == "__main__":
    pass
