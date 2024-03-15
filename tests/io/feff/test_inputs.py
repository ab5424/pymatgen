from __future__ import annotations

import os
import unittest

from numpy.testing import assert_allclose
from pytest import approx

from pymatgen.core import Molecule, Structure
from pymatgen.io.feff.inputs import Atoms, Header, Paths, Potential, Tags
from pymatgen.util.testing import TEST_FILES_DIR

FEFF_TEST_DIR = f"{TEST_FILES_DIR}/feff"

header_string = """* This FEFF.inp file generated by pymatgen
TITLE comment: From cif file
TITLE Source:  CoO19128.cif
TITLE Structure Summary:  Co2 O2
TITLE Reduced formula:  CoO
TITLE space group: (Cmc2_1), space number:  (36)
TITLE abc:  3.297078   3.297078   5.254213
TITLE angles: 90.000000  90.000000 120.000000
TITLE sites: 4
* 1 Co     0.666666     0.333332     0.496324
* 2 Co     0.333333     0.666667     0.996324
* 3 O     0.666666     0.333332     0.878676
* 4 O     0.333333     0.666667     0.378675"""


class TestHeader(unittest.TestCase):
    def test_init(self):
        filepath = f"{FEFF_TEST_DIR}/HEADER"
        header = Header.header_string_from_file(filepath)
        h_lines = header.splitlines()
        h_str = header_string.splitlines()
        for idx, line in enumerate(h_lines):
            assert line == h_str[idx]
        assert header_string.splitlines() == header.splitlines(), "Failed to read HEADER file"

    def test_from_str(self):
        header = Header.from_str(header_string)
        assert header.struct.reduced_formula == "CoO", "Failed to generate structure from HEADER string"

    def test_get_str(self):
        cif_file = f"{TEST_FILES_DIR}/CoO19128.cif"
        h = Header.from_cif_file(cif_file)
        head = str(h)
        assert (
            head.splitlines()[3].split()[-1] == header_string.splitlines()[3].split()[-1]
        ), "Failed to generate HEADER from structure"

    def test_as_dict_and_from_dict(self):
        file_name = f"{FEFF_TEST_DIR}/HEADER"
        header = Header.from_file(file_name)
        dct = header.as_dict()
        header2 = Header.from_dict(dct)
        assert str(header) == str(header2), "Header failed to and from dict test"


class TestFeffAtoms(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.structure = Structure.from_file(f"{TEST_FILES_DIR}/CoO19128.cif")
        cls.atoms = Atoms(cls.structure, "O", 12.0)

    def test_absorbing_atom(self):
        atoms_1 = Atoms(self.structure, 0, 10.0)
        atoms_2 = Atoms(self.structure, 2, 10.0)
        assert atoms_1.absorbing_atom == "Co"
        assert atoms_2.absorbing_atom == "O"

    def test_single_absorbing_atom(self):
        """
        When there is only one absorbing atom in the structure, it should not appear
        in the pot_dict to avoid an error.
        """
        # one Zn+2, 9 triflate, plus water
        xyz = f"{FEFF_TEST_DIR}/feff_radial_shell.xyz"
        mol = Molecule.from_file(xyz)
        mol.set_charge_and_spin(-7)
        atoms = Atoms(mol, "Zn", 9)
        # Zn should not appear in the pot_dict
        assert not atoms.pot_dict.get("Zn", False)

    def test_absorber_line(self):
        atoms_lines = self.atoms.get_lines()
        # x
        assert float(atoms_lines[0][0]) == approx(0.0)
        # y
        assert float(atoms_lines[0][1]) == approx(0.0)
        # z
        assert float(atoms_lines[0][2]) == approx(0.0)
        # ipot
        assert int(atoms_lines[0][3]) == 0
        # atom symbol
        assert atoms_lines[0][4] == "O"
        # distance
        assert float(atoms_lines[0][5]) == approx(0.0)
        # id
        assert int(atoms_lines[0][6]) == 0

    def test_distances(self):
        atoms_1 = self.atoms.get_lines()
        distances_1 = [float(a[5]) for a in atoms_1]
        atoms_2 = Atoms.atoms_string_from_file(f"{FEFF_TEST_DIR}/ATOMS")
        atoms_2 = atoms_2.splitlines()[3:]
        distances_2 = [float(a.split()[5]) for a in atoms_2]
        assert_allclose(distances_1, distances_2, rtol=1e-5)

    def test_atoms_from_file(self):
        filepath = f"{FEFF_TEST_DIR}/ATOMS"
        atoms = Atoms.atoms_string_from_file(filepath)
        assert atoms.splitlines()[3].split()[4] == "O", "failed to read ATOMS file"

    def test_get_str(self):
        header = Header.from_str(header_string)
        struct = header.struct
        central_atom = "O"
        atoms = Atoms(struct, central_atom, radius=10.0)
        atoms = str(atoms)
        assert atoms.splitlines()[3].split()[4] == central_atom, "failed to create ATOMS string"

    def test_as_dict_and_from_dict(self):
        file_name = f"{FEFF_TEST_DIR}/HEADER"
        header = Header.from_file(file_name)
        struct = header.struct
        atoms = Atoms(struct, "O", radius=10.0)
        dct = atoms.as_dict()
        atoms2 = Atoms.from_dict(dct)
        assert str(atoms) == str(atoms2), "Atoms failed to and from dict test"

    def test_cluster_from_file(self):
        self.atoms.write_file("ATOMS_test")
        mol_1 = Atoms.cluster_from_file("ATOMS_test")
        mol_2 = Atoms.cluster_from_file(f"{FEFF_TEST_DIR}/ATOMS")
        assert mol_1.formula == mol_2.formula
        assert len(mol_1) == len(mol_2)
        os.remove("ATOMS_test")

    def test_atom_num(self):
        filepath = f"{FEFF_TEST_DIR}/Pt37_atoms.inp.gz"
        atoms = Atoms.cluster_from_file(filepath)
        assert len(atoms) == 37
        assert atoms.formula == "Pt37"


class TestFeffTags(unittest.TestCase):
    def test_init(self):
        filepath = f"{FEFF_TEST_DIR}/PARAMETERS"
        parameters = Tags.from_file(filepath)
        parameters["RPATH"] = 10
        assert parameters["COREHOLE"] == "Fsr", "Failed to read PARAMETERS file"
        assert parameters["LDOS"] == [-30.0, 15.0, 0.1], "Failed to read PARAMETERS file"

    def test_diff(self):
        filepath1 = f"{FEFF_TEST_DIR}/PARAMETERS"
        parameters1 = Tags.from_file(filepath1)
        filepath2 = f"{FEFF_TEST_DIR}/PARAMETERS.2"
        parameters2 = Tags.from_file(filepath2)
        assert Tags(parameters1).diff(parameters2) == {
            "Different": {},
            "Same": {
                "CONTROL": [1, 1, 1, 1, 1, 1],
                "MPSE": [2],
                "OPCONS": "",
                "SCF": [6.0, 0, 30, 0.2, 1],
                "EXCHANGE": [0, 0.0, 0.0, 2],
                "S02": [0.0],
                "COREHOLE": "Fsr",
                "FMS": [8.5, 0],
                "XANES": [3.7, 0.04, 0.1],
                "EDGE": "K",
                "PRINT": [1, 0, 0, 0, 0, 0],
                "LDOS": [-30.0, 15.0, 0.1],
            },
        }

    def test_as_dict_and_from_dict(self):
        file_name = f"{FEFF_TEST_DIR}/PARAMETERS"
        tags = Tags.from_file(file_name)
        dct = tags.as_dict()
        tags2 = Tags.from_dict(dct)
        assert tags == tags2, "Parameters do not match to and from dict"

    def test_eels_tags(self):
        ans_1 = {
            "CONTROL": [1, 1, 1, 1, 1, 1],
            "COREHOLE": "Fsr",
            "EDGE": "K",
            "ELNES": {
                "ANGLES": "7.6 6.4",
                "BEAM_ENERGY": "200 1 0 1",
                "ENERGY": "4.0 .04 0.1",
                "MESH": "50 1",
                "POSITION": "0 0",
            },
            "EXCHANGE": [0, 0.0, 0.0, 2],
            "FMS": [7.5],
            "PRINT": [1, 0, 0, 0, 0, 0],
            "RPATH": [-1],
            "S02": [0.0],
            "SCF": [6, 0, 30, 0.2, 5],
        }
        tags_1 = Tags.from_file(f"{FEFF_TEST_DIR}/feff_eels_powder.inp")
        assert dict(tags_1) == ans_1
        ans_1["ELNES"]["BEAM_ENERGY"] = "200 0 1 1"
        ans_1["ELNES"]["BEAM_DIRECTION"] = "1 0 0"
        tags_2 = Tags.from_file(f"{FEFF_TEST_DIR}/feff_eels_x.inp")
        assert dict(tags_2) == ans_1


class TestFeffPot(unittest.TestCase):
    def test_init(self):
        filepath = f"{FEFF_TEST_DIR}/POTENTIALS"
        feff_pot = Potential.pot_string_from_file(filepath)
        dct, dr = Potential.pot_dict_from_str(feff_pot)
        assert dct["Co"] == 1, "Wrong symbols read in for Potential"
        assert dr == {0: "O", 1: "Co", 2: "O"}

    def test_single_absorbing_atom(self):
        """
        When there is only one absorbing atom in the structure, it should not appear
        in the pot_dict to avoid an error.
        """
        # one Zn+2, 9 triflate, plus water
        xyz = f"{FEFF_TEST_DIR}/feff_radial_shell.xyz"
        mol = Molecule.from_file(xyz)
        mol.set_charge_and_spin(-7)
        pot = Potential(mol, "Zn")
        # Zn should not appear in the pot_dict
        assert not pot.pot_dict.get("Zn", False)
        # Zn should only appear in the first row of the string representation
        assert str(pot).count("Zn") == 1

    def test_as_dict_and_from_dict(self):
        file_name = f"{FEFF_TEST_DIR}/HEADER"
        header = Header.from_file(file_name)
        struct = header.struct
        pot = Potential(struct, "O")
        dct = pot.as_dict()
        pot2 = Potential.from_dict(dct)
        assert str(pot) == str(pot2), "Potential to and from dict does not match"


class TestPaths(unittest.TestCase):
    def setUp(self):
        feo = Structure.from_dict(
            {
                "lattice": {
                    "a": 3.3960486211791285,
                    "alpha": 91.45136142952781,
                    "b": 3.410591877060444,
                    "beta": 89.27127081348024,
                    "c": 10.71766796897646,
                    "gamma": 120.14175587658389,
                    "matrix": [
                        [3.39597035, -0.00828486, 0.02151698],
                        [-1.70515997, 2.9534242, -0.04303398],
                        [0.06812465, -0.11799566, 10.71680189],
                    ],
                    "volume": 107.31813123502585,
                },
                "sites": [
                    {
                        "abc": [0.33497754, 0.66579918, 0.97174225],
                        "label": "Fe",
                        "properties": {
                            "coordination_no": 4,
                            "forces": [-0.01537896, -0.08731049, 0.04884326],
                        },
                        "species": [{"element": "Fe", "occu": 1}],
                        "xyz": [
                            0.06847928463257683,
                            1.8489508003914767,
                            10.392524897825345,
                        ],
                    },
                    {
                        "abc": [0.99661905, 0.00734083, 0.22366433],
                        "label": "Fe",
                        "properties": {
                            "coordination_no": 4,
                            "forces": [-0.01685376, -0.01008504, 0.05451912],
                        },
                        "species": [{"element": "Fe", "occu": 1}],
                        "xyz": [
                            3.387208508781326,
                            -0.0129676845693048,
                            2.4180946415046494,
                        ],
                    },
                    {
                        "abc": [0.00338095, 0.01072178, 0.72366433],
                        "label": "Fe",
                        "properties": {
                            "coordination_no": 4,
                            "forces": [0.01716078, 0.00955327, 0.05451912],
                        },
                        "species": [{"element": "Fe", "occu": 1}],
                        "xyz": [
                            0.04249863509042039,
                            -0.053751296415148794,
                            7.754978606437029,
                        ],
                    },
                    {
                        "abc": [0.66502246, 0.33082164, 0.47174225],
                        "label": "Fe",
                        "properties": {
                            "coordination_no": 4,
                            "forces": [0.08330257, -0.03033668, 0.04884326],
                        },
                        "species": [{"element": "Fe", "occu": 1}],
                        "xyz": [
                            1.7264300141777726,
                            0.9158834813430974,
                            5.055640939524896,
                        ],
                    },
                    {
                        "abc": [0.33062914, 0.66733572, 0.77744897],
                        "label": "O",
                        "properties": {
                            "coordination_no": 4,
                            "forces": [-0.07726687, -0.00523346, -0.05206924],
                        },
                        "species": [{"element": "O", "occu": 1}],
                        "xyz": [
                            0.03785603896498114,
                            1.8764506445041333,
                            8.310162619639584,
                        ],
                    },
                    {
                        "abc": [0.00312189, 0.99229908, 0.52714445],
                        "label": "O",
                        "properties": {
                            "coordination_no": 4,
                            "forces": [-0.06744419, 0.00047044, -0.05129314],
                        },
                        "species": [{"element": "O", "occu": 1}],
                        "xyz": [
                            -1.6455152924521734,
                            2.8684534947950637,
                            5.606667232944964,
                        ],
                    },
                    {
                        "abc": [0.99687811, 0.98917618, 0.02714445],
                        "label": "O",
                        "properties": {
                            "coordination_no": 4,
                            "forces": [0.03331469, 0.05864361, -0.05129314],
                        },
                        "species": [{"element": "O", "occu": 1}],
                        "xyz": [
                            1.7005140848662161,
                            2.9099949452040543,
                            0.2697833114717219,
                        ],
                    },
                    {
                        "abc": [0.66937086, 0.33670658, 0.27744897],
                        "label": "O",
                        "properties": {
                            "coordination_no": 4,
                            "forces": [0.04316575, 0.06429835, -0.05206924],
                        },
                        "species": [{"element": "O", "occu": 1}],
                        "xyz": [
                            1.7179261258365088,
                            0.9561539434765862,
                            2.9732786612521678,
                        ],
                    },
                ],
            }
        )
        atoms = Atoms(feo, 0, 10.0)
        self.paths = Paths(atoms, [[22, 16, 0], [250, 282, 250, 0]])

    def test_paths_string(self):
        lines = [
            "PATH",
            "---------------",
            "9999 3 1",
            "x y z ipot label",
            "-3.428830 -3.869922 -5.336884 1 Fe",
            "-5.133990 -0.916498 -5.379918 1 Fe",
            "0.000000 0.000000 0.000000 0 Fe",
            "9998 4 1",
            "x y z ipot label",
            "5.056157 2.964354 -2.082362 2 O",
            "8.473635 0.956940 2.742372 1 Fe",
            "5.056157 2.964354 -2.082362 2 O",
            "0.000000 0.000000 0.000000 0 Fe",
        ]
        answer = "\n".join(lines)
        assert answer == str(self.paths)
