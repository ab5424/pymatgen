#!/usr/bin/env python

"""
This module implements Compatibility corrections for mixing runs of different
functionals.
"""

from __future__ import division

__author__ = "Shyue Ping Ong, Anubhav Jain, Sai Jayaraman"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "1.0"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "Mar 19, 2012"


import os
import ConfigParser
from collections import defaultdict

from pymatgen.io.vaspio_set import MITVaspInputSet, MPVaspInputSet
from pymatgen.core.periodic_table import Element
from pymatgen.analysis.structure_analyzer import oxide_type

import abc


class Correction(object):
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def correct_entry(self, entry):
        """
        Process a single entry.

        Args:
            entry:
                A ComputedEntry object.

        Returns:
            An processed entry. None if entry is not compatible within the
            processing scheme.
        """
        return


class PotcarCorrection(Correction):
    """
    Check that POTCARs are valid.
    """
    def __init__(self, name):
        """
        Args:
            name:
                Name of POTCAR settings to use.
        """
        if name == "MP":
            input_set = MPVaspInputSet()
        elif name == "MIT":
            input_set = MITVaspInputSet()
        self.valid_potcars = set(input_set.potcar_settings.values())

    def correct_entry(self, entry):
        try:
            psp_settings = set([sym.split(" ")[1]
                                for sym
                                in entry.parameters["potcar_symbols"]])
        except KeyError:
            raise ValueError("PotcarCorrection can only be "
                             "checked for entries with a \"potcar_symbols\" in "
                             "entry.parameters")
        if not self.valid_potcars.issuperset(psp_settings):
            return None
        return entry


class GasCorrection(Correction):
    """
    Correct gas energies.
    """
    def __init__(self, name):
        """
        Args:
            name:
                Name of POTCAR settings to use.
        """
        module_dir = os.path.dirname(os.path.abspath(__file__))
        config = ConfigParser.SafeConfigParser()
        config.optionxform = str
        config.readfp(open(os.path.join(module_dir, "Compatibility.cfg")))
        cpd_energies = dict(
            config.items("{}CompoundEnergies".format(name)))
        self.cpd_energies = {k: float(v) for k, v in cpd_energies.items()}
        self.oxide_correction = {
            k: float(v) for k, v
            in config.items("{}OxideCorrection".format(name))}

    def correct_entry(self, entry):
        comp = entry.composition
        rform = entry.composition.reduced_formula
        if rform in self.cpd_energies:
            entry.entry_id = -comp.keys()[0].Z
            entry.correction += self.cpd_energies[rform] * comp.num_atoms \
                - entry.uncorrected_energy
            return entry

        correction = 0
        #Check for peroxide, superoxide, and ozonide corrections. For
        # now, ozonides will be ignored.
        if len(comp) >= 2 and Element("O") in comp:
            if "oxide_type" in entry.data:
                if entry.data["oxide_type"] in self.oxide_correction:
                    ox_corr = self.oxide_correction[
                        entry.data["oxide_type"]]
                    if entry.data["oxide_type"] is "ozonide":
                        correction += ox_corr * comp["O"] * 2 / 3
                    else:
                        correction += ox_corr * comp["O"] / 2
            elif hasattr(entry, "structure"):
                ox_type, nbonds = oxide_type(entry.structure, 1.1,
                                             return_nbonds=True)
                if ox_type in self.oxide_correction:
                    if ox_type is "ozonide":
                        correction += self.oxide_correction[ox_type] \
                            * nbonds * 2 / 3
                    else:
                        correction += self.oxide_correction[ox_type] * \
                            nbonds
            else:
                if rform in UCorrection.common_peroxides:
                    correction += self.oxide_correction["peroxide"] * \
                        comp["O"] / 2
                elif rform in UCorrection.common_superoxides:
                    correction += self.oxide_correction["superoxide"] * \
                        comp["O"] / 2
                elif rform in UCorrection.ozonides:
                    correction += self.oxide_correction["ozonide"] * \
                        comp["O"] / 2 * 3

        entry.correction += correction
        return entry

class UCorrection(Correction):
    """
    Correct U.
    """
    common_peroxides = ["Li2O2", "Na2O2", "K2O2", "Cs2O2", "Rb2O2", "BeO2",
                        "MgO2", "CaO2", "SrO2", "BaO2"]
    common_superoxides = ["LiO2", "NaO2", "KO2", "RbO2", "CsO2"]
    ozonides = ["LiO3", "NaO3", "KO3", "NaO5"]

    def __init__(self, name, compat_type):
        """
        Args:
            name:
                Name of POTCAR settings to use.
            compat_type:
                Two options, GGA or Advanced.  GGA means all GGA+U entries are
                excluded.  Advanced means mixing scheme is implemented to make
                entries compatible with each other, but entries which are
                supposed to be done in GGA+U will have the equivalent GGA
                entries excluded. For example, Fe oxides should have a U value
                under the Advanced scheme. A GGA Fe oxide run will therefore be
                excluded under the scheme.
        """
        module_dir = os.path.dirname(os.path.abspath(__file__))
        config = ConfigParser.SafeConfigParser()
        config.optionxform = str
        config.readfp(open(os.path.join(module_dir, "Compatibility.cfg")))
        if name == "MP":
            self.input_set = MPVaspInputSet()
        elif name == "MIT":
            self.input_set = MITVaspInputSet()
        else:
            raise ValueError("Invalid input set name {}".format(name))

        u_corrections = {}
        for el in self.input_set.incar_settings["LDAUU"].keys():
            sect_name = "{}{}UCorrections{}".format(name, compat_type, el)
            if sect_name in config.sections():
                corr = dict(config.items(sect_name))
                u_corrections[el] = {k: float(v) for k, v in corr.items()}

        self.u_corrections = u_corrections
        self.u_settings = self.input_set.incar_settings["LDAUU"]

        if compat_type == "GGA":
            self.u_corrections = {}
            self.u_settings = {}

    def correct_entry(self, entry):
        if entry.parameters.get("run_type", "GGA") == "HF":
            return None

        calc_u = entry.parameters.get("hubbards", None)
        calc_u = defaultdict(int) if calc_u is None else calc_u
        comp = entry.composition

        elements = sorted([el for el in comp.elements if comp[el] > 0],
                           key=lambda el: el.X)
        most_electroneg = elements[-1].symbol
        correction = 0

        ucorr = self.u_corrections.get(most_electroneg, {})
        usettings = self.u_settings.get(most_electroneg, {})

        for el in comp.elements:
            sym = el.symbol
            #Check for bad U values
            if calc_u.get(sym, 0) != usettings.get(sym, 0):
                return None
            if sym in ucorr:
                correction += float(ucorr[sym]) * comp[el]

        entry.correction += correction
        return entry


class Compatibility(object):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. This is a base class from which other specific compatibility
    schemes are implemented.

    For compatibility to be checked, the entry supplied have two additional
    restrictions in terms of its parameters key:

    1. Entry.parameters must contain a "hubbards" key which is a dict of all
       non-zero Hubbard U values used in the calculation. For example,
       if you ran a Fe2O3 calculation with Materials Project parameters,
       this would look like entry.parameters["hubbards"] = {"Fe": 5.3}
       If the "hubbards" key is missing, a GGA run is assumed.
    2. Entry.parameters must contain a "potcar_symbols" key that is a list of
       all POTCARs used in the run. Again, using the example of an Fe2O3 run
       using Materials Project parameters, this would look like
       entry.parameters["potcar_symbols"] = ['PAW_PBE Fe_pv 06Sep2000',
       'PAW_PBE O 08Apr2002'].

    It should be noted that ComputedEntries assimilated using the
    pymatgen.apps.borg package and obtained via the MaterialsProject REST
    interface using the pymatgen.matproj.rest package will automatically have
    these fields populated.
    """
    def __init__(self, corrections):
        """
        Args:
            corrections:
                List of corrections to apply.
        """
        self.corrections = corrections

    def process_entry(self, entry):
        """
        Process a single entry with the chosen Compatibility scheme.

        Args:
            entry:
                A ComputedEntry object.

        Returns:
            An adjusted entry if entry is compatible, otherwise None is
            returned.

        Raises:
            ValueError if entry do not contain "potcar_symbols" key.
        """
        if entry.parameters.get("run_type", "GGA") == "HF":
            return None

        entry.correction = 0
        for c in self.corrections:
            entry = c.correct_entry(entry)
            if entry is None:
                return None
        return entry

    def process_entries(self, entries):
        """
        Process a sequence of entries with the chosen Compatibility scheme.

        Args:
            entries - A sequence of entries.

        Returns:
            An list of adjusted entries.  Entries in the original list which
            are not compatible are excluded.
        """
        return filter(None, map(self.process_entry, entries))


class MaterialsProjectCompatibility(Compatibility):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MaterialsProject parameters (see pymatgen.io.vaspio_set.MPVaspInputSet).
    Using this compatibility scheme on runs with different parameters is not
    valid.
    """

    def __init__(self, compat_type="Advanced"):
        """
        Args:
            compat_type:
                Two options, GGA or Advanced.  GGA means all GGA+U entries are
                excluded.  Advanced means mixing scheme is implemented to make
                entries compatible with each other, but entries which are
                supposed to be done in GGA+U will have the equivalent GGA
                entries excluded. For example, Fe oxides should have a U value
                under the Advanced scheme. A GGA Fe oxide run will therefore be
                excluded under the scheme.
        """
        Compatibility.__init__(self,
                               [PotcarCorrection("MP"), GasCorrection("MP"),
                                UCorrection("MP", compat_type)])


class MITCompatibility(Compatibility):
    """
    This class implements the GGA/GGA+U mixing scheme, which allows mixing of
    entries. Note that this should only be used for VASP calculations using the
    MIT parameters (see pymatgen.io.vaspio_set MITVaspInputSet). Using
    this compatibility scheme on runs with different parameters is not valid.
    """

    def __init__(self, compat_type="Advanced"):
        """
        Args:
            compat_type:
                Two options, GGA or Advanced.  GGA means all GGA+U entries are
                excluded.  Advanced means mixing scheme is implemented to make
                entries compatible with each other, but entries which are
                supposed to be done in GGA+U will have the equivalent GGA
                entries excluded. For example, Fe oxides should have a U value
                under the Advanced scheme. A GGA Fe oxide run will therefore be
                excluded under the scheme.
        """
        Compatibility.__init__(self,
                               [PotcarCorrection("MIT"),
                                GasCorrection("MIT"),
                                UCorrection("MIT", compat_type)])
