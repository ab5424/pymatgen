# coding: utf-8
# Copyright (c) Pymatgen Development Team.
# Distributed under the terms of the MIT License.

"""
This module implements methods for writing LAMMPS input files.
"""


import os
import re
import shutil
import warnings
from string import Template

from monty.io import zopen
from monty.json import MSONable

from pymatgen.io.lammps.data import LammpsData

__author__ = "Kiran Mathew, Brandon Wood, Zhi Deng"
__copyright__ = "Copyright 2018, The Materials Virtual Lab"
__version__ = "1.0"
__maintainer__ = "Zhi Deng"
__email__ = "z4deng@eng.ucsd.edu"
__date__ = "Aug 1, 2018"


class GulpRun(MSONable):
    """
    Examples for various simple LAMMPS runs with given simulation box,
    force field and a few more settings. Experience LAMMPS users should
    consider using write_lammps_inputs method with more sophisticated
    templates.

    """

    template_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "templates")

    def __init__(self, structure, title=None):
        """
        Base constructor.

        Args:
            structure (Structure): String template for input script
                with placeholders. The format for placeholders has to
                be '$variable_name', e.g., '$temperature'
            filename (str): Filename for the input script.

        """
        self.structure = structure
        self.title = title

    def __str__(self):
        """
        CSSR.__str__ method is modified to padd 0's to the CSSR site data.
        The padding is to conform with the CSSR format supported Zeo++.
        The oxidation state is stripped from site.specie
        Also coordinate system is rotated from xyz to zxy
        """
        output = [
            "title",
            str(self.title) if self.title else str(self.structure.formula),
            "end",
            "cell",
            "{:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}".format(
                self.structure.lattice.c,
                self.structure.lattice.a,
                self.structure.lattice.b,
                self.structure.lattice.gamma,
                self.structure.lattice.alpha,
                self.structure.lattice.beta,
            ),
            "frac {}".format(len(self.structure))
        ]
        for site in self.structure.sites:
            # if not hasattr(site, 'charge'):
            #    charge = 0
            # else:
            #    charge = site.charge
            charge = site.charge if hasattr(site, "charge") else 0
            specie = site.specie.symbol
            # specie = site.species_string
            output.append(
                "{} core {:.4f} {:.4f} {:.4f} {:.4f}".format(
                    specie, site.c, site.a, site.b, charge
                )
            )

        return "\n".join(output)

    def write_inputs(self, filename, **kwargs):
        """
        Writes all input files (input script, and data if needed).
        Other supporting files are not handled at this moment.

        Args:
            filename (str): Directory to output the input files.
            **kwargs: kwargs supported by LammpsData.write_file.

        """
        with zopen(filename, "wt") as f:
            f.write(self.__str__())

    @classmethod
    def md(cls, data, force_field, temperature, nsteps, other_settings=None):
        r"""
        Example for a simple MD run based on template md.txt.

        Args:
            data (LammpsData or str): Data file as a LammpsData
                instance or path to an existing data file.
            force_field (str): Combined force field related cmds. For
                example, 'pair_style eam\npair_coeff * * Cu_u3.eam'.
            temperature (float): Simulation temperature.
            nsteps (int): No. of steps to run.
            other_settings (dict): other settings to be filled into
                placeholders.

        """
        template_path = os.path.join(cls.template_dir, "md.txt")
        with open(template_path) as f:
            script_template = f.read()
        settings = other_settings.copy() if other_settings is not None else {}
        settings.update({"force_field": force_field, "temperature": temperature, "nsteps": nsteps})
        script_filename = "in.md"
        return cls(
            script_template=script_template,
            settings=settings,
            data=data,
            script_filename=script_filename,
        )


