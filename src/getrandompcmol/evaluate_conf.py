"""
Evaluate a list of conformer ensembles
and setup a directory structure for subsequent QM calculations.
"""

from __future__ import annotations

import json
import os

from numpy.random import RandomState

from .miscelleanous import bcolors, chdir, create_directory


def eval_conf_ensemble(minnumconf: int, maxnumconf: int, mols: list[int]) -> None:
    """
    Evaluates a given conformer ensemble.

    Arguments:
    minnumconf: minimum number of conformers desired per ensemble
    maxnumconf: maximum number of conformers desired per ensemble
    mols: list of molecules to be evaluated
    """

    print(
        f"{bcolors.BOLD}Minimum number of conformers \
desired per ensemble: {minnumconf}{bcolors.ENDC}"
    )
    print(
        f"{bcolors.BOLD}Maximum number of conformers \
desired per ensemble: {maxnumconf}{bcolors.ENDC}"
    )
    print("")

    pwd = os.getcwd()
    # go through directories and evaluate the conformer ensemble
    conf_props = []
    for i in mols:
        chdir(str(i))
        # read the conformer.json file in the directory
        try:
            with open("conformer.json", encoding="UTF-8") as f:
                conformer = json.load(f)
        except FileNotFoundError as e:
            print(f"{bcolors.FAIL}Error: {e}{bcolors.ENDC}")
            print(
                f"{bcolors.FAIL}File 'conformer.json' not found \
in directory {i}.{bcolors.ENDC}"
            )
            chdir(pwd)
            raise SystemExit(1) from e

        # check if the conformer ensemble matches to the constraints
        if conformer["nconf"] < minnumconf:
            print(
                f"{bcolors.WARNING}Warning: \
Number of conformers ({conformer['nconf']}) in {i} is less than \
the minimum number of conformers ({minnumconf}). Skipping...{bcolors.ENDC}"
            )
        if conformer["nconf"] > maxnumconf:
            print(
                f"{bcolors.OKCYAN}Number of conformers ({conformer['nconf']}) \
in {i} is greater than the maximum number of conformers ({maxnumconf}). \
Truncating the ensemble randomly to {maxnumconf}...{bcolors.ENDC}"
            )
            seed = RandomState(1995)
            values = seed.randint(1, conformer["nconf"] + 1, maxnumconf)
            # sort the values in ascending order
            values.sort()
            # take only the conformers with indices in values
            # and recalculate the respective properties
            # add an index referring to the original conformer numbering to the conformer
            conformer["indices"] = values.tolist()
            conformer["energies"] = [conformer["energies"][j - 1] for j in values]
            conformer["nconf"] = maxnumconf
            # calculate energy range
            conformer["energy_range"] = max(conformer["energies"]) - min(
                conformer["energies"]
            )
            # calculate mean energy
            conformer["mean_energy"] = sum(conformer["energies"]) / conformer["nconf"]
        else:
            conformer["indices"] = list(range(1, conformer["nconf"] + 1))

        # read the relevant structures (the index) from the crest_conformers.xyz file
        # and write them to a new file
        try:
            with open("crest/crest_conformers.xyz", encoding="UTF-8") as f:
                lines = f.readlines()
                # write the number of atoms to the file
                for j in conformer["indices"]:
                    # read the structure of the conformer index
                    # and write it to the file "conformer["index"].xyz"
                    with open(f"{j}.xyz", "w", encoding="UTF-8") as g:
                        line = (j - 1) * (conformer["natoms"] + 2)
                        g.write(lines[line])
                        for k in range(conformer["natoms"] + 1):
                            line += 1
                            g.write(lines[line])
        except FileNotFoundError as e:
            print(f"{bcolors.FAIL}Error: {e}{bcolors.ENDC}")
            print(
                f"{bcolors.FAIL}File 'crest_conformers.xyz' not found \
in directory {i}/crest.{bcolors.ENDC}"
            )
            chdir(pwd)
            raise SystemExit(1) from e

        # append the conformer properties to the list
        conf_props.append(conformer)
        # print the single properties for the conformer of interest
        print(f"{i}:")
        print(f"Number of conformers: {conformer['nconf']}")
        print(f"Energy range: {conformer['energy_range']}")
        print(f"Mean energy: {conformer['mean_energy']}")
        print(f"Number of atoms: {conformer['natoms']}")
        print(f"Charge: {conformer['charge']}")
        if "indices" in conformer:
            print(f"Conformer indices: {conformer['indices']}")
        print("")
        chdir(pwd)
