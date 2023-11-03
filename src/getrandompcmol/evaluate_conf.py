"""
Evaluate a list of conformer ensembles
and setup a directory structure for subsequent QM calculations.
"""

from __future__ import annotations

import json
import os
import shutil

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
    conf_props: list[dict[str, str | int | float | list[float]]] = []
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
            chdir(pwd)
            continue
        if conformer["nconf"] > maxnumconf:
            print(
                f"{bcolors.OKCYAN}Number of conformers ({conformer['nconf']}) \
in {i} is greater than the maximum number of conformers ({maxnumconf}). \
Truncating the ensemble randomly to {maxnumconf}...{bcolors.ENDC}"
            )
            seed = RandomState(1995)
            # create a list of 'maxnumconf' random integers between 1 and conformer["nconf"]
            # WITHOUT redundancies
            values = seed.choice(conformer["nconf"], maxnumconf, replace=False) + 1
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

        # create a subdirectory for each conformer
        # and copy the relevant structure files to it (see below)
        if os.path.exists("index.conformers"):
            os.remove("index.conformers")
        for j in conformer["indices"]:
            direxist = create_directory(str(j))
            if direxist:
                print(
                    f"{bcolors.WARNING}Warning: \
Directory {j} already exists. Skipping...{bcolors.ENDC}"
                )
            if conformer["charge"] != 0:
                # copy the .CHRG file to the new directory
                if os.path.exists(".CHRG"):
                    shutil.copy2(".CHRG", str(j))
                else:
                    print(
                        f"{bcolors.FAIL}Fail: \
File '.CHRG' not found in directory {i} even though it is allocated.{bcolors.ENDC}"
                    )
                    chdir(pwd)
                    raise SystemExit(1)
            # append the conformer index (the dir name) to a list for later use
            with open("index.conformers", "a", encoding="UTF-8") as f:
                f.write(f"{j}\n")

        # read the relevant structures (the index) from the crest_conformers.xyz file
        # and write them to a new file
        try:
            with open("crest/crest_conformers.xyz", encoding="UTF-8") as f:
                lines = f.readlines()
                # write the number of atoms to the file
                for j in conformer["indices"]:
                    # read the structure of the conformer index
                    # and write it to the file "conformer["index"].xyz"
                    with open(f"{j}/{j}.xyz", "w", encoding="UTF-8") as g:
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

        # check if the name of the molecule is part of 'conformer'. If not, append it.
        if "name" not in conformer:
            conformer["name"] = i
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

    # write the compounds with eligible conformer ensembles to a file
    with open("compounds.conformers.txt", "w", encoding="UTF-8") as f:
        for confdict in conf_props:
            # search for this entry in the list of compounds and
            # write the name of the compound to the file
            trivialname = ""
            with open("compounds.txt", encoding="UTF-8") as g:
                lines = g.readlines()
                for line in lines:
                    if str(confdict["name"]) in line:
                        try:
                            trivialname = line.split()[1]
                        except IndexError:
                            trivialname = ""

            f.write(f"{confdict['name']} {trivialname}\n")
