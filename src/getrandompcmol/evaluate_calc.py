"""
Evaluate the output of the QM calculations based on the conformer ensemble directory structure.
"""

from __future__ import annotations

import json
import os
import warnings

import numpy as np
from scipy.stats import spearmanr

from .miscelleanous import bcolors, chdir

# define a general parameter for the evaluation of the conformer ensemble
HARTREE2KCAL = 627.5094740631  # Hartree
MINDIFF = 0.01  # kcal/mol


def get_calc_ensemble() -> dict[int, dict[str, list[int | float]]]:
    """
    Evaluates a given conformer ensemble.

    """

    pwd = os.getcwd()
    # check if the directory structure is correct
    if not os.path.exists("compounds.conformers.txt"):
        print(
            f"{bcolors.FAIL}Error: Directory structure is not correct. \
File 'compounds.conformers.txt' not found.{bcolors.ENDC}"
        )
        raise SystemExit(1)
    try:
        with open("compounds.conformers.txt", encoding="UTF-8") as f:
            lines = f.readlines()
            mols = [int(i.strip()) for i in lines]
    except FileNotFoundError as e:
        print(f"{bcolors.FAIL}Error: {e}{bcolors.ENDC}")
        print(f"{bcolors.FAIL}File 'compounds.conformers.txt' not found.{bcolors.ENDC}")
        raise SystemExit(1) from e

    energies: dict[int, dict[str, list[int | float]]] = {}
    for i in mols:
        chdir(str(i))
        try:
            with open("index.conformers", encoding="UTF-8") as f:
                lines = f.readlines()
                conformer = [int(j.strip()) for j in lines]
        except FileNotFoundError as e:
            print(f"{bcolors.FAIL}Error: {e}{bcolors.ENDC}")
            print(
                f"{bcolors.FAIL}File 'index.conformers' not found \
in directory {i}.{bcolors.ENDC}"
            )
            chdir(pwd)
            raise SystemExit(1) from e

        # read the conformer.json file in the directory
        try:
            with open("conformer.json", encoding="UTF-8") as f:
                conformer_prop = json.load(f)
        except FileNotFoundError as e:
            print(f"{bcolors.FAIL}Error: {e}{bcolors.ENDC}")
            print(
                f"{bcolors.FAIL}File 'conformer.json' not found \
in directory {i}.{bcolors.ENDC}"
            )
            chdir(pwd)
            raise SystemExit(1) from e
        # check if the conformer ensemble matches to the constraints
        if len(conformer) == 0:
            print(
                f"{bcolors.WARNING}Warning: \
Number of conformers ({len(conformer)}) in {i} is zero. Skipping...{bcolors.ENDC}"
            )
            chdir(pwd)
            continue
        else:
            print(
                f"{bcolors.OKGREEN}Number of effective conformers \
in {i} is {len(conformer)}{bcolors.ENDC}"
            )

        # check for charge constraints
        if conformer_prop["charge"] != 0:
            print(
                f"{bcolors.WARNING}Warning: \
Charge of molecule {i} is not zero. Skipping...{bcolors.ENDC}"
            )
            chdir(pwd)
            continue

        tmpgfn2: list[float] = []
        tmpgp3: list[float] = []
        tmpwb97xd4: list[float] = []
        tmpindex: list[int] = []
        indices: list[int] = []

        moldir = os.getcwd()

        # go through the conformers and read the energies for gfn2, gp3, and TZ
        continuewithmol = False
        for j in conformer:
            chdir(str(j))
            # check if the calculation directories exist
            if not os.path.exists("gfn2"):
                print(
                    f"{bcolors.FAIL}Error: Directory structure is not correct. \
Directory 'gfn2' not found in {j}.{bcolors.ENDC}"
                )
                continuewithmol = True
                break
            if not os.path.exists("gp3"):
                print(
                    f"{bcolors.FAIL}Error: Directory structure is not correct. \
Directory 'gp3' not found in {j}.{bcolors.ENDC}"
                )
                continuewithmol = True
                break
            if not os.path.exists("TZ"):
                print(
                    f"{bcolors.FAIL}Error: Directory structure is not correct. \
Directory 'TZ' not found in {j}.{bcolors.ENDC}"
                )
                continuewithmol = True
                break
            # read the energies from the output files
            tmpgfn2.append(read_energy_file("gfn2/energy"))
            tmpgp3.append(read_energy_file("gp3/energy"))
            tmpwb97xd4.append(read_energy_file("TZ/energy"))
            tmpindex.append(j)

            chdir(moldir)

        if continuewithmol:
            chdir(pwd)
            continue

        # sort the energies in ascending order with the TZ/wB97X-D4 energies as reference
        # the GFN2 and GP3 energies should be sorted with the same indices as the TZ energies.

        energies[i] = {}
        energies[i]["GFN2"] = []
        energies[i]["GP3"] = []
        energies[i]["wB97X-D4"] = []
        energies[i]["conformer_index"] = []
        energies[i]["conformer_lowest"] = []

        indices = sorted(range(len(tmpwb97xd4)), key=lambda k: tmpwb97xd4[k])
        tmpwb97xd4 = [tmpwb97xd4[k] for k in indices]
        tmpgfn2 = [tmpgfn2[k] for k in indices]
        tmpgp3 = [tmpgp3[k] for k in indices]
        tmpindex = [tmpindex[k] for k in indices]

        # check if energy values in the reference are closer to each other than
        # the minimum difference. For that, go through all combinations of energies
        # and check if the difference is smaller than the minimum difference.

        conformer_deleted = (
            False  # Initialize a flag tracking whether a conformer was deleted
        )
        if len(tmpwb97xd4) > 2:
            while True:
                for k in range(1, len(tmpwb97xd4)):
                    conformer_deleted = False  # Initialize a flag tracking
                    # whether a conformer was deleted
                    if len(tmpwb97xd4) <= 3:
                        break
                    if (
                        abs((tmpwb97xd4[k] - tmpwb97xd4[k - 1]) * HARTREE2KCAL)
                        < MINDIFF
                    ):
                        # delete the conformer with the higher energy
                        if tmpwb97xd4[k] > tmpwb97xd4[k - 1]:
                            del tmpwb97xd4[k]
                            del tmpgfn2[k]
                            del tmpgp3[k]
                            del tmpindex[k]
                            del conformer[k]
                        else:
                            del tmpwb97xd4[k - 1]
                            del tmpgfn2[k - 1]
                            del tmpgp3[k - 1]
                            del tmpindex[k - 1]
                            del conformer[k - 1]
                        conformer_deleted = True  # Set the flag to True
                        print(
                            f"{bcolors.WARNING}Warning: \
conformer {k} was deleted in {i} due to too close-lying energies.{bcolors.ENDC}"
                        )
                        break  # Break out of the for loop
                if not conformer_deleted:
                    # If no conformers were deleted in the current iteration,
                    # break out of the while loop
                    break

        ### DEV OUTPUT ###
        # print the tmp arrays next to each other for comparison
        for j in range(len(conformer)):
            print(
                f"{tmpindex[j]:4d} {tmpwb97xd4[j]:10.6f} {tmpgfn2[j]:10.6f} \
{tmpgp3[j]:10.6f}"
            )

        for j in range(len(conformer)):
            if j == 0:
                energies[i]["conformer_lowest"].append(tmpindex[j])
            else:
                energies[i]["wB97X-D4"].append(
                    (tmpwb97xd4[j] - tmpwb97xd4[0]) * HARTREE2KCAL
                )
                energies[i]["GFN2"].append((tmpgfn2[j] - tmpgfn2[0]) * HARTREE2KCAL)
                energies[i]["GP3"].append((tmpgp3[j] - tmpgp3[0]) * HARTREE2KCAL)
                energies[i]["conformer_index"].append(tmpindex[j])

        # add the conformer_prop infos to the energies dict
        energies[i]["natoms"] = conformer_prop["natoms"]
        energies[i]["charge"] = conformer_prop["charge"]
        energies[i]["nconf"] = conformer_prop["nconf"]

        chdir(pwd)

    # write the old and new indices to a JSON file
    with open("energies.json", "w", encoding="UTF-8") as f:
        json.dump(energies, f, indent=4)
    return energies


def eval_calc_ensemble(confe: dict[int, dict[str, list[int | float]]]) -> None:
    """
    Evaluates a given conformer ensemble in dictionary format.
    """

    # total number of data points and total number of molecules
    ndatapoints = 0
    nmolecules = len(confe.keys())

    # calculate the Spearman rank correlation coefficient for the GFN2 and GP3 energies
    # the order of the wB97X-D4 energies is used as reference
    spearmanrcc: dict[str, list[float]] = {}
    spearmanrcc["GFN2"] = []
    spearmanrcc["GP3"] = []
    # 1st: iterate over the molecules via the keys of the dictionary items
    print(
        f"\n{bcolors.OKCYAN}Spearman rank correlation coefficients \
for the energy ranking in each selected conformer ensemble:{bcolors.ENDC}"
    )
    print(
        5 * " "
        + f"{bcolors.OKCYAN}Mol:"
        + 5 * " "
        + "GFN2"
        + 4 * " "
        + f"GP3{bcolors.ENDC}"
    )
    for i in confe.keys():
        # 2nd: calculate the Spearman rank correlation coefficient

        gfn2scc = 0.0
        gp3scc = 0.0
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("error")
                gfn2scc = spearmanr(confe[i]["wB97X-D4"], confe[i]["GFN2"])[0]
        except ValueError as e:
            print(f"{bcolors.FAIL}Error: {e}{bcolors.ENDC}")
            print(
                f"{bcolors.FAIL}Error in molecule {i}: \
Spearman rank correlation coefficient could not be calculated.{bcolors.ENDC}"
            )
            raise SystemExit(1) from e
        except RuntimeWarning as e:
            print(f"{bcolors.WARNING}ERROR in {i}: {e}\nSkipping...{bcolors.ENDC}")
            continue
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("error")
                gp3scc = spearmanr(confe[i]["wB97X-D4"], confe[i]["GP3"])[0]
        except ValueError as e:
            print(f"{bcolors.FAIL}Error: {e}{bcolors.ENDC}")
            print(
                f"{bcolors.FAIL}Error in molecule {i}: \
Spearman rank correlation coefficient could not be calculated.{bcolors.ENDC}"
            )
            raise SystemExit(1) from e
        except RuntimeWarning as e:
            print(f"{bcolors.WARNING}ERROR in {i}: {e}\nSkipping...{bcolors.ENDC}")
            continue

        if np.isnan(gfn2scc) or np.isnan(gp3scc):
            print("Skipping...")
            continue
        spearmanrcc["GFN2"].append(gfn2scc)
        spearmanrcc["GP3"].append(gp3scc)
        print(
            3 * " "
            + f"{i:8d}: {spearmanrcc['GFN2'][-1]:6.3f} {spearmanrcc['GP3'][-1]:6.3f}"
        )
        # 3rd: add the counter for the number of data points
        ndatapoints += len(confe[i]["wB97X-D4"])

    # calculate the mean and standard deviation of the Spearman rank correlation coefficient
    print(
        f"{bcolors.OKCYAN}       Mean and \
StdDev of the Spearman rank correlation coefficient:{bcolors.ENDC}"
    )
    print(
        f"{bcolors.BOLD}GFN2: {np.mean(spearmanrcc['GFN2']):6.3f}   \
{np.std(spearmanrcc['GFN2']):6.3f}{bcolors.ENDC}"
    )
    print(
        f"{bcolors.BOLD}GP3:  {np.mean(spearmanrcc['GP3']):6.3f}   \
{np.std(spearmanrcc['GP3']):6.3f}{bcolors.ENDC}"
    )

    # calculate the number of cases for which the Spearman rank correlation coefficient
    # is better with GP3 than with GFN2 and vice versa
    better_scc: dict[str, int] = {}
    better_scc["GFN2"] = 0
    better_scc["GP3"] = 0
    better_scc["equal"] = 0
    for i in range(len(spearmanrcc["GFN2"])):
        if spearmanrcc["GFN2"][i] > spearmanrcc["GP3"][i]:
            better_scc["GFN2"] += 1
        elif spearmanrcc["GFN2"][i] < spearmanrcc["GP3"][i]:
            better_scc["GP3"] += 1
        else:
            better_scc["equal"] += 1
    print(
        f"{bcolors.OKCYAN}       Number of cases for which \
the Spearman rank correlation coefficient is better:{bcolors.ENDC}"
    )
    print(f"{bcolors.BOLD}GFN2:  {better_scc['GFN2']:6d}{bcolors.ENDC}")
    print(f"{bcolors.BOLD}GP3:   {better_scc['GP3']:6d}{bcolors.ENDC}")
    print(f"{bcolors.BOLD}equal: {better_scc['equal']:6d}{bcolors.ENDC}")

    # calculate the RMSD of the GFN2 and GP3 energies with respect to the wB97X-D4 energies
    # for the whole data set
    rmsd: dict[str, list[float]] = {}
    rmsd["GFN2"] = []
    rmsd["GP3"] = []
    print(
        f"\n{bcolors.OKCYAN}Root mean square deviations \
of the energies in each selected conformer ensemble:{bcolors.ENDC}"
    )
    print(
        5 * " "
        + f"{bcolors.OKCYAN}Mol:"
        + 5 * " "
        + "GFN2"
        + 4 * " "
        + f"GP3{bcolors.ENDC}"
    )
    for i in confe.keys():
        rmsd["GFN2"].append(
            np.sqrt(
                np.mean(np.square(np.subtract(confe[i]["wB97X-D4"], confe[i]["GFN2"])))
            )
        )
        rmsd["GP3"].append(
            np.sqrt(
                np.mean(np.square(np.subtract(confe[i]["wB97X-D4"], confe[i]["GP3"])))
            )
        )
        print(3 * " " + f"{i:8d}: {rmsd['GFN2'][-1]:6.3f} {rmsd['GP3'][-1]:6.3f}")

    # calculate the mean and standard deviation of the RMSD
    print(
        f"{bcolors.OKCYAN}       Mean and \
StdDev of the RMSD:{bcolors.ENDC}"
    )
    print(
        f"{bcolors.BOLD}GFN2: {np.mean(rmsd['GFN2']):6.3f}   \
{np.std(rmsd['GFN2']):6.3f}{bcolors.ENDC}"
    )
    print(
        f"{bcolors.BOLD}GP3:  {np.mean(rmsd['GP3']):6.3f}   \
{np.std(rmsd['GP3']):6.3f}{bcolors.ENDC}"
    )

    # calculate the number of cases
    # for which the RMSD is better with GP3 than with GFN2 and vice versa
    better_rmsd: dict[str, int] = {}
    better_rmsd["GFN2"] = 0
    better_rmsd["GP3"] = 0
    better_rmsd["equal"] = 0
    for i in range(len(rmsd["GFN2"])):
        if rmsd["GFN2"][i] > rmsd["GP3"][i]:
            better_rmsd["GP3"] += 1
        elif rmsd["GFN2"][i] < rmsd["GP3"][i]:
            better_rmsd["GFN2"] += 1
        else:
            better_rmsd["equal"] += 1
    print(
        f"{bcolors.OKCYAN}       Number of cases for which \
the RMSD is better:{bcolors.ENDC}"
    )
    print(f"{bcolors.BOLD}GFN2:  {better_rmsd['GFN2']:6d}{bcolors.ENDC}")
    print(f"{bcolors.BOLD}GP3:   {better_rmsd['GP3']:6d}{bcolors.ENDC}")
    print(f"{bcolors.BOLD}equal: {better_rmsd['equal']:6d}{bcolors.ENDC}")

    # print the total number of data points and the total number of molecules considered
    print(
        f"\n{bcolors.BOLD}Total number of data points: \
{ndatapoints}{bcolors.ENDC}"
    )
    print(
        f"{bcolors.BOLD}Total number of molecules: \
{nmolecules}{bcolors.ENDC}"
    )


def read_energy_file(file: str) -> float:
    """
    Reads the energy from the 'energy' output file of a QM calculation.
    """
    try:
        with open(file, encoding="UTF-8") as f:
            lines = f.readlines()
            energy = float(lines[1].strip().split()[1])
        return energy
    except FileNotFoundError as e:
        print(f"{bcolors.FAIL}Error: File {file} not found.{bcolors.ENDC}")
        raise SystemExit(1) from e
    except ValueError as e:
        print(f"{bcolors.FAIL}Error: {e} not a float.{bcolors.ENDC}")
        raise SystemExit(1) from e
    # except for the case that the file is empty
    except IndexError as e:
        print(f"{bcolors.FAIL}Error: File {file} is empty.{bcolors.ENDC}")
        raise SystemExit(1) from e
    except Exception as e:
        print(f"{bcolors.FAIL}Error: {e} - general error.{bcolors.ENDC}")
        raise SystemExit(1) from e
