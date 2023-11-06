#!/usr/bin/env python

"""
This script generates random molecules from PubChem and optimizes them with xTB.
"""
from __future__ import annotations

import argparse as ap
import json
import os
import shutil
from multiprocessing import Pool

from numpy.random import RandomState

from .evaluate_conf import eval_conf_ensemble
from .miscelleanous import bcolors, chdir, create_directory
from .qmcalc import crest_protonate, crest_sampling, get_sdf, xtb_opt, xtb_sp

HLGAP_THRESHOLD = 0.2


def main(arguments: ap.Namespace) -> None:
    """
    Main function of the script.
    """
    # get number of compounds from command line
    numcomp = arguments.numcomp
    opt = arguments.opt

    # set the seed
    print(f"Generating random numbers between 1 and {arguments.maxcid:d} ...")
    seed = RandomState(arguments.seed)
    values = seed.choice(range(1, arguments.maxcid), size=100 * numcomp, replace=False)
    print("Done.")

    pwd = os.getcwd()

    existing_dirs = False
    # run PubGrep for each value and set up a list with successful downloads
    comp: list[int] = []
    molname: list[str] = []
    #### Hard-code some CIDs for testing ####
    # values = []
    # values.append(5264)
    # values.append(3207)
    ########################################
    for i in values:
        # Check if the directory exists and skip the CID if it does
        if os.path.exists(f"{i}/{i}.sdf"):
            existing_dirs = True
            print(
                f"{bcolors.OKCYAN}\nDirectory {i} with SDF file already exists.{bcolors.ENDC}"
            )
            with open(f"{i}/{i}.sdf", encoding="UTF-8") as f:
                lines = f.readlines()
                nat = int(lines[3].split()[0])
                print(" " * 3 + f"# of atoms: {nat:8d}")
                if int(nat) > arguments.maxnumat:
                    print(
                        f"{bcolors.WARNING}Number of \
atoms in {i}.sdf is larger than {arguments.maxnumat} - \
skipping CID {i}.{bcolors.ENDC}"
                    )
                    # rm the folder
                    shutil.rmtree(f"{i}")
                    continue
        else:
            # > Run PubGrep...
            print(f"\nDownloading CID {i:7d} ...")
            pg_error = get_sdf(str(i))
            if pg_error != "":
                continue
            try:
                with open(f"{pwd}/{i}.sdf", encoding="UTF-8") as f:
                    lines = f.readlines()
                    nat = int(lines[3].split()[0])
                    print(" " * 3 + f"# of atoms: {nat:8d}")
                    if int(nat) > arguments.maxnumat:
                        print(
                            f"{bcolors.WARNING}Number of \
atoms in {i}.sdf is larger than {arguments.maxnumat} - \
skipping CID {i}.{bcolors.ENDC}"
                        )
                        # rm the sdf file
                        os.remove(f"{i}.sdf")
                        continue
            except FileNotFoundError:
                print(
                    f"{bcolors.WARNING}File {i}.sdf not found - \
skipping CID {i}.{bcolors.ENDC}"
                )
                if os.path.exists(f"{i}.sdf"):
                    os.remove(f"{i}.sdf")
                continue

            direxists = create_directory(str(i))
            shutil.copy2(f"{pwd}/{i}.sdf", f"{pwd}/{i}/{i}.sdf")

        if opt:
            chdir(str(i))
            # copy the sdf file to the new directory
            # run xTB optimization
            print(" " * 3 + f"Running xTB optimization for CID {i} ...")
            error = xtb_opt(str(i))
            chdir(pwd)
            if error != "":
                if os.path.exists(f"{i}.sdf"):
                    os.remove(f"{i}.sdf")
                continue
            hlgap = 0.0
            with open(f"{str(i)}/xtb.out", encoding="UTF-8") as f:
                lines = f.readlines()
                for line in lines:
                    if "HOMO-LUMO GAP" in line:
                        hlgap = float(line.split()[3])
                        break
            print(" " * 3 + f"HOMO-LUMO gap: {hlgap:5.2f}")
            if hlgap < HLGAP_THRESHOLD:
                print(
                    f"{bcolors.WARNING} HOMO-LUMO gap with GFN2-xTB below 0.2 eV - \
skipping CID {i}.{bcolors.ENDC}"
                )
                continue
            print(
                f"{bcolors.OKGREEN}Structure of \
{i} successfully generated and optimized.{bcolors.ENDC}"
            )
        else:
            print(
                f"{bcolors.OKGREEN}Structure of \
{i} successfully generated.{bcolors.ENDC}"
            )

        # grep the name of the molecule from found.results (first entry in first line)
        if os.path.exists("found.results"):
            with open("found.results", encoding="UTF-8") as f:
                first_line = f.readline()
                molname.append(first_line.split()[0])
                print(f"Compound name: {bcolors.BOLD}{molname[-1]:s}{bcolors.ENDC}")

        # > remove the sdf file from the main directory if it exists
        if os.path.exists(f"{i}.sdf"):
            os.remove(f"{i}.sdf")

        # > append the CID to the list of successful downloads
        comp.append(int(i))
        print("[" + str(len(comp)) + "/" + str(numcomp) + "]")

        # > break the loop if the number of successful downloads is equal to numcomp
        if len(comp) >= numcomp:
            break

    # clean up
    files_to_remove = [
        "iupac",
        "list.tmp",
        "pubchem_data",
        "found.results",
        "xtb_3d.out",
        ".sccnotconverged",
    ]
    for i in files_to_remove:
        if os.path.exists(i):
            os.remove(i)

    # print number of successful downloads
    print("\nNumber of successful downloads: ", len(comp))
    print("Compounds: ", comp)

    if len(comp) == 0:
        print(
            f"{bcolors.BOLD}No compounds were downloaded. \
Exiting the program.{bcolors.ENDC}"
        )
        exit(0)

    # write the list of successful downloads to a file
    if not existing_dirs:
        with open("compounds.txt", "w", encoding="UTF-8") as f:
            for i in comp:
                f.write(str(i) + " " + molname[comp.index(i)] + "\n")
    else:
        print(
            f"{bcolors.BOLD}Some directories already existed. \
The list of successful downloads was written only with CIDs.{bcolors.ENDC}"
        )
        with open("compounds.txt", "w", encoding="UTF-8") as f:
            for i in comp:
                f.write(str(i) + "\n")

    print("")
    print(f"{bcolors.HEADER}### CREST CONFORMER SEARCH ###{bcolors.ENDC}")

    if arguments.crest:
        # select between CREST run modes
        if arguments.crest == "protonate":
            # determine the number of cores for protonation
            protonated_cids: list[int] = []
            totalcores = os.cpu_count()
            if totalcores is None:
                n_threads = 4
            else:
                n_threads = min(int(totalcores), 4)
            crest_protonate_options: dict[str, int | float | str] = {
                "nthreads": n_threads,
            }

            print(f"{bcolors.BOLD}Running CREST protonations ...{bcolors.ENDC}")
            for i in comp:
                # run CREST protonation
                error = crest_protonate(pwd, str(i), crest_protonate_options)
                if error != "":
                    print(
                        f"{bcolors.FAIL}CREST protonation failed with error statement - \
skipping CID {i}.{bcolors.ENDC}"
                    )
                    print(" " * 3 + "ERROR: ", error)
                    continue
                # quick xtb calculation to check if HL gap is still reasonable
                chdir(str(i))
                error = xtb_sp("opt_proto.xyz")
                hlgap = 0.0
                try:
                    with open("xtb.out", encoding="UTF-8") as f:
                        lines = f.readlines()
                        for line in lines:
                            if "HOMO-LUMO GAP" in line:
                                hlgap = float(line.split()[3])
                                break
                except FileNotFoundError:
                    print(
                        f"{bcolors.FAIL} GFN2-xTB single-point failed - \
skipping protonated CID {i}.{bcolors.ENDC}"
                    )
                    continue
                chdir(pwd)
                if hlgap < HLGAP_THRESHOLD:
                    print(
                        f"{bcolors.WARNING}   HOMO-LUMO gap \
with GFN2-xTB ({hlgap}) below {HLGAP_THRESHOLD:4.2f} eV - \
skipping protonated CID {i}.{bcolors.ENDC}"
                    )
                    continue
                print(
                    f"CREST protonation for \
{bcolors.OKGREEN}CID {i}{bcolors.ENDC} successfully finished."
                )
                print(" " * 3 + f"HOMO-LUMO gap: {hlgap:5.2f}")
                protonated_cids.append(i)

            if len(protonated_cids) == 0:
                print(
                    f"{bcolors.BOLD}No protonations were successful. \
Exiting the program.{bcolors.ENDC}"
                )
            reductions = len(comp) - len(protonated_cids)
            if reductions > 0:
                print(
                    f"{bcolors.BOLD}CREST protonations were not successful for \
{reductions} compounds. Post-processing only with remaining molecules.{bcolors.ENDC}"
                )
            # further processing only for protonated CIDs
            comp = protonated_cids
            print(f"{bcolors.OKBLUE}CREST protonations finished.{bcolors.ENDC}\n")

        # get number of cores
        totalcores = os.cpu_count()
        if totalcores is None:
            totalcores = 4
        else:
            totalcores = int(totalcores)
        num_cores = min(totalcores, len(comp))
        n_threads = totalcores // num_cores

        crest_options: dict[str, int | float | str] = {
            "nthreads": n_threads,
            "mddump": 250,
            "mdlen": "x0.75",
        }
        if arguments.crest == "protonate":
            crest_options["strucfile"] = "opt_proto.xyz"
        elif arguments.crest == "normal":
            crest_options["strucfile"] = "opt.xyz"
        else:
            raise ValueError(
                f"{bcolors.FAIL}Unknown CREST run mode: {arguments.crest}{bcolors.ENDC}"
            )
        sum_cores = num_cores * n_threads
        print(f"Number of detected cores on this machine: {totalcores}")
        print(3 * " " + f"Number of cores per process: {n_threads}")
        print(3 * " " + f"Number of parallel processes: {num_cores}")
        print(
            f"{bcolors.BOLD}Running CREST sampling with {sum_cores} cores.{bcolors.ENDC}"
        )
        print("Finished conformer search for compounds: ", end="", flush=True)
        results: list[dict[str, str | int | float | list[float]]] = []
        with Pool(processes=num_cores) as p:
            results = p.starmap(
                crest_sampling,
                zip([pwd] * len(comp), comp, [crest_options] * len(comp)),
            )
        print(f"{bcolors.OKBLUE}done.{bcolors.ENDC}")
        for i, o in zip(comp, results):
            print(f"Number of conformers for CID {i}: {o['nconf']}")
            # get first element of the list of energies
            if isinstance(o["energies"], list):
                if len(o["energies"]) > 1:
                    energy_range = max(o["energies"]) - min(o["energies"])
                    print(
                        3 * " "
                        + f"Energy range of conformers: {energy_range:.3f} kcal/mol"
                    )
                    # Add the "energy_range" and "mean_energy" to the "o" dictionary
                    o["energy_range"] = energy_range
                else:
                    o["energy_range"] = 0.0

                if len(o["energies"]) > 0:
                    mean_energy = sum(o["energies"]) / len(o["energies"])
                    print(
                        3 * " "
                        + f"Mean energy of conformers:  {mean_energy:.3f} kcal/mol"
                    )
                    o["mean_energy"] = mean_energy
                else:
                    o["mean_energy"] = 0.0

                # Write the entire "o" dictionary to a JSON file in the directory
                with open(f"{i}/conformer.json", "w", encoding="UTF-8") as f:
                    json.dump(o, f, indent=4)

        if arguments.evalconf:
            eval_conf_ensemble(arguments.evalconf[0], arguments.evalconf[1], comp)

    else:
        print(f"{bcolors.BOLD}CREST conformer search was not requested.{bcolors.ENDC}")
