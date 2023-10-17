#!/usr/bin/env python

"""
This script generates random molecules from PubChem and optimizes them with xTB.
"""
from __future__ import annotations

import argparse as ap
import os
import shutil
import subprocess
from multiprocessing import Pool

from numpy.random import RandomState

from .miscelleanous import bcolors, chdir, checkifinpath, create_directory
from .qmcalc import crest_sampling, xtbopt


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
    values = seed.randint(1, arguments.maxcid, size=100 * numcomp)
    print("Done.")

    pwd = os.getcwd()

    existing_dirs = False
    # run PubGrep for each value and set up a list with successful downloads
    comp: list[int] = []
    molname: list[str] = []
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
            print(f"\nDownloading CID {i:7d} ...")
            try:
                pgout = subprocess.run(
                    ["PubGrep", "--input", "cid", str(i)],
                    check=True,
                    capture_output=True,
                    timeout=30,
                )
            except subprocess.TimeoutExpired as exc:
                print(f"Process timed out.\n{exc}")
                continue
            except subprocess.CalledProcessError as exc:
                print("Status : FAIL", exc.returncode, exc.output)
                continue

            # print the return code
            if pgout.returncode != 0:
                print("Return code:", pgout.returncode)

            if pgout.stderr.decode("utf-8") == "":
                print(" " * 3 + f"Downloaded {i} successfully.")
            else:
                if "abnormal termination" in pgout.stderr.decode("utf-8"):
                    print(
                        " " * 3
                        + f"""{bcolors.WARNING}xTB error in conversion process - \
skipping CID {i}.{bcolors.ENDC}"""
                    )
                    continue
                elif not "normal termination" in pgout.stderr.decode("utf-8"):
                    print(
                        " " * 3
                        + f"{bcolors.WARNING}Unknown PubGrep/xTB conversion error - \
skipping CID {i}.{bcolors.ENDC}"
                    )
                    errmess = "PubGrep_error" + str(i) + ".err"
                    with open(errmess, "w", encoding="UTF-8") as f:
                        f.write(pgout.stderr.decode("utf-8"))
                    continue
                else:
                    print(
                        " " * 3
                        + f"Downloaded {i:7d} successfully after xTB conversion."
                    )

            # print the first entry of the fourth line compund of interest in sdf format
            with open(f"{i}.sdf", encoding="UTF-8") as f:
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

            direxists = create_directory(str(i))
            shutil.copy2(f"{pwd}/{i}.sdf", f"{pwd}/{i}/{i}.sdf")

        if opt:
            chdir(str(i))
            # copy the sdf file to the new directory
            # run xTB optimization
            print(" " * 3 + f"Running xTB optimization for CID {i} ...")
            error = xtbopt(str(i))
            if error != "":
                continue
            chdir(pwd)
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
        comp.append(i)
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
    ]
    for i in files_to_remove:
        if os.path.exists(i):
            os.remove(i)

    # print number of successful downloads
    print("\nNumber of successful downloads: ", len(comp))
    print("Compounds: ", comp)

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
        crest_options: dict[str, int | float | str] = {
            "nthreads": 8,
            "mddump": 250,
            "mdlen": "x0.75",
        }
        num_cores = 2
        sum_cores = num_cores * crest_options["nthreads"]
        print(
            f"{bcolors.BOLD}Running CREST sampling with {sum_cores} cores.{bcolors.ENDC}"
        )
        with Pool(processes=num_cores) as p:
            p.starmap(
                crest_sampling,
                zip([pwd] * len(comp), comp, [crest_options] * len(comp)),
            )


def console_entry_point() -> int:
    # parse arguments
    parser = ap.ArgumentParser(description="Generate random molecules from PubChem")
    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        default=2,
        help="Number of compounds to generate",
        required=True,
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s {version}".format(version="0.1.0"),
    )
    parser.add_argument(
        "--opt",
        action="store_true",
        help="Optimize the generated structures with xTB",
        required=False,
        default=False,
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2009,
        help="Seed for the random number generator",
        required=False,
    )
    parser.add_argument(
        "--maxcid",
        type=int,
        default=100000,
        help="Maximum CID to generate random numbers",
        required=False,
    )
    parser.add_argument(
        "--maxnumat",
        type=int,
        default=35,
        help="Maximum number of atoms in a molecule",
        required=False,
    )
    parser.add_argument(
        "--crest",
        action="store_true",
        help="Conformer sampling with CREST",
        required=False,
        default=False,
    )

    args = parser.parse_args()

    # check for inconsistencies between arguments
    if (not args.opt) and args.crest:
        print(f"{bcolors.WARNING}You cannot use both --opt and --crest.{bcolors.ENDC}")

    # check if dependencies are installed
    checkifinpath("PubGrep")
    if args.opt:
        checkifinpath("xtb")
        checkifinpath("mctc-convert")
    if args.crest:
        checkifinpath("crest")
    main(args)
    return 0


if __name__ == "__main__":
    console_entry_point()
