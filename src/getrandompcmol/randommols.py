#!/usr/bin/env python

"""
This script generates random molecules from PubChem and optimizes them with xTB.
"""
from __future__ import annotations

import argparse as ap
import json
import os
import shutil
import subprocess
from multiprocessing import Pool

from numpy.random import RandomState

from .evaluate_calc import create_res_dir, eval_calc_ensemble, get_calc_ensemble
from .evaluate_conf import eval_conf_ensemble
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
    #### Hard-code the CIDs for testing ####
    # values = []
    # values.append(792349)
    # values.append(47090)
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
            print(f"\nDownloading CID {i:7d} ...")
            try:
                pgout = subprocess.run(
                    ["PubGrep", "--input", "cid", str(i), "--output", "sdf"],
                    check=True,
                    capture_output=True,
                    timeout=30,
                )
            except subprocess.TimeoutExpired as exc:
                print(f"Process timed out.\n{exc}")
                if os.path.exists(f"{i}.sdf"):
                    os.remove(f"{i}.sdf")
                continue
            except subprocess.CalledProcessError as exc:
                print("Status : FAIL", exc.returncode, exc.output)
                if os.path.exists(f"{i}.sdf"):
                    os.remove(f"{i}.sdf")
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
                    if os.path.exists(f"{i}.sdf"):
                        os.remove(f"{i}.sdf")
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
                    if os.path.exists(f"{i}.sdf"):
                        os.remove(f"{i}.sdf")
                    continue
                else:
                    print(" " * 3 + f"Unknown PubGrep error for {i:7d}.")
                    if os.path.exists(f"{i}.sdf"):
                        os.remove(f"{i}.sdf")
                    continue

            # print the first entry of the fourth line compund of interest in sdf format
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
            error = xtbopt(str(i))
            chdir(pwd)
            if error != "":
                if os.path.exists(f"{i}.sdf"):
                    os.remove(f"{i}.sdf")
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


def console_entry_point() -> int:
    """
    Entry point for the console script.
    """
    # parse arguments
    parser = ap.ArgumentParser(description="Generate random molecules from PubChem")
    parser.add_argument(
        "-n",
        "--numcomp",
        type=int,
        default=2,
        help="Number of compounds to generate",
        required=False,
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
    parser.add_argument(
        "--evalconfonly",
        # store two integers as lower and upper limit for the number of conformers
        nargs=2,
        type=int,
        help="ONLY evaluate the conformer ensemble. \
Provide the lower and upper limit for the number of conformers.",
        required=False,
        default=False,
    )
    parser.add_argument(
        "--evalconf",
        # store two integers as lower and upper limit for the number of conformers
        nargs=2,
        type=int,
        help="Evaluate the conformer ensemble. \
Provide the lower and upper limit for the number of conformers.",
        required=False,
        default=False,
    )
    parser.add_argument(
        "--evalcalconly",
        help="Only evaluate the QM calculations. \
This option is only useful if you have already generated the conformer ensemble. \
Give 'wipe' as the second argument to delete the conformer directories.",
        required=False,
        default=False,
        nargs=1,
        type=str,
    )

    args = parser.parse_args()

    # check for inconsistencies between arguments
    if (not args.opt) and args.crest:
        print(
            f"{bcolors.FAIL}You cannot leave --opt out and use --crest.{bcolors.ENDC}"
        )
        raise SystemExit(1)
    if (not args.opt) and args.evalconf:
        print(
            f"{bcolors.FAIL}You cannot leave --opt out and use --evalconf.{bcolors.ENDC}"
        )
        raise SystemExit(1)
    if (not args.crest) and args.evalconf:
        print(
            f"{bcolors.FAIL}You cannot leave --crest out and use --evalconf.{bcolors.ENDC}"
        )
        raise SystemExit(1)
    # raise an error if --evalconfonly is used with any other argument
    if args.evalconfonly and (args.opt or args.crest or args.evalconf):
        print(
            f"{bcolors.FAIL}You cannot use --evalconfonly with any other argument.{bcolors.ENDC}"
        )
        raise SystemExit(1)

    # check if dependencies are installed
    checkifinpath("PubGrep")
    checkifinpath("xtb")
    if args.opt:
        checkifinpath("mctc-convert")
    if args.crest:
        checkifinpath("crest")

    if args.evalconfonly:
        # get the list of directories from compounds.txt
        dirs = []
        try:
            with open("compounds.txt", encoding="UTF-8") as f:
                lines = f.readlines()
                for line in lines:
                    dirs.append(int(line.strip().split()[0]))
        except FileNotFoundError as e:
            print(f"{bcolors.FAIL}Error: {e}{bcolors.ENDC}")
            print(
                f"{bcolors.FAIL}File 'compounds.txt' not found \
and no compound directories provided.{bcolors.ENDC}"
            )
            raise SystemExit(1) from e

        eval_conf_ensemble(args.evalconfonly[0], args.evalconfonly[1], dirs)
        # exit the program
        return 0

    if args.evalcalconly:
        calcenergies = get_calc_ensemble()
        eval_calc_ensemble(calcenergies)
        wipe = False
        if args.evalcalconly[0] == "wipe":
            wipe = True
        create_res_dir(calcenergies, wipe)
        return 0

    main(args)
    return 0


if __name__ == "__main__":
    console_entry_point()
