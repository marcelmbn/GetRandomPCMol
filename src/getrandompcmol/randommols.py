#!/usr/bin/env python

"""
This script generates random molecules from PubChem and optimizes them with xTB.
"""
from __future__ import annotations

import argparse as ap
import os
import subprocess

from numpy.random import RandomState

from .miscelleanous import chdir, create_directory


def main(arguments: ap.Namespace) -> None:
    """
    Main function of the script.
    """
    # get number of compounds from command line
    numcomp = arguments.numcomp

    ### GENERAL PARAMETERS ###
    maxcid = 10000000
    maxnumat = 50

    # set the seed
    print(f"Generating random numbers between 1 and {maxcid:d} ...")
    seed = RandomState(2009)
    values = seed.randint(1, maxcid, size=3 * numcomp)
    print("Done.")

    pwd = os.getcwd()

    # run PubGrep for each value and set up a list with successful downloads
    comp: list[int] = []
    molname: list[str] = []
    for i in values:
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
                    + f"""{bcolors.WARNING}xTB error in conversion process -
 skipping CID %7d {bcolors.ENDC}"""
                    % i
                )
                continue
            elif not "normal termination" in pgout.stderr.decode("utf-8"):
                print(
                    " " * 3
                    + f"""{bcolors.WARNING}Unknown PubGrep/xTB conversion error -
 skipping CID %7d{bcolors.ENDC}"""
                    % i
                )
                errmess = "PubGrep_error" + str(i) + ".err"
                with open(errmess, "w", encoding="UTF-8") as f:
                    f.write(pgout.stderr.decode("utf-8"))
                continue
            else:
                print(" " * 3 + f"Downloaded {i:7d} successfully after xTB conversion.")

        # print the first entry of the fourth line compund of interest in sdf format
        with open(f"{i}.sdf", encoding="UTF-8") as f:
            lines = f.readlines()
            nat = int(lines[3].split()[0])
            print(" " * 3 + f"# of atoms: {nat:8d}")
            if int(nat) > maxnumat:
                print(
                    f"{bcolors.WARNING}Number of \
atoms in {i}.sdf is larger than {maxnumat} -\
skipping CID {i:7d}{bcolors.ENDC}"
                )
                continue

        direxists = create_directory(str(i))
        chdir(str(i))
        try:
            pgout = subprocess.run(
                ["xtb", f"{i}.sdf", "--opt"],
                check=True,
                capture_output=True,
                timeout=120,
            )
            with open("xtb.out", "w", encoding="UTF-8") as f:
                f.write(pgout.stdout.decode("utf-8"))
        except subprocess.TimeoutExpired as exc:
            print(" " * 3 + f"Process timed out.\n{exc}")
            chdir(pwd)
            continue
        except subprocess.CalledProcessError as exc:
            print(
                " " * 3 + f"{bcolors.FAIL}Status : FAIL{bcolors.ENDC}", exc.returncode
            )
            with open("xtb_error.out", "w", encoding="UTF-8") as f:
                f.write(exc.output.decode("utf-8"))
            with open("xtb_error.err", "w", encoding="UTF-8") as f:
                f.write(pgout.stderr.decode("utf-8"))
            print(
                f"{bcolors.WARNING}xTB optimization failed - skipping CID %7d{bcolors.ENDC}"
                % i
            )
            chdir(pwd)
            continue

        # load fourth entry of a line with ":: total charge" of xtb.out into a variable
        chrg = 0
        with open("xtb.out", encoding="UTF-8") as f:
            lines = f.readlines()
            for line in lines:
                if ":: total charge" in line:
                    chrg = round(float(line.split()[3]))
        print(" " * 3 + f"Total charge: {chrg:6d}")
        # write chrg to a file called .CHRG
        with open(".CHRG", "w", encoding="UTF-8") as f:
            f.write(str(chrg) + "\n")

        try:
            pgout = subprocess.run(
                ["mctc-convert", "xtbopt.sdf", "struc.xyz"],
                check=True,
                capture_output=True,
                timeout=120,
            )
        except subprocess.TimeoutExpired as exc:
            print(" " * 3 + f"Process timed out.\n{exc}")
            chdir(pwd)
            continue
        except subprocess.CalledProcessError as exc:
            print(" " * 3 + "Status : FAIL", exc.returncode, exc.output)
            # write the error output to a file
            with open("mctc-convert_error.err", "w", encoding="UTF-8") as f:
                f.write(pgout.stderr.decode("utf-8"))
            print(
                f"{bcolors.WARNING}mctc-convert failed - skipping CID %7d{bcolors.ENDC}"
                % i
            )
            chdir(pwd)
            continue

        chdir(pwd)

        # grep the name of the molecule from found.results (first entry in first line)
        with open("found.results", encoding="UTF-8") as f:
            first_line = f.readline()
            molname.append(first_line.split()[0])
            print(" " * 3 + f"Compound name: {molname[-1]:s}")

        print(
            f"{bcolors.OKGREEN}Structure of\
 %7d successfully generated and optimized.{bcolors.ENDC}"
            % i
        )
        comp.append(i)
        print("[" + str(len(comp)) + "/" + str(numcomp) + "]")
        if len(comp) >= numcomp:
            break

    # print number of successful downloads
    print("\nNumber of successful downloads: ", len(comp))
    print("Compounds: ", comp)

    # write the list of successful downloads to a file
    with open("compounds.txt", "w", encoding="UTF-8") as f:
        for i in comp:
            f.write(str(i) + " " + molname[comp.index(i)] + "\n")
    os.remove("found.results")


class bcolors:
    """
    Class for colorizing the output.
    """

    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


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
    args = parser.parse_args()
    main(args)
    return 0


if __name__ == "__main__":
    console_entry_point()
