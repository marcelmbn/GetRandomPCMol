"""
Module for the QM calculation of the random PCMOLs.
"""

from __future__ import annotations

import os
import shutil
import subprocess

from .miscelleanous import bcolors, chdir, create_directory


def xtbopt(name: str) -> str:
    """
    Function to run xTB optimization and convert the output to xyz format.
    """
    error = ""
    pgout = None
    try:
        pgout = subprocess.run(
            ["xtb", f"{name}.sdf", "--opt"],
            check=True,
            capture_output=True,
            timeout=120,
        )
        with open("xtb.out", "w", encoding="UTF-8") as f:
            f.write(pgout.stdout.decode("utf-8"))
        with open("xtb.err", "w", encoding="UTF-8") as f:
            f.write(pgout.stderr.decode("utf-8"))
    except subprocess.TimeoutExpired as exc:
        error = " " * 3 + f"Process timed out.\n{exc}"
        print(error)
        return error
    except subprocess.CalledProcessError as exc:
        print(" " * 3 + f"{bcolors.FAIL}Status : FAIL{bcolors.ENDC}", exc.returncode)
        with open("xtb_error.out", "w", encoding="UTF-8") as f:
            f.write(exc.output.decode("utf-8"))
        error = f"{bcolors.WARNING}xTB optimization failed - skipping CID {name}.{bcolors.ENDC}"
        print(error)
        return error

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
            ["mctc-convert", "xtbopt.sdf", "opt.xyz"],
            check=True,
            capture_output=True,
            timeout=120,
        )
    except subprocess.TimeoutExpired as exc:
        error = " " * 3 + f"Process timed out.\n{exc}"
        print(error)
        return error
    except subprocess.CalledProcessError as exc:
        print(" " * 3 + "Status : FAIL", exc.returncode, exc.output)
        # write the error output to a file
        with open("mctc-convert_error.err", "w", encoding="UTF-8") as f:
            f.write(pgout.stderr.decode("utf-8"))
        error = (
            f"{bcolors.WARNING}mctc-convert failed - skipping CID {name}.{bcolors.ENDC}"
        )
        print(error)
        return error

    return error


def crest_sampling(
    homedir: str, name: str, crestsettings: dict[str, int | float | str]
) -> dict[str, int]:
    """
    Function to run CREST sampling.
    """
    error = ""
    pgout = None
    chdir(name)
    direxist = create_directory("crest")
    shutil.copy2("opt.xyz", "crest/")
    # if exist, copy the .CHRG file to the crest directory
    if os.path.exists(".CHRG"):
        shutil.copy2(".CHRG", "crest/")

    chdir("crest")

    print(f"\nRunning CREST sampling for CID {name} ...")
    error = ""
    try:
        pgout = subprocess.run(
            [
                "crest",
                "opt.xyz",
                "--squick",
                "--T",
                str(crestsettings["nthreads"]),
                "--mddump",
                str(crestsettings["mddump"]),
                "--mdlen",
                str(crestsettings["mdlen"]),
            ],
            check=True,
            capture_output=True,
            timeout=120,
        )
        with open("crest.out", "w", encoding="UTF-8") as f:
            f.write(pgout.stdout.decode("utf-8"))
        with open("crest.err", "w", encoding="UTF-8") as f:
            f.write(pgout.stderr.decode("utf-8"))
    except subprocess.TimeoutExpired as exc:
        error = " " * 3 + f"Process timed out.\n{exc}"
    except subprocess.CalledProcessError as exc:
        print(" " * 3 + f"{bcolors.FAIL}Status : FAIL{bcolors.ENDC}", exc.returncode)
        with open("crest_error.out", "w", encoding="UTF-8") as f:
            f.write(exc.output.decode("utf-8"))
        error = f"{bcolors.WARNING}CREST conformer search failed - \
skipping CID {name}.{bcolors.ENDC}"

    conformer_prop: dict[str, int] = {}
    # parse crest.out and get the number of conformers
    # the relevant line is "number of unique conformers for further calc"
    try:
        with open("crest.out", encoding="UTF-8") as f:
            lines = f.readlines()
            for line in lines:
                if "number of unique conformers for further calc" in line:
                    conformer_prop["nconf"] = int(line.split()[7])
                    print(
                        f"{bcolors.BOLD}CID: {name}{bcolors.ENDC}: "
                        + f"# of conformers: {bcolors.BOLD}{conformer_prop['nconf']}{bcolors.ENDC}"
                    )
                    break
    except FileNotFoundError:
        error = f"{bcolors.FAIL}CREST conformer search failed - \
skipping CID {name}.{bcolors.ENDC}"

    print(
        f"{bcolors.OKGREEN}Conformer ensemble of \
{name} successfully generated and optimized.{bcolors.ENDC}"
    )
    chdir(homedir)

    return conformer_prop
