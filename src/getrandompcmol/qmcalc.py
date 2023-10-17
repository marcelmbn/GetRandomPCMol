"""
Module for the QM calculation of the random PCMOLs.
"""

from __future__ import annotations

import subprocess

from .miscelleanous import bcolors


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
