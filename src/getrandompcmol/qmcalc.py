"""
Module for the QM calculation of the random PCMOLs.
"""

from __future__ import annotations

import os
import shutil
import subprocess

from .miscelleanous import bcolors, chdir, create_directory


def xtb_sp(name: str) -> str:
    """
    Function to run xTB single point calculation.

    Arguments:
    name: CID of the molecule to be optimized
    """
    error = ""
    pgout = None
    try:
        pgout = subprocess.run(
            ["xtb", name],
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
    return error


def xtb_opt(name: str) -> str:
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
                break
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
) -> dict[str, str | int | float | list[float]]:
    """
    Function to run CREST sampling.
    """
    error = ""
    pgout = None
    chdir(name)
    direxist = create_directory("crest")
    shutil.copy2(str(crestsettings["strucfile"]), "crest/")
    # if exist, copy the .CHRG file to the crest directory
    if os.path.exists(".CHRG"):
        shutil.copy2(".CHRG", "crest/")

    chdir("crest")
    conformer_prop: dict[str, str | int | float | list[float]] = {}
    # initialize conformer_prop with default values
    conformer_prop["name"] = name
    # obtain molecular charge from .CHRG
    if os.path.exists(".CHRG"):
        with open(".CHRG", encoding="UTF-8") as f:
            lines = f.readlines()
            conformer_prop["charge"] = int(lines[0].strip())
    else:
        conformer_prop["charge"] = 0
    # obtain number of atoms from opt.xyz
    with open(str(crestsettings["strucfile"]), encoding="UTF-8") as f:
        lines = f.readlines()
        conformer_prop["natoms"] = int(lines[0].strip())
    conformer_prop["energies"] = []
    conformer_prop["nconf"] = 0

    error = ""
    try:
        pgout = subprocess.run(
            [
                "crest",
                str(crestsettings["strucfile"]),
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
        )
        with open("crest.out", "w", encoding="UTF-8") as f:
            f.write(pgout.stdout.decode("utf-8"))
        with open("crest.err", "w", encoding="UTF-8") as f:
            f.write(pgout.stderr.decode("utf-8"))
    except subprocess.CalledProcessError as exc:
        print(
            f"{bcolors.FAIL}Status : FAIL for {name} with code {exc.returncode}{bcolors.ENDC}, ",
            end="",
            flush=True,
        )
        with open("crest_error.out", "w", encoding="UTF-8") as f:
            f.write(exc.output.decode("utf-8"))
        chdir(homedir)
        return conformer_prop

    # parse crest.out and get the number of conformers
    # the relevant line is "number of unique conformers for further calc"
    try:
        with open("crest.out", encoding="UTF-8") as f:
            lines = f.readlines()
            for line in lines:
                if "number of unique conformers for further calc" in line:
                    conformer_prop["nconf"] = int(line.split()[7])
                    break
    except FileNotFoundError:
        print(
            f"{bcolors.FAIL}CREST conformer search failed - \
skipping CID {name}.{bcolors.ENDC}"
        )
        chdir(homedir)
        return conformer_prop
    try:
        with open("crest.energies", encoding="UTF-8") as f:
            lines = f.readlines()
            if isinstance(conformer_prop["energies"], list):
                for line in lines:
                    conformer_prop["energies"].append(float(line.split()[1]))
    except FileNotFoundError:
        print(
            f"{bcolors.FAIL}CREST conformer search failed - \
skipping CID {name}.{bcolors.ENDC}"
        )
        conformer_prop["nconf"] = 0
        chdir(homedir)
        return conformer_prop

    print(f"{name}, ", end="", flush=True)
    chdir(homedir)

    return conformer_prop


def crest_protonate(
    homedir: str, name: str, crestsettings: dict[str, int | float | str]
) -> str:
    """
    Function to run CREST sampling.
    """
    error = ""
    pgout = None
    chdir(name)
    crest_dir = "protonation"
    direxist = create_directory(crest_dir)
    shutil.copy2("opt.xyz", f"{crest_dir}/")
    # if exist, copy the .CHRG file to the crest directory
    init_charge = 0
    if os.path.exists(".CHRG"):
        shutil.copy2(".CHRG", f"{crest_dir}/")
        with open(".CHRG", encoding="UTF-8") as f:
            lines = f.readlines()
            init_charge = int(lines[0].strip())
    # obtain number of atoms from opt.xyz
    nat = 0
    with open("opt.xyz", encoding="UTF-8") as f:
        lines = f.readlines()
        nat = int(lines[0].strip())

    chdir(crest_dir)

    error = ""
    try:
        pgout = subprocess.run(
            [
                "crest",
                "opt.xyz",
                "--protonate",
                "--T",
                str(crestsettings["nthreads"]),
            ],
            check=True,
            capture_output=True,
        )
        with open("crest.out", "w", encoding="UTF-8") as f:
            f.write(pgout.stdout.decode("utf-8"))
        with open("crest.err", "w", encoding="UTF-8") as f:
            f.write(pgout.stderr.decode("utf-8"))
    except subprocess.CalledProcessError as exc:
        print(
            f"{bcolors.FAIL}Status : FAIL for {name} with code {exc.returncode}{bcolors.ENDC}, ",
            end="",
            flush=True,
        )
        with open("crest_error.out", "w", encoding="UTF-8") as f:
            f.write(exc.output.decode("utf-8"))
        chdir(homedir)
        error = f"CREST protonation failed - skipping CID {name}."
        return error

    try:
        with open("crest.out", encoding="UTF-8") as f:
            lines = f.readlines()
            for line in lines:
                if "molecular fragmentation" in line:
                    error = (
                        f"CREST protonation failed - skipping CID {name}. "
                        + "Molecule fragmented."
                    )
                    chdir(homedir)
                    return error
    except FileNotFoundError:
        error = (
            f"CREST protonation failed - skipping CID {name}. "
            + "File 'crest.out' not found."
        )
        chdir(homedir)
        return error

    # read first nat+2 lines from protonated.xyz and write it to opt_proto.xyz
    try:
        with open("protonated.xyz", encoding="UTF-8") as f:
            lines = f.readlines()
            with open("opt_proto.xyz", "w", encoding="UTF-8") as g:
                print(f"{nat+1}\n", file=g)
                for i in range(2, nat + 3):
                    g.write(lines[i])
    except FileNotFoundError:
        error = (
            f"CREST protonation failed - skipping CID {name}. "
            + "File 'protonated.xyz' not found."
        )
        chdir(homedir)
        return error

    # write new .CHRG file with initial charge + 1
    with open(".CHRG", "w", encoding="UTF-8") as g:
        g.write(str(init_charge + 1) + "\n")
    shutil.copy2("opt_proto.xyz", f"{homedir}/{name}/")
    shutil.copy2(".CHRG", f"{homedir}/{name}/")
    chdir(homedir)

    return error


def get_sdf(cid: str) -> str:
    """
    Function to download a compound from PubChem in sdf format.
    """

    error = ""
    try:
        pgout = subprocess.run(
            ["PubGrep", "--input", "cid", cid, "--output", "sdf"],
            check=True,
            capture_output=True,
            timeout=30,
        )
    except subprocess.TimeoutExpired as exc:
        print(f"Process timed out.\n{exc}")
        if os.path.exists(f"{cid}.sdf"):
            os.remove(f"{cid}.sdf")
        error = f"ERROR - PubGrep timed out for CID {cid}."
        return error
    except subprocess.CalledProcessError as exc:
        print("Status : FAIL", exc.returncode, exc.output)
        if os.path.exists(f"{cid}.sdf"):
            os.remove(f"{cid}.sdf")
        error = f"ERROR - PubGrep failed for CID {cid}."
        return error

    # print the return code
    if pgout.returncode != 0:
        print("Return code:", pgout.returncode)

    if (pgout.stderr.decode("utf-8") == "") or (
        "normal termination" in pgout.stderr.decode("utf-8")
    ):
        print(" " * 3 + f"Downloaded {cid} successfully.")
        if not os.path.exists(f"{cid}.sdf"):
            print(
                " " * 3
                + f"{bcolors.WARNING}Warning: \
File {cid}.sdf not found even though it is allocated. Skipping...{bcolors.ENDC}"
            )
            error = f"ERROR - PubGrep failed for CID {cid}."
            return error
    else:
        if "abnormal termination" in pgout.stderr.decode("utf-8"):
            print(
                " " * 3
                + f"""{bcolors.WARNING}xTB error in conversion process - \
skipping CID {cid}.{bcolors.ENDC}"""
            )
            error = f"ERROR - PubGrep xTB error in conversion process for CID {cid}."
            return error
        elif not "normal termination" in pgout.stderr.decode("utf-8"):
            print(
                " " * 3
                + f"{bcolors.WARNING}Unknown PubGrep/xTB conversion error - \
skipping CID {cid}.{bcolors.ENDC}"
            )
            error = f"Unknown PubGrep/xTB conversion error for CID {cid}."
        else:
            print(" " * 3 + f"Unknown PubGrep error for {cid:7d}.")
            return error
        if os.path.exists(f"{cid}.sdf"):
            os.remove(f"{cid}.sdf")

    return error
