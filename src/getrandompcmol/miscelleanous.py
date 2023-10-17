"""
Miscelleanous functions that are used in the project.
"""

from __future__ import annotations

import os
import subprocess as sp


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


def create_directory(name: str) -> bool:
    """
    Creates a directory with the given name if it does not exist already.
    """
    if not os.path.exists(name):
        os.mkdir(name)
        exist = False
    else:
        exist = True

    # check if the new directory exists and raise an error and stop execution if not
    try:
        if not os.path.exists(name):
            raise FileNotFoundError(f"Directory {name} does not exist.")
    except FileNotFoundError as e:
        print(f"Error: {e}")
        raise SystemExit(1) from e
    return exist


# define a function which goes back to the original working directory if it is called
def chdir(dirname: str) -> None:
    """
    Change the active directory.
    """
    try:
        os.chdir(str(dirname))
        # print("Current working directory: {0}".format(os.getcwd()))
    except FileNotFoundError:
        print(f"Directory: {dirname} does not exist")
    except NotADirectoryError:
        print(f"{dirname} is not a directory")
    except PermissionError:
        print(f"You do not have permissions to change to {dirname}")

    return None


def checkifinpath(executable: str) -> None:
    try:
        sp.run(["which", executable], stdout=sp.DEVNULL, stderr=sp.DEVNULL, check=True)
    except sp.CalledProcessError as exc:
        raise FileNotFoundError(f"'{executable}' is not in PATH") from exc
    return None
