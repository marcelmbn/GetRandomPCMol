"""
Command line interface for the getrandompcmol package.
"""
from __future__ import annotations

import argparse as ap

from .evaluate_calc import create_res_dir, eval_calc_ensemble, get_calc_ensemble
from .evaluate_conf import eval_conf_ensemble
from .main import main
from .miscelleanous import bcolors, checkifinpath


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
        help="Conformer sampling with CREST. \
Provide a keyword for a run mode: 'normal', 'protonate'. \
The default is 'normal'.",
        required=False,
        type=str,
        choices=["normal", "protonate"],
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
        type=str,
        choices=["wipe", "keep"],
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

    # raise an error if --evalcalconly is used with any other argument
    if args.evalcalconly and (args.opt or args.crest or args.evalconf):
        print(
            f"{bcolors.FAIL}You cannot use --evalcalconly with any other argument.{bcolors.ENDC}"
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
