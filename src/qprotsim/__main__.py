import argparse
import warnings
import pathlib
import sys

warnings.filterwarnings("ignore")

def main():
    parser = argparse.ArgumentParser(
        description=(
            "qprotsim: A package for studying thermodynamic properties "
            "of proteins using Quantum Computing methods."
        ),
        epilog="Use '<command> --help' for more information on a specific command.",
    )

    subparsers = parser.add_subparsers(
        dest="command",
        title="Available Commands",
        metavar=""
    )

    # -------------------------
    # PREPROCESS SUBCOMMAND
    # -------------------------
    preprocess_parser = subparsers.add_parser(
        "preprocess",
        help="Preprocess a protein structure to extract molecular graph features.",
        description=(
            "Preprocess a single protein file (PDB) and contruct the molecular graph useful for quantum state initialisation"
        )
    )

    preprocess_parser.add_argument(
        "-f", "--pdb",
        type=pathlib.Path,
        required=True,
        help="Path to the input PDB file."
    )

    preprocess_parser.add_argument(
        "-s", "--selection",
        type=str,
        default="protein",
        help="Atom selection (default: 'protein')."
    )

    preprocess_parser.add_argument(
        "-c", "--chain",
        type=str,
        default="A",
        help="Chain ID to extract (default: A)."
    )

    preprocess_parser.add_argument(
        "-o", "--output",
        type=pathlib.Path,
        required=True,
        help="Path to save the preprocessed output file."
    )

    # -------------------------
    # ANALYZE SUBCOMMAND
    # -------------------------
    analyze_parser = subparsers.add_parser(
        "analyze",
        aliases=["analyse", "analysis", "analyze-protein"],
        help="Run the analysis pipeline and generate a report from YAML config.",
        description="Execute the full pipeline, generate plots, and create a reportlab PDF.",
    )
    analyze_parser.add_argument(
        "-c", "--config",
        type=pathlib.Path,
        required=True,
        help="Path to YAML config file.",
    )
    analyze_parser.add_argument(
        "-o", "--output-dir",
        type=pathlib.Path,
        default=None,
        help="Optional output directory for plots and CSV.",
    )
    analyze_parser.add_argument(
        "-r", "--report",
        type=pathlib.Path,
        default=None,
        help="Optional report PDF output path.",
    )
    analyze_parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Enable verbose debug output (default: False).",
    )

    # -------------------------
    # BASIN SUBCOMMAND
    # -------------------------
    basin_parser = subparsers.add_parser(
        "basin",
        help="Run basin analysis.",
        description="Run basin analysis from PDB, samples, and energy parameters.",
    )
    # basin_parser.add_argument(
    #     "-i", "--pdb",
    #     type=pathlib.Path,
    #     required=True,
    #     help="Path to the input PDB file.",
    # )
    basin_parser.add_argument(
        "-z", "--z_samples",
        type=pathlib.Path,
        required=True,
        help="Path to Z-samples input.",
    )
    basin_parser.add_argument(
        "--E0",
        type=pathlib.Path,
        default=None,
        help="Optional path to E0 .npy file. If omitted, values are derived from circuit.energies_ref().",
    )
    basin_parser.add_argument(
        "--E1",
        type=pathlib.Path,
        default=None,
        help="Optional path to E1 .npy file. If omitted, values are derived from circuit.energies_ref().",
    )
    basin_parser.add_argument(
        "-b", "--basins",
        type=float,
        nargs="+",
        required=True,
        help="Basin boundaries as floats (must contain an even number of values).",
    )
    basin_parser.add_argument(
        "-o", "--out",
        type=pathlib.Path,
        required=True,
        help="Path to output.",
    )
    basin_parser.add_argument(
        "--flip_m",
        type=int,
        default=2000,
        help="Flip count parameter (default: 2000).",
    )

    # -------------------------
    # Parse arguments
    # -------------------------
    args = parser.parse_args()

    # If no command is given → show top-level help
    if args.command is None:
        parser.print_help()
        return

    # -------------------------
    # COMMAND ROUTING
    # -------------------------
    if args.command == "preprocess":
        run_preprocess(args)
    elif args.command in ("analyze", "analyse", "analysis", "analyze-protein"):
        run_analyze(args)
    elif args.command == "basin":
        if len(args.basins) % 2 != 0:
            parser.error("--basins must contain an even number of floats.")
        run_basin(args)

def run_preprocess(args):
    print("\nRunning preprocess with:")
    print(f" PDB file:   {args.pdb}")
    print(f" Selection:  {args.selection}")
    print(f" Chain:      {args.chain}")
    print(f" Output:     {args.output}\n")


def run_analyze(args):
    from qprotsim.report_analysis import run_analysis_from_yaml, run_batch_analysis_from_yaml, load_yaml_config
    import qprotsim.info_theory as info_theory

    try:
        # Set debug mode if requested
        info_theory.DEBUG = bool(args.debug)
        if args.debug:
            print("[DEBUG MODE ENABLED]\n")

        # Check if batch mode
        config = load_yaml_config(args.config)
        is_batch = config.get("batch_mode", False)

        if is_batch:
            print("\n[INFO] Batch mode detected - processing multiple proteins\n")
            reports = run_batch_analysis_from_yaml(
                config_path=args.config,
                output_dir=args.output_dir,
            )
            print(f"\n[COMPLETE] Generated {len(reports)} reports:")
            for name, path in reports.items():
                print(f"  - {name}: {path}")
        else:
            output_dir = args.output_dir
            report_filename = None
            if args.report is not None:
                report_filename = args.report.name
                if output_dir is None and args.report.parent != pathlib.Path("."):
                    output_dir = args.report.parent
            report_path = run_analysis_from_yaml(
                config_path=args.config,
                output_dir=output_dir,
                report_filename=report_filename,
            )
            print(f"[saved] report -> {report_path}")

    except FileNotFoundError as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"\n{'='*70}", file=sys.stderr)
        print(f"ERROR: Configuration Error", file=sys.stderr)
        print(f"{'='*70}", file=sys.stderr)
        print(f"{e}", file=sys.stderr)
        print(f"{'='*70}\n", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"\n{'='*70}", file=sys.stderr)
        print(f"ERROR: Unexpected error occurred", file=sys.stderr)
        print(f"{'='*70}", file=sys.stderr)
        print(f"{type(e).__name__}: {e}", file=sys.stderr)
        print(f"{'='*70}\n", file=sys.stderr)
        if args.debug:
            print("\nFull traceback:", file=sys.stderr)
            raise
        sys.exit(1)


def run_basin(args):
    from qprotsim.info_theory import run_basin as run_basin_analysis

    return run_basin_analysis(args)

if __name__ == "__main__":
    main()
