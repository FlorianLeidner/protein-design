import json
import argparse
from pathlib import Path
from typing import Optional, Union

import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer


def load_and_parameterize(
    pdb_path: Union[str, Path],
    ph: float = 7.0,
    add_missing_atoms: bool = True,
    add_missing_residues: bool = False,
    remove_heterogens: bool = False,
    constraints: app.HBonds = app.HBonds,
    force_field: str = "amber99sbildn.xml",
    rebuild_hydrogens: bool = True,
    normalize_charmm_termini: bool = True,
) -> mm.app.Simulation:
    """Load a PDB file, repair it with PDBFixer, and build an OpenMM Simulation
    parameterized with the AMBER99SB-ILDN force field ready for minimization.

    The system is set up in *vacuum* (no solvent, no periodic boundary conditions).
    All pairwise non-bonded interactions are evaluated (NoCutoff).

    Parameters
    ----------
    pdb_path:
        Path to the input PDB file.
    ph:
        pH used by PDBFixer when adding missing hydrogens.
    add_missing_atoms:
        Whether to ask PDBFixer to add missing heavy atoms / terminals.
    add_missing_residues:
        Whether to ask PDBFixer to model missing internal residue loops.
        Default is False because loop building can unexpectedly change a structure.
    remove_heterogens:
        Strip non-protein HETATM records (water, ligands …) before minimization.
    constraints:
        Bond-constraint scheme (default: constrain all bonds involving hydrogen).
    force_field:
        OpenMM force-field XML file name or path, e.g. ``amber99sbildn.xml``.
    rebuild_hydrogens:
        Delete all input hydrogens and add hydrogens with OpenMM Modeller using
        the requested force field. This is recommended when converting a PDB
        generated with a different force field, for example CHARMM -> AMBER.
    normalize_charmm_termini:
        Rename common CHARMM terminal oxygen names OT1/OT2 to AMBER/PDB-style
        O/OXT before template matching.

    Returns
    -------
    openmm.app.Simulation
        Fully initialized simulation object, positions loaded, ready to minimize.
    """
    pdb_path = Path(pdb_path)
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    print(f"Loading and fixing structure: {pdb_path}")

    # ---- PDBFixer: repair common PDB issues --------------------------------
    fixer = PDBFixer(filename=str(pdb_path))

    fixer.findMissingResidues()

    # PDBFixer stores missing residues in fixer.missingResidues.
    if not add_missing_residues:
        fixer.missingResidues = {}

    if remove_heterogens:
        fixer.removeHeterogens(keepWater=False)
    if add_missing_atoms:
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()

    topology = fixer.topology
    positions = fixer.positions

    if normalize_charmm_termini:
        _normalize_charmm_terminal_atom_names(topology)

    print(f"Structure has {topology.getNumAtoms()} atoms across {sum(1 for _ in topology.residues())} residues")

    # ---- Force field -------------------------------------------------------
    print(f"Using force field: {force_field}")
    forcefield = app.ForceField(force_field)

    # Use the target force field, not the source PDB's force-field convention,
    # to define protonation and terminal residues

    if rebuild_hydrogens:
        topology, positions = _rebuild_hydrogens(
            topology, positions, forcefield, ph
        )

    _report_unmatched_residues(forcefield, topology)

    # ---- System: vacuum, no PBC --------------------------------------------
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.NoCutoff,        # full electrostatics in vacuum
        constraints=constraints,
        removeCMMotion=True,
    )

    # ---- Integrator (only needed by Simulation; irrelevant for minimization) --
    integrator = mm.LangevinMiddleIntegrator(
        300 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )

    # ---- Platform: prefer CUDA  -------------------------------------------
    platform, platform_props = _select_platform()

    print(f"Using OpenMM platform: {platform.getName()}")

    simulation = app.Simulation(topology, system, integrator, platform, platform_props)
    simulation.context.setPositions(positions)

    return simulation

def minimize_energy(
    simulation: app.Simulation,
    tolerance: unit.Quantity = 10.0 * unit.kilojoule_per_mole / unit.nanometer,
    max_iterations: int = 0,
    output_path: Optional[Union[str, Path]] = None,
) -> tuple[mm.unit.Quantity, app.Topology]:
    """
    Run energy minimization on a parameterized simulation.

    Parameters
    ----------
    simulation:
        An OpenMM Simulation object with positions already set
        (e.g. the return value of :func:`load_and_parameterize`).
    tolerance:
        Convergence criterion: minimization stops when the RMS force on every
        atom falls below this value. Lower values give tighter convergence at
        the cost of more iterations.
    max_iterations:
        Maximum number of minimization steps. ``0`` means iterate until the
        tolerance is reached (recommended).
    output_path:
        If given, write the minimized structure to this PDB file.

    Returns
    -------
    positions : openmm.unit.Quantity
        Minimized atomic positions (Nx3 array in nanometers).
    topology : openmm.app.Topology
        Topology object (unchanged from input, provided for convenience).
    """
    state_before = simulation.context.getState(getEnergy=True)
    e_before = state_before.getPotentialEnergy().in_units_of(unit.kilocalorie_per_mole)
    print(f"Potential energy before minimization: {e_before._value:.3f} kcal/mol")

    iters_str = max_iterations if max_iterations > 0 else "unlimited"
    tol_val = tolerance.value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
    print(f"Starting minimization (tolerance={tol_val:.1f} kJ/mol/nm, max_iterations={iters_str})")

    simulation.minimizeEnergy(tolerance=tolerance, maxIterations=max_iterations)

    state_after = simulation.context.getState(getEnergy=True, getPositions=True)
    e_after = state_after.getPotentialEnergy().in_units_of(unit.kilocalorie_per_mole)
    print(f"Potential energy after  minimization: {e_after._value:.3f} kcal/mol")
    print(f"Energy change: {e_after._value - e_before._value:.3f} kcal/mol")

    positions = state_after.getPositions()

    if output_path is not None:
        _write_pdb(output_path, simulation.topology, positions)
        _write_report(output_path, e_before, e_after)

    return positions, simulation.topology


def minimize_pdb(
    pdb_path: Union[str, Path],
    output_path: Optional[Union[str, Path]] = None,
    force_field: str = "amber99sbildn.xml",
    ph: float = 7.0,
    remove_heterogens: bool = True,
    add_missing_residues: bool = False,
    rebuild_hydrogens: bool = True,
    normalize_charmm_termini: bool = True,
    tolerance: unit.Quantity = 10.0 * unit.kilojoule_per_mole / unit.nanometer,
    max_iterations: int = 0,
) -> tuple[mm.unit.Quantity, app.Topology]:
    """Load, parameterize, minimize, and optionally write one PDB structure."""
    simulation = load_and_parameterize(
        pdb_path,
        ph=ph,
        remove_heterogens=remove_heterogens,
        add_missing_residues=add_missing_residues,
        force_field=force_field,
        rebuild_hydrogens=rebuild_hydrogens,
        normalize_charmm_termini=normalize_charmm_termini,
    )
    return minimize_energy(
        simulation,
        tolerance=tolerance,
        max_iterations=max_iterations,
        output_path=output_path,
    )


def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser(
        description="Energy-minimize one or more PDB files using OpenMM."
    )
    parser.add_argument(
        "input_structures",
        nargs="+",
        type=Path,
        help="One or more input PDB files to minimize.",
    )
    parser.add_argument(
        "-f",
        "--force-field",
        default="amber99sbildn.xml",
        help="OpenMM force-field XML file name or path. Default: amber99sbildn.xml",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path("minimized_structures"),
        help="Directory where minimized PDB files will be written.",
    )
    parser.add_argument(
        "--ph",
        type=float,
        default=7.0,
        help="pH used when adding missing hydrogens. Default: 7.0",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=10.0,
        help="Minimization tolerance in kJ/mol/nm. Default: 10.0",
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=0,
        help="Maximum minimization iterations. Use 0 for unlimited. Default: 0",
    )
    parser.add_argument(
        "--add-missing-residues",
        action="store_true",
        help="Ask PDBFixer to model missing residue loops.",
    )

    parser.add_argument(
        "--remove-hetatms",
        dest="remove_hetatms",
        action="store_true",
        default=False,
        help="Remove HETATM records before minimization.",
    )
    parser.add_argument(
        "--keep-input-hydrogens",
        action="store_true",
        help=(
            "Do not rebuild hydrogens with the target force field. "
            "By default, input hydrogens are removed and regenerated for AMBER/target FF compatibility."
        ),
    )
    parser.add_argument(
        "--no-normalize-charmm-termini",
        action="store_true",
        help="Do not rename CHARMM terminal OT1/OT2 atoms to O/OXT.",
    )

    args = parser.parse_args()

    if not args.outdir.exists():
        print(f"\n=== Creating Output directory: {args.outdir} ===")
        args.outdir.parent.mkdir(parents=True, exist_ok=True)
    elif not args.outdir.isdir():
        parser.error(f"{args.outdir} is a file, expected directory")

    for filename in args.input_structures:
        if not Path(filename).exists():
            parser.error(f"{filename} not found")

    return args


def main() -> int:
    """CLI entry point: parameterize and minimize all requested input PDB files."""
    args = parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    tolerance = args.tolerance * unit.kilojoule_per_mole / unit.nanometer

    for pdb_path in args.input_structures:

        output_path = args.outdir / f"{pdb_path.stem}_minimized.pdb"

        print(f"\n=== Minimizing {pdb_path} -> {output_path} ===")
        minimize_pdb(
            pdb_path,
            output_path=output_path,
            force_field=args.force_field,
            ph=args.ph,
            remove_heterogens=args.remove_hetatms,
            add_missing_residues=args.add_missing_residues,
            rebuild_hydrogens=not args.keep_input_hydrogens,
            normalize_charmm_termini=not args.no_normalize_charmm_termini,
            tolerance=tolerance,
            max_iterations=args.max_iterations,
        )

    return 0


# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------


def _normalize_charmm_terminal_atom_names(topology: app.Topology) -> None:
    """
    Rename common CHARMM terminal oxygen names to PDB/AMBER names in-place.

    CHARMM C-termini often contain OT1/OT2.  AMBER/OpenMM protein templates
    generally expect the backbone carbonyl oxygen to be O and the terminal
    oxygen to be OXT.  Atom names are not the only matching criterion, but
    normalizing them makes PDBFixer/Modeller output much less ambiguous.
    """

    for chain in topology.chains():
        residues = list(chain.residues())
        if not residues:
            continue
        c_terminal = residues[-1]
        atoms_by_name = {atom.name: atom for atom in c_terminal.atoms()}
        if "OT1" in atoms_by_name and "O" not in atoms_by_name:
            atoms_by_name["OT1"].name = "O"
        if "OT2" in atoms_by_name and "OXT" not in atoms_by_name:
            atoms_by_name["OT2"].name = "OXT"


def _rebuild_hydrogens(
    topology: app.Topology,
    positions: unit.Quantity,
    forcefield: app.ForceField,
    ph: float,
) -> tuple[app.Topology, unit.Quantity]:
    """Strip hydrogens and add them back using the requested force field."""
    modeller = app.Modeller(topology, positions)
    hydrogens = [atom for atom in modeller.topology.atoms() if atom.element == app.element.hydrogen]
    if hydrogens:
        print(f"Removing {len(hydrogens)} input hydrogens before target-FF protonation")
        modeller.delete(hydrogens)

    print(f"Adding hydrogens with target force field at pH {ph}")
    modeller.addHydrogens(forcefield, pH=ph)
    return modeller.topology, modeller.positions


def _report_unmatched_residues(forcefield: app.ForceField, topology: app.Topology) -> None:
    """Print residues that still fail template matching, when supported by OpenMM."""
    try:
        unmatched = forcefield.getUnmatchedResidues(topology)
    except AttributeError:
        return
    if unmatched:
        print("Residues still unmatched by the selected force field:")
        for residue in unmatched:
            atom_names = ",".join(atom.name for atom in residue.atoms())
            print(f"  {residue.chain.id}:{residue.id} {residue.name} atoms=[{atom_names}]")


def _select_platform() -> tuple[mm.Platform, dict]:
    """Return the fastest available OpenMM platform and its properties."""
    preference = ["CUDA", "OpenCL", "CPU"]
    for name in preference:
        try:
            platform = mm.Platform.getPlatformByName(name)
            props = {"Precision": "mixed"} if name in ("CUDA", "OpenCL") else {}
            return platform, props
        except mm.OpenMMException:
            continue
    # Fallback should always be available
    return mm.Platform.getPlatformByName("Reference"), {}


def _write_pdb(path: Path, topology: app.Topology, positions) -> None:
    """Write a PDB file from an OpenMM topology and positions."""

    with open(path, "w") as fh:
        app.PDBFile.writeFile(topology, positions, fh)
    print(f"Minimized structure written to: {path}")


def _write_report(path: Path, e_before: float, e_after: float) -> None:
    """Log job data to json file"""

    report = {"e_before": f"{e_before:.3f}",
              "e_after": f"{e_after:.3f}",
              "delta_e": f"{e_before-e_after:.3f}",
              "unit": f"kcal/mol",
              "coordinates": f"{str(path)}"
            }

    path = path.with_suffix(".json")

    path.write_text(json.dumps(report))
    


if __name__ == "__main__":
    main()
