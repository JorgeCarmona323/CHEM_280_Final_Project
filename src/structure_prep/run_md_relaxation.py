"""
MD relaxation of stitched IDP+stable domain structures using OpenMM.

Restrains Cα atoms in the stable domain while allowing the IDP region
to relax freely.

Usage:
    python run_md_relaxation.py \
        --input data/structures/tau_stitched/stitched_tau_monomer.pdb \
        --output data/structures/tau_relaxed/tau_monomer_relaxed.pdb \
        --idp_residues 244-368 \
        --steps 50000
"""

import argparse
import os


def run_relaxation(
    input_pdb: str,
    output_pdb: str,
    idp_residues: tuple,
    steps: int = 50000,
    temperature: float = 300.0,
    forcefield: str = "amber14",
):
    try:
        from openmm import app, unit, LangevinMiddleIntegrator, CustomExternalForce
        from openmm.app import PDBFile, ForceField, Modeller, Simulation
        import openmm as mm
    except ImportError:
        raise ImportError(
            "OpenMM not installed. Install with: conda install -c conda-forge openmm"
        )

    idp_start, idp_end = idp_residues
    print(f"Loading structure: {input_pdb}")
    pdb = PDBFile(input_pdb)

    # Force field
    if forcefield == "amber14":
        ff = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")
    else:
        ff = ForceField("charmm36.xml", "charmm36/water.xml")

    # Add hydrogens and solvate
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.addHydrogens(ff)
    modeller.addSolvent(ff, model="tip3p", padding=1.0 * unit.nanometers)

    print("Building system...")
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometers,
        constraints=app.HBonds,
    )

    # Restrain Cα of stable domain residues (outside IDP range)
    restraint = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", 1000.0 * unit.kilojoules_per_mole / unit.nanometers**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    for atom in modeller.topology.atoms():
        resnum = atom.residue.index + 1  # 1-indexed
        if atom.name == "CA" and not (idp_start <= resnum <= idp_end):
            pos = modeller.positions[atom.index]
            restraint.addParticle(atom.index, [pos.x, pos.y, pos.z])

    system.addForce(restraint)

    # Integrator
    integrator = LangevinMiddleIntegrator(
        temperature * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )

    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)

    # Energy minimize
    print("Energy minimizing...")
    simulation.minimizeEnergy(maxIterations=1000)

    # NVT equilibration
    print(f"Running {steps} MD steps (~{steps * 0.002 / 1000:.0f} ps)...")
    simulation.step(steps)

    # Write output
    os.makedirs(os.path.dirname(output_pdb) or ".", exist_ok=True)
    positions = simulation.context.getState(getPositions=True).getPositions()

    with open(output_pdb, "w") as f:
        PDBFile.writeFile(simulation.topology, positions, f)

    print(f"Relaxed structure written to {output_pdb}")


def main():
    parser = argparse.ArgumentParser(description="OpenMM MD relaxation for stitched IDP structures")
    parser.add_argument("--input", required=True, help="Input PDB (stitched structure)")
    parser.add_argument("--output", required=True, help="Output relaxed PDB")
    parser.add_argument("--idp_residues", default="244-368",
                        help="IDP residue range (unrestrained), e.g. 244-368")
    parser.add_argument("--steps", type=int, default=50000, help="MD steps (default: 50000 = 100ps)")
    parser.add_argument("--temperature", type=float, default=300.0, help="Temperature in Kelvin")
    parser.add_argument("--forcefield", default="amber14", choices=["amber14", "charmm36"])
    args = parser.parse_args()

    idp_start, idp_end = map(int, args.idp_residues.split("-"))
    run_relaxation(
        input_pdb=args.input,
        output_pdb=args.output,
        idp_residues=(idp_start, idp_end),
        steps=args.steps,
        temperature=args.temperature,
        forcefield=args.forcefield,
    )


if __name__ == "__main__":
    main()
