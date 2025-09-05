#!/bin/bash

# GROMACS MD Simulation Setup Script for Protein-Lysozyme Complexes
# Usage: ./setup_md.sh complex_number input_pdb_file

COMPLEX_NUM=$1
INPUT_PDB=$2

if [ $# -ne 2 ]; then
    echo "Usage: $0 <complex_number> <input_pdb_file>"
    exit 1
fi

# Set working directory
WORKDIR="/home/faruk/md_simulations/complex_${COMPLEX_NUM}"
cd $WORKDIR

echo "Setting up MD simulation for Complex $COMPLEX_NUM"
echo "Working directory: $WORKDIR"

# Load GROMACS module
module load GROMACS

# 1. Clean and prepare the PDB structure
echo "Step 1: Cleaning PDB structure..."
gmx pdb2gmx -ignh -f $INPUT_PDB -o input/processed.gro -p topology/topol.top -i topology/posre.itp << EOF
1
1
EOF

# 2. Define simulation box
echo "Step 2: Defining simulation box..."
gmx editconf -f input/processed.gro -o input/newbox.gro -c -d 1.0 -bt cubic

# 3. Solvate the system
echo "Step 3: Adding water molecules..."
gmx solvate -cp input/newbox.gro -cs spc216.gro -o input/solv.gro -p topology/topol.top

# 4. Add ions for neutralization
echo "Step 4: Adding ions..."
# Create ions.mdp file
cat > input/ions.mdp << 'MDP_EOF'
; ions.mdp - used as input into grompp to generate ions.tpr
integrator      = steep         ; Algorithm (steep = steepest descent minimization)
emtol           = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep          = 0.01          ; Energy minimization step size
nsteps          = 50000         ; Maximum number of (minimization) steps to perform
energygrps      = System        ; Which energy group(s) to write to disk
MDP_EOF

# Generate tpr file for adding ions
gmx grompp -f input/ions.mdp -c input/solv.gro -p topology/topol.top -o input/ions.tpr

# Add ions (choosing solvent group, typically group 13)
gmx genion -s input/ions.tpr -o input/solv_ions.gro -p topology/topol.top -pname NA -nname CL -neutral << EOF
13
EOF

echo ""
echo "========================================="
echo "System preparation completed for Complex $COMPLEX_NUM"
echo "========================================="
echo "Files created:"
echo "  - Topology: topology/topol.top"
echo "  - Final coordinates: input/solv_ions.gro"
echo "  - Position restraints: topology/posre.itp"
echo ""
echo "Next steps:"
echo "1. Run: ./create_mdp_files.sh $COMPLEX_NUM"
echo "2. Check the files and adjust parameters if needed"
echo "3. Submit the job: sbatch run_md_complex_${COMPLEX_NUM}.slurm"
echo ""
