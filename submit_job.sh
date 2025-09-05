#!/bin/bash

# Create SLURM job scripts for GROMACS MD simulations
# Usage: ./create_slurm_jobs.sh

echo "Creating SLURM job scripts for MD simulations"

# Master submission script
cat > submit_all_complexes.sh << 'MASTER_EOF'
#!/bin/bash

# Submit all 5 complexes to SLURM queue
# Each complex will be assigned to a different GPU node

for i in {1..5}; do
    echo "Submitting Complex $i"
    sbatch run_md_complex_${i}.slurm
    sleep 2  # Small delay between submissions
done

echo "All jobs submitted. Check status with: squeue"
MASTER_EOF

chmod +x submit_all_complexes.sh

# Create individual SLURM scripts for each complex
for COMPLEX_NUM in {1..5}; do

cat > run_md_complex_${COMPLEX_NUM}.slurm << SLURM_EOF
#!/bin/bash
#SBATCH --job-name=md_complex_${COMPLEX_NUM}
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28        # Use all CPU cores on one node
#SBATCH --gres=gpu:4              # Use all 4 GPUs per node
#SBATCH --time=24:00:00           # Adjust based on your simulation length
#SBATCH --output=complex_${COMPLEX_NUM}/complex_${COMPLEX_NUM}_%j.out
#SBATCH --error=complex_${COMPLEX_NUM}/complex_${COMPLEX_NUM}_%j.err

# Load required modules
module load GROMACS

# Set working directory
WORKDIR="/usr/home/md_simulations/complex_${COMPLEX_NUM}"
cd \$WORKDIR

echo "Starting MD simulation for Complex ${COMPLEX_NUM}"
echo "Job ID: \$SLURM_JOB_ID"
echo "Node: \$SLURMD_NODENAME"
echo "Working directory: \$PWD"
echo "Start time: \$(date)"

# Set OpenMP threads (should equal cpus-per-task)
export OMP_NUM_THREADS=28

# Step 1: Energy minimization
echo "Step 1: Energy minimization"
gmx grompp -f input/minim.mdp -c input/solv_ions.gro -p topology/topol.top -o equilibration/em.tpr
gmx mdrun -v -deffnm equilibration/em -gpu_id 0123

# Step 2: NVT equilibration
echo "Step 2: NVT equilibration"
gmx grompp -f input/nvt.mdp -c equilibration/em.gro -r equilibration/em.gro -p topology/topol.top -o equilibration/nvt.tpr
gmx mdrun -v -deffnm equilibration/nvt -gpu_id 0123

# Step 3: NPT equilibration
echo "Step 3: NPT equilibration"
gmx grompp -f input/npt.mdp -c equilibration/nvt.gro -r equilibration/nvt.gro -t equilibration/nvt.cpt -p topology/topol.top -o equilibration/npt.tpr
gmx mdrun -v -deffnm equilibration/npt -gpu_id 0123

# Step 4: Production MD
echo "Step 4: Production MD simulation"
gmx grompp -f input/md.mdp -c equilibration/npt.gro -t equilibration/npt.cpt -p topology/topol.top -o production/md_0_1.tpr
gmx mdrun -v -deffnm production/md_0_1 -gpu_id 0123

echo "Simulation completed for Complex ${COMPLEX_NUM}"
echo "End time: \$(date)"

# Basic analysis
echo "Running basic analysis..."
cd analysis

# RMSD analysis
echo "1 4 \n 4 4" | gmx rms -s ../production/md_0_1.tpr -f ../production/md_0_1.xtc -o rmsd.xvg -tu ns

# Radius of gyration
echo "1" | gmx gyrate -s ../production/md_0_1.tpr -f ../production/md_0_1.xtc -o gyrate.xvg

# Potential energy
gmx energy -f ../production/md_0_1.edr -o potential.xvg << ENERGY_EOF
10
0
ENERGY_EOF

echo "Basic analysis completed for Complex ${COMPLEX_NUM}"
echo "Check the analysis/ directory for results"

SLURM_EOF

done

echo "SLURM job scripts created:"
for i in {1..5}; do
    echo "  - run_md_complex_${i}.slurm"
done
echo "  - submit_all_complexes.sh (master submission script)"
