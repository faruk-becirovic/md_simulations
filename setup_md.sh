#!/bin/bash

# User-aware setup for GROMACS MD simulations
# This script detects the correct user and paths

echo "=== GROMACS MD Simulation Setup ==="
echo ""

# Check current user and setup paths accordingly
CURRENT_USER=$(whoami)
echo "Current user: $CURRENT_USER"

if [ "$CURRENT_USER" = "root" ]; then
    echo "Warning: You're running as root user"
    echo "It's recommended to use the 'faruk' user for this work"
    echo ""
    echo "Switch to faruk user with:"
    echo "su - faruk"
    echo ""
    echo "Or if you want to continue as root, we'll use /home/faruk as the working directory"
    WORK_BASE="/home/faruk"
elif [ "$CURRENT_USER" = "faruk" ]; then
    echo "✓ Running as faruk user - this is recommended"
    WORK_BASE="/home/faruk"
else
    echo "Running as user: $CURRENT_USER"
    WORK_BASE="/home/$CURRENT_USER"
fi

echo "Working base directory: $WORK_BASE"
echo ""

# Check if base directory exists and is writable
if [ ! -d "$WORK_BASE" ]; then
    echo "Error: Base directory $WORK_BASE does not exist!"
    exit 1
fi

if [ ! -w "$WORK_BASE" ]; then
    echo "Error: Cannot write to $WORK_BASE"
    echo "Check permissions or switch to the correct user"
    exit 1
fi

# Set up the simulation directory
SIMULATION_DIR="$WORK_BASE/md_simulations"

echo "Creating simulation directory structure..."
mkdir -p "$SIMULATION_DIR"
cd "$SIMULATION_DIR"

# Create directories for each complex
for i in {1..5}; do
    mkdir -p "complex_$i"/{input,topology,equilibration,production,analysis}
    echo "✓ Created complex_$i directories"
done

echo ""
echo "Directory structure created successfully!"
echo "Simulation directory: $SIMULATION_DIR"
echo ""

# Create a configuration file for other scripts to use
cat > config.sh << EOF
#!/bin/bash
# Configuration file for MD simulation scripts
WORK_BASE="$WORK_BASE"
SIMULATION_DIR="$SIMULATION_DIR"
CURRENT_USER="$CURRENT_USER"
EOF

echo "Configuration saved to: $SIMULATION_DIR/config.sh"
echo ""

# Check SLURM configuration based on user
echo "=== SLURM Configuration Check ==="

# Check if we can access SLURM
if command -v sinfo &> /dev/null; then
    echo "✓ SLURM is available"
    echo "Available partitions:"
    sinfo -s
    echo ""
    echo "Available nodes:"
    sinfo -N
else
    echo "✗ SLURM commands not available"
    echo "You may need to:"
    echo "1. Login to the login node (192.168.10.102)"
    echo "2. Check if SLURM is properly configured"
fi

echo ""
echo "=== Next Steps ==="
echo "1. Place your PDB files in: $SIMULATION_DIR/"
echo "   Name them: complex1.pdb, complex2.pdb, ..., complex5.pdb"
echo ""
echo "2. Update the setup scripts with the correct paths:"
echo "   - All scripts now use: $SIMULATION_DIR"
echo ""
echo "3. Run setup for each complex:"
echo "   cd $SIMULATION_DIR"
echo "   ./setup_md.sh 1 complex1.pdb"
echo ""
echo "4. If running as root, consider switching to faruk user:"
echo "   su - faruk"
echo "   cd $SIMULATION_DIR"

# Check NFS mount points mentioned in documentation
echo ""
echo "=== Storage Check ==="
if [ -d "/usr/local" ]; then
    echo "✓ /usr/local exists (mentioned in documentation for software)"
else
    echo "⚠ /usr/local not found - check NFS mounts"
fi

# Check if we're on the right node
echo ""
echo "=== Node Information ==="
echo "Hostname: $(hostname)"
echo "IP address: $(hostname -I | awk '{print $1}')"

case "$(hostname -I | awk '{print $1}')" in
    "192.168.10.102")
        echo "✓ You're on the login node - correct for setup"
        ;;
    "192.168.10.10")
        echo "ℹ You're on the main server"
        ;;
    "192.168.10.11"|"192.168.10.12"|"192.168.10.13")
        echo "ℹ You're on a compute node"
        ;;
    *)
        echo "ℹ Unknown node - check your connection"
        ;;
esac

echo ""
echo "Setup complete!"
