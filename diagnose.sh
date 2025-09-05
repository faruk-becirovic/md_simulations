#!/bin/bash

# Diagnostic script for GROMACS MD setup
# Usage: ./diagnose.sh [complex_number]

COMPLEX_NUM=$1

echo "=== GROMACS MD Setup Diagnostics ==="
echo ""

# Check GROMACS installation
echo "1. Checking GROMACS installation:"
if command -v gmx &> /dev/null; then
    gmx --version | head -5
    echo "✓ GROMACS is available"
else
    echo "✗ GROMACS not found in PATH"
    echo "Try: module load GROMACS"
    echo "Available modules:"
    module avail 2>&1 | grep -i gromacs || echo "No GROMACS modules found"
fi

echo ""

# Check current directory and files
echo "2. Current directory and files:"
echo "Current directory: $(pwd)"
echo "Contents:"
ls -la

echo ""

# If complex number provided, check that complex specifically
if [ -n "$COMPLEX_NUM" ]; then
    WORKDIR="/usr/home/faruk/md_simulations/complex_${COMPLEX_NUM}"
    echo "3. Checking Complex $COMPLEX_NUM setup:"
    echo "Working directory: $WORKDIR"

    if [ -d "$WORKDIR" ]; then
        cd "$WORKDIR"
        echo "Directory structure:"
        find . -type f | sort

        echo ""
        echo "Checking critical files:"

        # Check topology
        if [ -f "topology/topol.top" ]; then
            echo "✓ topology/topol.top exists"
            echo "  Number of lines: $(wc -l < topology/topol.top)"
        else
            echo "✗ topology/topol.top missing"
        fi

        # Check coordinates
        if [ -f "input/solv_ions.gro" ]; then
            echo "✓ input/solv_ions.gro exists"
            echo "  Number of atoms: $(head -2 input/solv_ions.gro | tail -1 | awk '{print $1}')"
        else
            echo "✗ input/solv_ions.gro missing"

            # Check intermediate files
            if [ -f "input/solv.gro" ]; then
                echo "  → input/solv.gro exists (solvation completed)"
            else
                echo "  → input/solv.gro missing (solvation failed)"
            fi

            if [ -f "input/newbox.gro" ]; then
                echo "  → input/newbox.gro exists (box definition completed)"
            else
                echo "  → input/newbox.gro missing (box definition failed)"
            fi

            if [ -f "input/processed.gro" ]; then
                echo "  → input/processed.gro exists (pdb2gmx completed)"
            else
                echo "  → input/processed.gro missing (pdb2gmx failed)"
            fi
        fi

        # Check MDP files
        echo ""
        echo "MDP files:"
        for mdp in input/*.mdp; do
            if [ -f "$mdp" ]; then
                echo "✓ $mdp exists"
            fi
        done

    else
        echo "✗ Complex $COMPLEX_NUM directory doesn't exist"
        echo "Expected: $WORKDIR"
    fi
else
    echo "3. Checking overall setup:"
    if [ -d "/usr/home/faruk/md_simulations" ]; then
        echo "✓ Main simulation directory exists"
        echo "Complexes found:"
        ls -d /usr/home/faruk/md_simulations/complex_* 2>/dev/null || echo "No complexes found"
    else
        echo "✗ Main simulation directory missing"
        echo "Expected: /usr/home/faruk/md_simulations"
    fi
fi

echo ""

# Check common issues
echo "4. Common troubleshooting:"
echo ""
echo "If pdb2gmx fails:"
echo "  - Check PDB format: grep '^ATOM' yourfile.pdb | head -5"
echo "  - Remove HETATM except ligands: grep -v '^HETATM' yourfile.pdb > clean.pdb"
echo "  - Check for missing heavy atoms"
echo ""
echo "If solvation fails:"
echo "  - Check if processed.gro has correct format"
echo "  - Verify box was created properly"
echo ""
echo "If ion addition fails:"
echo "  - List available groups: echo 'q' | gmx genion -s input/ions.tpr"
echo "  - Choose SOL or water group (usually 13, sometimes different)"
echo ""

# System information
echo "5. System information:"
echo "Hostname: $(hostname)"
echo "User: $(whoami)"
echo "Available space: $(df -h . | tail -1 | awk '{print $4}')"
echo "Memory: $(free -h | grep Mem | awk '{print $7}') available"

echo ""
echo "=== End Diagnostics ==="
