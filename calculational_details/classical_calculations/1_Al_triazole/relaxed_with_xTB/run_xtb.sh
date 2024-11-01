#!/bin/bash

# Exit on error
set -e

# Set up environment variables
echo "Setting up environment variables..."
export PATH=$HOME/xtb/bin:$PATH
export XTBPATH=$HOME/xtb/share/xtb
export LD_LIBRARY_PATH=$HOME/xtb/lib:$LD_LIBRARY_PATH

# Configure stack size and OpenMP settings
ulimit -s unlimited
export OMP_STACKSIZE=4G
export OMP_NUM_THREADS=48,1
export OMP_MAX_ACTIVE_LEVELS=1
export MKL_NUM_THREADS=48

# Create calculation directory
CALC_DIR="xtb_calc_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$CALC_DIR"
cd "$CALC_DIR"

# Create .CHRG file (neutral system)
echo "0" > .CHRG

# Create xcontrol file with enhanced settings for periodic metallic system
cat > .xcontrol << EOF
$set
  maxiter=1000
  broydamp=0.2
  mixing=adiis
  econv=1e-4
  acc=1.0
  scfconv=1e-4
  electronic_temperature=10000.0
  max_scc_iter=1000
$end

$periodic
  cell periodic=2
$end

$scf
  mixer=damp
  mixing=0.1
  maxiterations=1000
$end

$opt
   maxcycle=500
   optlevel=crude
   microcycle=100
   maxdispl=0.1
   hlow=0.01
$end
EOF

# Copy input file
cp ../"$1" ./input.xyz

# Run optimization in background
{
    echo "XTB calculation started at $(date)"
    echo "Input file: $1"
    echo "Working directory: $CALC_DIR"
    echo "----------------------------------------"

    # First run with very high temperature for better initial guess
    echo "Running initial high-temperature calculation..."
    xtb input.xyz --gfn 1 --etemp 10000.0 --input .xcontrol

    if [ -f "xtbrestart" ]; then
        echo "----------------------------------------"
        echo "Initial run completed, starting optimization..."
        # Then run optimization with still high temp
        xtb input.xyz --opt crude --gfn 1 --etemp 10000.0 --input .xcontrol --cycles 1000
    fi

    echo "----------------------------------------"
    if [ -f "xtbopt.xyz" ]; then
        echo "Optimization completed successfully at $(date)"
        echo "Optimized geometry saved in xtbopt.xyz"
        cp xtbopt.xyz ../optimized_structure.xyz
    else
        echo "Warning: Optimization may have failed. Check output files."
        if [ -f "xtblast.xyz" ]; then
            echo "Last geometry available in xtblast.xyz"
            cp xtblast.xyz ../last_structure.xyz
        fi
    fi
} > "xtb_${1%.*}.log" 2>&1 &

# Get the background process ID
bg_pid=$!

# Create a PID file
echo $bg_pid > xtb.pid

echo "XTB calculation started in background with PID: $bg_pid"
echo "Log file: xtb_${1%.*}.log"
echo "Working directory: $CALC_DIR"
echo "To monitor the progress, use: tail -f xtb_${1%.*}.log"
echo "To check if process is still running: ps -p $bg_pid"