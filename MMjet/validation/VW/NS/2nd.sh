#!/bin/bash


########### SLURM CONFIGURATION ############

# Name of the job
#SBATCH --job-name=2nd_ns
#SBATCH --output=slurm-%j.out
#SBATCH --no-requeue

# Number of nodes, number of cores per node
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --ntasks-per-core=1

# Execution time
#SBATCH --time=460:00:00

# Quality of service
#SBATCH --qos=default-users

############################################


#----------------------------
# Config file name
CONFIG_FILE=INL.cfg


#----------------------------
# Mpi options
MPI_FLAGS='--bind-to core --report-bindings'

if test -d "ADAPT"; then
        rm -r ADAPT
fi

############################################

JOBID=$SLURM_JOB_ID

set -e
set -u

DIR=ADAPT

mkdir $DIR

# Creating config files: copying and removing leading blank lines
sed -n '/<DARWIN CONFIGURATION>/{:a;n;/<END OF DARWIN CONFIGURATION>/b;p;ba}' $CONFIG_FILE | sed '/./,$!d' > ${DIR}/darwin.cfg
sed -n '/<PREPRO CONFIGURATION>/{:a;n;/<END OF PREPR0 CONFIGURATION>/b;p;ba}' $CONFIG_FILE | sed '/./,$!d' > ${DIR}/prepro.cfg
sed -n '/<SU2 CONFIGURATION>/{:a;n;/<END OF SU2 CONFIGURATION>/b;p;ba}' $CONFIG_FILE | sed '/./,$!d' > ${DIR}/su2.cfg

# Get other configuratons
source <(sed -n '/<OTHER>/{:a;n;/<END OF OTHER>/b;p;ba}' $CONFIG_FILE)

# Copy cfg file for Darwin sparsekit
cp spkit.cfg ${DIR}

cd $DIR

# Get problem name, grid name and maximum adaptation level
PROB_NAME=$(sed '2q;d' darwin.cfg | awk '{print $1}')
GRID_NAME=$(sed '3q;d' darwin.cfg | awk '{print $1}')
MAX_ADAPT_LEVEL=$(sed '7q;d' darwin.cfg | awk '{print $1}')

# Copy initial grid
for f in ../${INIT_GRID_DIR}/*.${GRID_NAME}
do cp $f $(basename $f)_0; done
cp ../${INIT_GRID_DIR}/grid.${GRID_NAME}.su2 admesh.${GRID_NAME}.su2


# Set mesh file name and output file in su2.cfg
sed -i 's/^\s*MESH_FILENAME=.*$/MESH_FILENAME= 'admesh.${GRID_NAME}.su2'/' su2.cfg
sed -i 's/^\s*RESTART_FILENAME=.*$/RESTART_FILENAME= 'solution.dat'/' su2.cfg

# Initial solution
echo	
echo
echo '======================================================================='
echo '======================================================================='
echo '                       SU2 - INITIAL SOLUTION'
echo; echo; echo; echo "MPI - REPORT BINDINGS:"	

if [ "$INITIALIZE_FLOW_0" = YES ]; then
	cp ../$RESTART_FLOW_0 .
	sed -i 's/^\s*RESTART_SOL=.*$/RESTART_SOL= YES/' su2.cfg
	sed -i 's/^\s*SOLUTION_FILENAME=.*$/SOLUTION_FILENAME= '$RESTART_FLOW_0'/' su2.cfg
	sed -i 's/^\s*READ_BINARY_RESTART=.*$/READ_BINARY_RESTART= YES/' su2.cfg
else
	sed -i 's/^\s*RESTART_SOL=.*$/RESTART_SOL= NO/' su2.cfg
fi
if [ "$ADAPT_MESH_SIZE" = YES ]; then
	REF_LENGTH=$INIT_MESH_SIZE
	sed -i 's/^\s*REF_ELEM_LENGTH=.*$/REF_ELEM_LENGTH= '$REF_LENGTH'/' su2.cfg
fi

mpirun -n $SLURM_NTASKS $MPI_FLAGS SU2_CFD su2.cfg

# Set restart input file name in su2.cfg and restart

#sed -i 's/^\s*READ_BINARY_RESTART=.*$/READ_BINARY_RESTART= NO/' su2.cfg
sed -i 's/^\s*RESTART_SOL=.*$/RESTART_SOL= YES/' su2.cfg # THIS IS FOR RESTART AT HIGHER LEVELS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
sed -i 's/^\s*ITER=.*$/ITER= 10000/' su2.cfg
sed -i 's/^\s*MUSCL_FLOW=.*$/MUSCL_FLOW= YES/' su2.cfg
sed -i 's/^\s*SOLUTION_FILENAME=.*$/SOLUTION_FILENAME= 'solution_interpolated.dat'/' su2.cfg
#sed -i 's/^\s*FREESTREAM_PRESSURE=.*$/FREESTREAM_PRESSURE=1.125e6/' su2.cfg
#sed -i 's/^\s*MARKER_RIEMANN=.*$/MARKER_RIEMANN= ( 1, TOTAL_CONDITIONS_PT, 1.7e6, 525.9, 1.0, 0.0, 0.0, 4, TOTAL_CONDITIONS_PT, 3.2e5, 300, 1.0, 0.0, 0.0, 5, STATIC_PRESSURE, 3.2e5, 0.0, 0.0, 0.0, 0.0 )/' su2.cfg

cp solution.dat solution_interpolated.dat

# Get vtk of initial solution
sed -i 's/,/ /g' solution.csv
cp solution.csv solution.dat

if test -f "flow.vtu"; then
        mv flow.vtu ../flow_0.vtu
fi


# Adaptation loop

for ADAPT_LEVEL in $(seq 1 $MAX_ADAPT_LEVEL) 
do 

	sed -i '4s/.*/'${GRID_NAME}'_'$((ADAPT_LEVEL-1))'         ! Grid name/' prepro.cfg

	echo	
	echo
	echo '======================================================================='
	echo '======================================================================='
	echo '                       PREPROCESSING - LEVEL' $ADAPT_LEVEL

	prepro

	if [ ! -f node_pair.${GRID_NAME}_$((ADAPT_LEVEL-1)) ]; then
		exit
	fi 


	sed -i '8s/.*/'$((ADAPT_LEVEL-1))'	       ! Starting Level/' darwin.cfg

	echo	
	echo
	echo '======================================================================='
	echo '======================================================================='
	echo '                     MESH ADAPTATION - LEVEL' $ADAPT_LEVEL

	cp solution.dat solution.${PROB_NAME}
	darwin

	if [ ! -f grid.${GRID_NAME}_$ADAPT_LEVEL ]; then
		exit
	fi 

	echo	
	echo
	echo '======================================================================='
	echo '======================================================================='
	echo '                            SU2 - LEVEL' $ADAPT_LEVEL

	
	echo; echo; echo; echo "MPI - REPORT BINDINGS:"	

	if [ "$ADAPT_MESH_SIZE" = YES ]; then
		REF_LENGTH=$(awk "BEGIN {print ${REF_LENGTH}*0.7}")
		sed -i 's/^\s*REF_ELEM_LENGTH=.*$/REF_ELEM_LENGTH= '${REF_LENGTH}'/' su2.cfg
	fi

	mpirun -n $SLURM_NTASKS $MPI_FLAGS SU2_CFD su2.cfg

	cp solution.dat solution_interpolated.dat

	sed -i 's/,/ /g' solution.csv
	cp solution.csv solution.dat

	if test -f "flow.vtu"; then
		mv flow.vtu ../flow_$ADAPT_LEVEL.vtu
	fi
done

mv ../slurm-${SLURM_JOB_ID}.out .
mv history_2nd.csv ../

cd ../

# rm -rf ADAPT

############################################
############################################
# END OF Run.sh

