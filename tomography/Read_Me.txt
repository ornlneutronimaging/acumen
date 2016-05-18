%% Work Flow %%


1) Convert fits file on mac to dat.  First env var are set, code is compiled and exectuted.  It is assumed that all projection data is in ROOT_CT_DIR/projections

ROOT_CT_DIR="/Users/2ya/Documents/Working/Work/python/tomography/fortran/main_codes/data"
EXP_CT_NAME="CT_Aluminum_ironcenter_0090"
export ROOT_CT_DIR EXP_CT_NAME
cd $ROOT_CT_DIR
cd $EXP_CT_NAME
rm -rf reconstruction
cd ../..
mpif90 fits_to_hpc.F90 -o fits_to_hpc -lcfitsio
mpirun -np 1 ./fits_to_hpc 
EXP_CT_NAME="Derek_injec_0040"
cd $ROOT_CT_DIR
cd $EXP_CT_NAME
rm -rf reconstruction
cd ../..
mpirun -np 1 ./fits_to_hpc 
EXP_CT_NAME="TURBINECT_0010"
cd $ROOT_CT_DIR
cd $EXP_CT_NAME
rm -rf reconstruction
cd ../..
mpirun -np 1 ./fits_to_hpc 
EXP_CT_NAME="TURBINECT_0180"
cd $ROOT_CT_DIR
cd $EXP_CT_NAME
rm -rf reconstruction
cd ../..
mpirun -np 1 ./fits_to_hpc 

sftp 2ya@dtn.ccs.ornl.gov
cd /lustre/atlas1/cli106/scratch/2ya/tomo/tomo_real/data
cd CT_Aluminum_ironcenter_0090/reconstruction


2) Convert projections to sinograms.  Following is work flow of logging on to TITAN and running code.  Assumes projections ftp to TITAN

qsub -I -lnodes=32 -lwalltime=02:00:00 -X -A mat134
cd /lustre/atlas1/cli106/scratch/2ya/tomo/tomo_real

ROOT_CT_DIR="/lustre/atlas1/cli106/scratch/2ya/tomo/tomo_real/data"
EXP_CT_NAME="CT_Aluminum_ironcenter_0090"
export ROOT_CT_DIR EXP_CT_NAME
module load PrgEnv-pgi
module load cudatoolkit
export CRAY_CUDA_PROXY=1

ftn -o projections_to_sinogram_hpc projections_to_sinogram_hpc.F90 -Minfo=all -ta=tesla,cc35,cuda7.0,ptxinfo -acc -Minfo=accel 
aprun -n 2 -N 1 ./projections_to_sinogram_hpc

3) Run Filtered Back Projection

ftn -o fbp_hpc fbp_hpc.F90 -Minfo=all -ta=tesla,cc35,cuda7.0,ptxinfo -acc -Minfo=accel 
aprun -n 32 -N 1 ./fbp_hpc


4) Run TV  (note, can run at sub resolution for r and theta by changing nsub and nsub_t in input_parm.txt)

ftn -o tv_hpc tv_hpc.F90 -Minfo=all -ta=tesla,cc35,cuda7.0,ptxinfo -acc -Minfo=accel 
aprun -n 8 -N 1 ./tv_hpc

Note:  If you just want to run a slice;

ftn -o tv_hpc_slice tv_hpc_slice.F90 -Minfo=all -ta=tesla,cc35,cuda7.0,ptxinfo -acc -Minfo=accel 
aprun -n 1 -N 1 ./tv_hpc_slice