cd ../makefiles
make clean
make config=production poisson=hypre bc=periodic precision=double bs=32 -j
cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 2 -maxBPD 32 -test advection -ic 0 > results_advection_ic0
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 2 -maxBPD 128 -test advection -ic 1 > results_advection_ic1
export OMP_NUM_THREADS=48;mpirun -np 1 ../makefiles/test -minBPD 8 -maxBPD 64 -test poisson -ic 1 > results_poisson_solver2_MG_ic1_varCoeffs



#=============================
# Passed - primary importance
#=============================
#passed - with double precision
./testDiffusion.sh
cd ../makefiles
make clean
make config=production poisson=hypre bc=periodic precision=double bs=32 -j
cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test diffusion > results_diffusion
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test gravity > results_gravity
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test penalization -minDT 1e-4 -maxDT 1e-4 > results_penalization
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test translation -minDT 1e-2 -maxDT 1e-2 -ic 0 > results_translation_ic0_body
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test translation -minDT 1e-2 -maxDT 1e-2 -ic 1 > results_translation_ic1_fromFlow
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test rotation -minDT 1 -maxDT 1 -ic 0 > results_rotation_ic0_body
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test rotation -minDT 1 -maxDT 1 -ic 1 > results_rotation_ic1_fromFlow

cd ../makefiles
make clean
make config=production poisson=hypre bc=periodic precision=double bs=32 -j
cd ../launch/
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -minDT 1 -maxDT 1 -solver 2 -ic 0 > results_pressure_solver2_MG_ic0_periodic
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -minDT 1 -maxDT 1 -solver 2 -ic 1 > results_pressure_solver2_MG_ic1_vortex

#passed? - can't do convergence in physical space for spectral method
cd ../makefiles
make clean
make config=production bc=periodic precision=single bs=32 -j
cd ../launch/
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test pressure -minDT 1 -maxDT 1 -solver 0 -ic 0 > results_pressure_solver0_FFTW_ic0
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test pressure -minDT 1 -maxDT 1 -solver 0 -ic 1 > results_pressure_solver0_FFTW_ic1
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test pressure -minDT 1 -maxDT 1 -solver 1 -ic 0 > results_pressure_solver1_FFTWspectral_ic0
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test pressure -minDT 1 -maxDT 1 -solver 1 -ic 1 > results_pressure_solver1_FFTWspectral_ic1

cd ../makefiles
make clean
make config=production poisson=hypre bc=mixed precision=double bs=32 -j
cd ../launch/
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -minDT 1 -maxDT 1 -solver 2 -ic 2 > results_pressure_solver2_MG_ic2_mixedBC
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test poisson -ic 0 > results_poisson_ic0_constCoeffs


#==============================
# Failed? secondary importance
#==============================

#failed
#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=periodic precision=double bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 3 -ic 0 > results_pressure_solver3_ic0
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 3 -ic 1 > results_pressure_solver3_ic1

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=mixed precision=double bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 3 -ic 2 > results_pressure_solver3_ic2