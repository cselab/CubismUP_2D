

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=periodic precision=double bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 2 -maxBPD 32 -test advection -ic 0 > results_advection_ic0
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 8 -test advection -ic 1 > results_advection_ic1
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 2 -maxBPD 128 -test advection -ic 1 > results_advection_ic1
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 64 -test poisson -solver 0 -ic 0 > results_poisson_solver0_ic0
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 8 -test poisson -solver 0 -ic 1 > results_poisson_solver0_ic1
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 64 -test poisson -solver 1 -ic 0 > results_poisson_solver1_ic0
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 64 -test poisson -solver 1 -ic 1 > results_poisson_solver1_ic1
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test poisson -solver 2 -ic 1 > results_poisson_solver2_ic1










#=============================
# Passed - primary importance
#=============================
#passed - with double precision
cd ../makefiles
make clean
make config=production poisson=hypre bc=periodic precision=double bs=32 -j
cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 1 -maxBPD 64 -test diffusion > results_diffusion
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 32 -test gravity > results_gravity
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 32 -test penalization > results_penalization
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 32 -test translation -ic 0 > results_translation_ic0
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 32 -test translation -ic 1 > results_translation_ic1
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 32 -test rotation -ic 0 > results_rotation_ic0
export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 32 -test rotation -ic 1 > results_rotation_ic1

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=periodic precision=double bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 3 -ic 0 > results_pressure_solver3_ic0
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 3 -ic 1 > results_pressure_solver3_ic1

#passed? - can't do convergence in physical space for spectral method
#cd ../makefiles
#make clean
#make config=production bc=periodic precision=single bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 64 -test pressure -solver 0 -ic 0 > results_pressure_solver0_ic0
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 64 -test pressure -solver 0 -ic 1 > results_pressure_solver0_ic1
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 64 -test pressure -solver 1 -ic 0 > results_pressure_solver1_ic0
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 64 -test pressure -solver 1 -ic 1 > results_pressure_solver1_ic1

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=mixed precision=double bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 3 -ic 2 > results_pressure_solver3_ic2
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test poisson -solver 2 -ic 0 > results_poisson_solver2_ic0


#==============================
# Failed? secondary importance
#==============================

#failed? - fortran
#cd ../makefiles
#make clean
#make config=production bc=periodic precision=single bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 64 -test pressure -solver 2 -ic 0 > results_pressure_solver2_ic0
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 64 -test pressure -solver 2 -ic 1 > results_pressure_solver2_ic1

#cd ../makefiles
#make clean
#make config=production poisson=fortran bc=mixed precision=single bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 4 -maxBPD 64 -test pressure -solver 2 -ic 2 > results_pressure_solver2_ic2

#failed
#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=periodic precision=double bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 4 -ic 0 > results_pressure_solver4_ic0
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 4 -ic 1 > results_pressure_solver4_ic1

#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=mixed precision=double bs=32 -j
#cd ../launch/
#export OMP_NUM_THREADS=48;../makefiles/test -minBPD 8 -maxBPD 64 -test pressure -solver 4 -ic 2 > results_pressure_solver4_ic2