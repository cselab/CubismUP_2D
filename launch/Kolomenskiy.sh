cd ../makefiles
make clean
make config=debug poisson=hypre bc=mixed precision=double particles=false -j
cd ../launch
#export OMP_NUM_THREADS=48;mpirun -np 32 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_Kolomenskiy_CFL0.1_bpd16 -CFL 0.1 -bpdx 16 -bpdy 16 -radius .0125 -tend 2. -rhoS 2 -ypos .85 -nu 0.0001551888 -sim falling -fdump 100
export OMP_NUM_THREADS=48;mpirun -np 32 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_Kolomenskiy_CFL0.1_bpd16 -CFL 0.1 -bpdx 16 -bpdy 16 -radius .05 -tend 2. -rhoS .001 -ypos .65 -nu 0.001 -sim falling -fdump 100