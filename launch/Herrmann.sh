cd ../makefiles/
make clean
make poisson=hypre bc=mixed config=production multiphase=true particles=true -j
cd -
export OMP_NUM_THREADS=48
#bsub -W 168:00 -n 48 -o Herrmann -J Herrmann mpirun -np 32 ../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/Herrmann -CFL 0.5 -bpdx 16 -bpdy 64 -tend .9 -sim rti -tdump 0.01
../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/Herrmann -CFL 0.5 -LCFL 0.1 -bpdx 2 -bpdy 8 -tend .9 -sim rti -tdump 0.01 -verbose