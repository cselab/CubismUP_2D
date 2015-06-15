cd ../makefiles/
make clean
make poisson=hypre bc=mixed config=production multiphase=true particles=false -j
cd -
export OMP_NUM_THREADS=48
../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/BubbleUp -CFL 0.5 -LCFL 0.1 -bpdx 8 -bpdy 8 -sim bubble -rhoS 0.9 -radius 0.025 -tdump .1 -CFL 0.1