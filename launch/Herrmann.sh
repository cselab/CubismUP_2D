#cd ../makefiles/
#make clean
#make poisson=hypre bc=mixed config=production multiphase=true particles=false -j
#cd -
#export OMP_NUM_THREADS=48
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/HerrmannFDSplit_2 -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx 2 -bpdy 8 -tend .9 -sim rti -fdump 40 -nsteps 3600 -split
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/HerrmannFDSplit_4 -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx 4 -bpdy 16 -tend .9 -sim rti -fdump 80 -nsteps 7200 -split


cd ../makefiles/
make clean
make poisson=hypre bc=mixed config=production multiphase=true particles=false -j
#make poisson=hypre bc=mixed config=debug multiphase=true particles=false -j
cd -
export OMP_NUM_THREADS=48
../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/HerrmannFD_2 -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx 2 -bpdy 8 -tend .9 -sim rti -fdump 40 -nsteps 3600
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/HerrmannFD_4 -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx 4 -bpdy 16 -tend .9 -sim rti -fdump 80 -nsteps 7200
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/HerrmannFD_8 -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx 8 -bpdy 32 -tend .9 -sim rti -fdump 160 -nsteps 14400
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/HerrmannFD_16 -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx 16 -bpdy 64 -tend .9 -sim rti -fdump 320 -nsteps 28800

#cd ../makefiles/
#make clean
#make poisson=hypre bc=mixed config=production multiphase=true particles=true -j
#make poisson=hypre bc=mixed config=debug multiphase=true particles=true -j
#cd -
#export OMP_NUM_THREADS=48
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/HerrmannParticles_2 -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx 2 -bpdy 8 -tend .9 -sim rti -fdump 40 -nsteps 3600
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/HerrmannParticles_4 -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx 4 -bpdy 16 -tend .9 -sim rti -fdump 80 -nsteps 72800
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/HerrmannParticles_8 -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx 8 -bpdy 32 -tend .9 -sim rti -fdump 160 -nsteps 14400
#../makefiles/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/HerrmannParticles_16 -CFL 0.5 -LCFL 0.1 -nu 0.0023096 -rhoS 7.2314 -bpdx 16 -bpdy 64 -tend .9 -sim rti -fdump 320 -nsteps 28800