cd ../makefiles/;make clean;make bc=periodic config=production multiphase=true precision=single vertexcentered=true poisson=split-fftw particles=false densitydiff=true rk2=false -j;cd -
export OMP_NUM_THREADS=48;bsub -n 48 -W 08:00 ../makefiles/test -minBPD 8 -maxBPD 8 -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/ThinShearLayer_030915/ShearLayer_rho7.2314_diffRho1e-5_CFL0.05_FD
export OMP_NUM_THREADS=48;bsub -n 48 -W 08:00 ../makefiles/test -minBPD 16 -maxBPD 16 -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/ThinShearLayer_030915/ShearLayer_rho7.2314_diffRho1e-5_CFL0.05_FD
export OMP_NUM_THREADS=48;bsub -n 48 -W 08:00 ../makefiles/test -minBPD 32 -maxBPD 32 -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/ThinShearLayer_030915/ShearLayer_rho7.2314_diffRho1e-5_CFL0.05_FD

cd ../makefiles/;make clean;make bc=periodic config=production multiphase=true precision=single vertexcentered=true poisson=split-fftw particles=true densitydiff=true rk2=false -j;cd -
export OMP_NUM_THREADS=48;bsub -n 48 -W 08:00 ../makefiles/test -minBPD 8 -maxBPD 8 -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/ThinShearLayer_030915/ShearLayer_rho7.2314_diffRho1e-5_CFL0.05_LCFL0.1_particles
export OMP_NUM_THREADS=48;bsub -n 48 -W 08:00 ../makefiles/test -minBPD 16 -maxBPD 16 -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/ThinShearLayer_030915/ShearLayer_rho7.2314_diffRho1e-5_CFL0.05_LCFL0.1_particles
export OMP_NUM_THREADS=48;bsub -n 48 -W 08:00 ../makefiles/test -minBPD 32 -maxBPD 32 -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/ThinShearLayer_030915/ShearLayer_rho7.2314_diffRho1e-5_CFL0.05_LCFL0.1_particles
