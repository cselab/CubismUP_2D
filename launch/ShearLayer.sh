for P in false true
#for P in true
do
	cd ../makefiles/
	make clean
	make bc=periodic config=production multiphase=true precision=single vertexcentered=true poisson=split-fftw constnu=false particles=${P} -j
#make bc=periodic config=production multiphase=true precision=single vertexcentered=true poisson=hypre particles=${P} -j
	cd -
	export OMP_NUM_THREADS=48

	for B in 8 16 32
	do
		#bsub -W 168:00 -n 48 ../makefiles/test -minBPD ${B} -maxBPD ${B} -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/ThinShearLayer_210815/ShearLayer${B}_rho7.2314_hypre_CFL0.01_particles${P}
		bsub -W 10:00 -n 48 ../makefiles/test -minBPD ${B} -maxBPD ${B} -test shearlayer -file /cluster/scratch_xp/public/cconti/CubismUP/ThinShearLayer_210815/ShearLayer${B}_rho7.2314_fftw_diffRho1e-6_CFL0.1_particles${P}
		#../makefiles/test -minBPD 1 -maxBPD 16 -test shearlayer -file ../data/testShearLayer_Convergence_
	done
done
