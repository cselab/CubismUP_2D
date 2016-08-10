for p in false
do
for rk in false
do
	cd ../makefiles/
	make clean
	make bc=periodic config=production multiphase=false precision=single vertexcentered=true poisson=fftw particles=${p} densitydiff=false vardensity=false rk2=${rk} -j
	cd -

	BASENAME=Thesis_ThinShearLayer_24022016_Particles${p}_RK2${rk}_constRho
	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Thesis/${BASENAME}
	mkdir ${FOLDER}

	for BPD in 4 8 16 32
	do
		cp ../makefiles/test ${FOLDER}/test_BPD${BPD}_Particles${P}_RK2${rk}
		cd ${FOLDER}
		export OMP_NUM_THREADS=48;bsub -n 48 -W 08:00 -o ${BASENAME}_BPD${BPD}_Particles${P}_RK2${rk} -J ${BASENAME}_BPD${BPD}_Particles${P}_RK2${rk} ./test_BPD${BPD}_Particles${P}_RK2${rk} -minBPD ${BPD} -maxBPD ${BPD} -test shearlayer -file ShearLayer_CFL0.1_constRho_BPD${BPD} -CFL 0.1
		cd -
	done


#	cd ../makefiles/
#	make clean
#	make bc=periodic config=production multiphase=true precision=single vertexcentered=true poisson=split-fftw particles=${p} densitydiff=true vardensity=true rk2=${rk} -j
#	cd -

#	BASENAME=Thesis_ThinShearLayer_20012016_Particles${p}_RK2${rk}_variableRho
#	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Thesis/${BASENAME}
#	mkdir ${FOLDER}

#	for BPD in 4 8 16 32
#	do
#		cp ../makefiles/test ${FOLDER}/test_BPD${BPD}_Particles${P}_RK2${rk}
#		cd ${FOLDER}
#		export OMP_NUM_THREADS=48;bsub -n 48 -W 08:00 -o ${BASENAME}_BPD${BPD}_Particles${P}_RK2${rk} -J ${BASENAME}_BPD${BPD}_Particles${P}_RK2${rk} ./test_BPD${BPD}_Particles${P}_RK2${rk} -minBPD ${BPD} -maxBPD ${BPD} -test shearlayer -file ShearLayer_CFL0.01_variableRho_BPD${BPD} -CFL 0.01
#		cd -
#	done
done
done
