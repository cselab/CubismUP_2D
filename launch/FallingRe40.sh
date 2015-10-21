cd ../makefiles
make clean
make config=production poisson=split-fftw bc=mixed precision=double particles=false -j
cd ../launch
for C in 0.05
#0.05 0.1 0.2
do
	for L in 1e6
#1e6 1e7
	do
		for BPD in 8 16 32
#for BPD in 32
		do
			NAME=Falling_DCT_Re40_2110_diff4th_Lambda${L}_bpd${BPD}_CFL${C}
			FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${NAME}
			mkdir ${FOLDER}
			cp ../makefiles/simulation ${FOLDER}
            export OMP_NUM_THREADS=48
           bsub -n 48 -W 10:00 -o ${NAME} -J ${NAME} ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL ${C} -LCFL ${C} -bpdx ${BPD} -bpdy $[${BPD}*4] -radius .025 -tend 2. -rhoS 2. -ypos .9 -nu 0.000846515308813225 -sim falling -tdump 0.025 -lambda ${L}
#bsub -n 48 -W 10:00 -o ${NAME} -J ${NAME} ${FOLDER}/simulation -sim falling -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -restart -serialization /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME}

		done
	done
done


#cd ../makefiles
#make clean
#make config=production poisson=hypre bc=mixed precision=single particles=false rk2=false -j
#cd ../launch
#for C in 0.05
#do
#    for L in 1e6
#    do
#        for BPD in 8 16 32
#        do
#            NAME=Falling_MG_Re40_1910_SP_Lambda${L}_bpd${BPD}_CFL${C}
#            FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${NAME}
#            mkdir ${FOLDER}
#            cp ../makefiles/simulation ${FOLDER}
#            export OMP_NUM_THREADS=48
#            bsub -n 48 -W 10:00 -o ${NAME} -J ${NAME} mpirun -np 32 ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL ${C} -LCFL ${C} -bpdx ${BPD} -bpdy $[${BPD}*4] -radius .025 -tend 3. -rhoS 2. -ypos .9 -nu 0.000846515308813225 -sim falling -tdump 0.025 -lambda ${L}
#        done
#    done
#done
