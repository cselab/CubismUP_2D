for V in true false
do
	for MATTIA in false
#true
	do

		cd ../makefiles
		make clean
		make config=debug poisson=hypre bc=mixed precision=double particles=false constnu=${V} mattia=${MATTIA} -j
		cd ../launch

		for L in 1e6
		do
			for BPD in 4 8 16 32
#4 8 16 32
			do
				NAME=Coquerelle_Mattia${MATTIA}_Nu${V}_Lambda${L}_bpd${BPD}
				FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${NAME}
				mkdir ${FOLDER}
				cp ../makefiles/simulation ${FOLDER}
				export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${NAME} -J ${NAME} mpirun -np 1 ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL 0.1 -LCFL 0.1 -bpdx ${BPD} -bpdy $[${BPD}*3] -radius 0.0625 -tend 20. -rhoS 1.25 -ypos .9 -nu 0.00353553390593274 -sim falling -tdump 0.01 -lambda ${L}
			done
		done

	done
done