for V in true false
do
#	for MATTIA in false true
	for AVG in true
	do
		MATTIA=false
		cd ../makefiles
		make clean
#		make config=production poisson=hypre bc=mixed precision=double particles=false constnu=${V} mattia=${MATTIA} -j
		make config=production poisson=hypre bc=mixed precision=double particles=false constnu=${V} mattia=${MATTIA} avgu=${AVG} -j
		cd ../launch

		for L in 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10
#1e4
		do
			for BPD in 8
#16 32
#4 8 16 32
			do
				NAME=Kolomenskiy_Mattia${MATTIA}_Avgu${AVG}_Nu${V}_Lambda${L}_bpd${BPD}
				FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${NAME}
				mkdir ${FOLDER}
				cp ../makefiles/simulation ${FOLDER}
				export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${NAME} -J ${NAME} mpirun -np 1 ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL 0.1 -LCFL 0.1 -bpdx ${BPD} -bpdy $[${BPD}*4] -radius .0125 -tend 2. -rhoS 2. -ypos .9 -nu 0.000155188776172763 -sim falling -tdump 0.01 -lambda ${L}
			done
		done

	done
done