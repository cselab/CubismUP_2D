cd ../makefiles
make clean
make config=production poisson=hypre bc=mixed precision=single particles=false -j
cd ../launch
for C in 0.05 0.1 0.2
do
	for L in 1e6 1e7
	do
		for BPD in 8 16 32
		do
			NAME=Falling_Re40_0209_reordered_Lambda${L}_bpd${BPD}_CFL${C}
			FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${NAME}
			mkdir ${FOLDER}
			cp ../makefiles/simulation ${FOLDER}
			export OMP_NUM_THREADS=48
			bsub -n 48 -W 168:00 -o ${NAME} -J ${NAME} mpirun -np 32 ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL ${C} -LCFL ${C} -bpdx ${BPD} -bpdy $[${BPD}*4] -radius .025 -tend 10 -rhoS 2. -ypos .9 -nu 0.0006 -sim falling -tdump 0.01 -lambda ${L}
		done
	done
done