cd ../makefiles
make clean
make config=production poisson=hypre bc=mixed particles=false -j
cd ../launch

for BPD in 16 32 64 128
do
	BASENAME=FallingCylinder_Namkoong_long_bpdy${BPD}
	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
	export OMP_NUM_THREADS=48;bsub -n 48 -W 168:00 -o ${BASENAME} -J ${BASENAME} mpirun -np 1 ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}/${BASENAME} -CFL 0.01 -bpdx $[${BPD}/8] -bpdy ${BPD} -tend 50. -rhoS 1.01 -ypos .85 -nu 0.0000044194 -sim falling -shape disk -radius 0.0078125 -tdump 0.1
done
