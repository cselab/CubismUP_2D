module load gcc
cd ../makefiles
make clean
make config=production poisson=hypre bc=mixed -j
cd ../launch
for BPD in 32 128
do
	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_Namkoong_bpdy${BPD}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
	export OMP_NUM_THREADS=48;bsub -n 48 -W 168:00 -o FallingCylinder_Namkoong_bpdy${BPD} mpirun -np 32 ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_Namkoong_bpdy${BPD}/Namkoong_bpdy${BPD} -CFL 0.01 -bpdx $[${BPD}/4] -bpdy ${BPD} -tend 10. -rhoS 1.01 -ypos .85 -nu 0.0000044194 -sim falling -shape disk -radius 0.0078125
done