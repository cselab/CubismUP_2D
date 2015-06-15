cd ../makefiles
make clean
make config=production poisson=hypre bc=mixed particles=false -j
cd ../launch
for BPD in 16
#32 64 128
do
	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_Namkoong_FD_bpdy${BPD}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
	export OMP_NUM_THREADS=48;bsub -n 48 -W 168:00 -o FallingCylinder_Namkoong_FD_bpdy${BPD} mpirun -np 1 ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_Namkoong_FD_bpdy${BPD}/Namkoong_FD_bpdy${BPD} -CFL 0.01 -bpdx $[${BPD}/8] -bpdy ${BPD} -tend 10. -rhoS 1.01 -ypos .85 -nu 0.0000044194 -sim falling -shape disk -radius 0.0078125
done
