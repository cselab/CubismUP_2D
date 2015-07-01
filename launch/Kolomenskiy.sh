cd ../makefiles
make clean
make config=debug poisson=hypre bc=mixed precision=double particles=false -j
cd ../launch

for BPD in 16 32 64
do
	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_Kolomenskiy_bpd${BPD}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
	export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o Kolomenskiy_${BPD} -J Kolomenskiy_${BPD} mpirun -np 32 ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/FallingCylinder_Kolomenskiy_bpd${BPD}/FallingCylinder_Kolomenskiy_bpd${BPD} -CFL 0.1 -LCFL 0.1 -bpdx ${BPD} -bpdy ${BPD} -radius .0125 -tend 2. -rhoS 2. -ypos .9 -nu 0.000155188776172763 -sim falling -tdump 0.01
done