module purge
module load gcc open_mpi
cd ../makefiles
make clean
make config=production poisson=hypre bc=mixed precision=double -j
cd ../launch

for BPD in 32
#64 128
do
	FOLDER_1=/cluster/scratch_xp/public/cconti/CubismUP/Fields_split_bpd${BPD}_ic1
	FOLDER_2=/cluster/scratch_xp/public/cconti/CubismUP/Fields_split_bpd${BPD}_ic2
	FOLDER_3=/cluster/scratch_xp/public/cconti/CubismUP/Fields_split_bpd${BPD}_ic3
	FOLDER_4=/cluster/scratch_xp/public/cconti/CubismUP/Fields_split_bpd${BPD}_ic4
	mkdir ${FOLDER_1}
	mkdir ${FOLDER_2}
	mkdir ${FOLDER_3}
	mkdir ${FOLDER_4}
	cp ../makefiles/simulation ${FOLDER_1}
	cp ../makefiles/simulation ${FOLDER_2}
	cp ../makefiles/simulation ${FOLDER_3}
	cp ../makefiles/simulation ${FOLDER_4}
	export OMP_NUM_THREADS=48;bsub -n 48 -W 168:00 -o Fields_double_hypresplit_bpd${BPD}_ic1_0 -J Fields_double_hypresplit_bpd${BPD}_ic1_0 mpirun -np 32 ${FOLDER_1}/simulation -file ${FOLDER_1}/Fields_double_hypresplit_bpd${BPD}_ic1 -CFL 0.05 -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS 1.1 -ypos .9 -nu 0.0000783022988168292 -sim falling -shape ellipse -semiAxisX .00625 -semiAxisY .025 -angle 0.1 -tdump 0.025 -split
	export OMP_NUM_THREADS=48;bsub -n 48 -W 168:00 -o Fields_double_hypresplit_bpd${BPD}_ic2_0 -J Fields_double_hypresplit_bpd${BPD}_ic2_0 mpirun -np 32 ${FOLDER_2}/simulation -file ${FOLDER_2}/Fields_double_hypresplit_bpd${BPD}_ic2 -CFL 0.05 -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS 1.538 -ypos .9 -nu 0.0000141627381528123 -sim falling -shape ellipse -semiAxisX .005 -semiAxisY .025 -angle 0.1 -tdump 0.025 -split
	export OMP_NUM_THREADS=48;bsub -n 48 -W 168:00 -o Fields_double_hypresplit_bpd${BPD}_ic3_0 -J Fields_double_hypresplit_bpd${BPD}_ic3_0 mpirun -np 32 ${FOLDER_3}/simulation -file ${FOLDER_3}/Fields_double_hypresplit_bpd${BPD}_ic3 -CFL 0.05 -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS 1.1 -ypos .9 -nu 0.00000783022988168292 -sim falling -shape ellipse -semiAxisX .00625 -semiAxisY .025 -angle 0.1 -tdump 0.025 -split
	export OMP_NUM_THREADS=48;bsub -n 48 -W 168:00 -o Fields_double_hypresplit_bpd${BPD}_ic4_0 -J Fields_double_hypresplit_bpd${BPD}_ic4_0 mpirun -np 32 ${FOLDER_4}/simulation -file ${FOLDER_4}/Fields_double_hypresplit_bpd${BPD}_ic4 -CFL 0.05 -bpdx ${BPD} -bpdy ${BPD} -tend 200. -rhoS 1.248 -ypos .9 -nu 0.00000985240110463296 -sim falling -shape ellipse -semiAxisX .0125 -semiAxisY .025 -angle 0.1 -tdump 0.025 -split

	for N in {1..5}
	do
		export OMP_NUM_THREADS=48;bsub -n 48 -W 168:00 -w "ended(Fields_double_hypresplit_bpd${BPD}_ic1_$[${N}-1])" -o Fields_double_hypresplit_bpd${BPD}_ic1_$N -J Fields_double_hypresplit_bpd${BPD}_ic1_$N mpirun -np 32 ${FOLDER_1}/simulation -serialization ${FOLDER_1}/Fields_double_hypresplit_bpd${BPD}_ic1 -restart -sim falling
		export OMP_NUM_THREADS=48;bsub -n 48 -W 168:00 -w "ended(Fields_double_hypresplit_bpd${BPD}_ic2_$[${N}-1])" -o Fields_double_hypresplit_bpd${BPD}_ic2_$N -J Fields_double_hypresplit_bpd${BPD}_ic2_$N mpirun -np 32 ${FOLDER_2}/simulation -serialization ${FOLDER_1}/Fields_double_hypresplit_bpd${BPD}_ic2 -restart -sim falling
		export OMP_NUM_THREADS=48;bsub -n 48 -W 168:00 -w "ended(Fields_double_hypresplit_bpd${BPD}_ic3_$[${N}-1])" -o Fields_double_hypresplit_bpd${BPD}_ic3_$N -J Fields_double_hypresplit_bpd${BPD}_ic3_$N mpirun -np 32 ${FOLDER_3}/simulation -serialization ${FOLDER_1}/Fields_double_hypresplit_bpd${BPD}_ic3 -restart -sim falling
		export OMP_NUM_THREADS=48;bsub -n 48 -W 168:00 -w "ended(Fields_double_hypresplit_bpd${BPD}_ic4_$[${N}-1])" -o Fields_double_hypresplit_bpd${BPD}_ic4_$N -J Fields_double_hypresplit_bpd${BPD}_ic4_$N mpirun -np 32 ${FOLDER_4}/simulation -serialization ${FOLDER_1}/Fields_double_hypresplit_bpd${BPD}_ic4 -restart -sim falling
	done
done