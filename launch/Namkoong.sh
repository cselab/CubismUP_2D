cd ../makefiles
make clean
make config=production poisson=split-fftw bc=mixed particles=false rk2=true precision=double -j
cd ../launch

#unclear how it should be run: ratio 1.01 or ratio 1.014?
for BPD in 16 32 64 128 256
do
	BASENAME=FallingCylinder_Namkoong_2110_RK2_CFL0.1_diff4th_bpdy${BPD}
	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
	export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME} -J ${BASENAME} ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}/${BASENAME} -CFL 0.1 -bpdx $[${BPD}/8] -bpdy ${BPD} -tend 30. -rhoS 1.01 -ypos .85 -nu 0.0000044442057811955 -sim falling -shape disk -radius 0.0078125 -lambda 1e6
done
