cd ../makefiles
make clean
make config=production poisson=split-fftw bc=mixed particles=false rk2=true precision=double -j
cd ../launch

#unclear how it should be run: ratio 1.01 or ratio 1.014?
for BPD in 64
#128 256
#16 32
do
	BASENAME=FallingCylinder_Namkoong_2610_RK2_CFL0.01_Dumps_diff2nd_bpdy${BPD}
	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
    cp Namkoong.sh ${FOLDER}
    cd ${FOLDER}
	export OMP_NUM_THREADS=48;bsub -n 48 -W 40:00 -o ${BASENAME} -J ${BASENAME} ./simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${BASENAME}/${BASENAME} -CFL 0.01 -bpdx $[${BPD}/8] -bpdy ${BPD} -tend 25. -rhoS 1.01 -ypos .85 -nu 0.0000044442057811955 -sim falling -shape disk -radius 0.0078125 -lambda 1e6 -tdump 0.025
    cd -
done
