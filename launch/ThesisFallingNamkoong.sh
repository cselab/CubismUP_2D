cd ../makefiles
make clean
make config=production poisson=split-fftw bc=mixed particles=false rk2=true precision=double dlm=true -j
cd ../launch

#unclear how it should be run: ratio 1.01 or ratio 1.014?
for D in 1 0.1 10
do
for R in 1.01
do
for CFL in 0.1 0.01
do
for BPD in 64
do
	BASENAME=FallingCylinder_Namkoong_24012016_DLM${D}_RK2_CFL${CFL}_bpdy${BPD}_R${R}
	FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Thesis/${BASENAME}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
    cp ThesisFallingNamkoong.sh ${FOLDER}
    cd ${FOLDER}
	export OMP_NUM_THREADS=48
	bsub -n 48 -W 10:00 -o ${BASENAME} -J ${BASENAME} ./simulation -file ${FOLDER}/${BASENAME} -CFL ${CFL} -bpdx $[${BPD}/8] -bpdy ${BPD} -tend 25. -rhoS ${R} -ypos .85 -nu 0.0000044442057811955 -sim falling -shape disk -radius 0.0078125 -dlm ${D}
    cd -
done
for BPD in 128
do
BASENAME=FallingCylinder_Namkoong_24012016_DLM${D}_RK2_CFL${CFL}_bpdy${BPD}_R${R}
FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Thesis/${BASENAME}
mkdir ${FOLDER}
cp ../makefiles/simulation ${FOLDER}
cp ThesisFallingNamkoong.sh ${FOLDER}
cd ${FOLDER}
export OMP_NUM_THREADS=48
bsub -n 48 -W 40:00 -o ${BASENAME} -J ${BASENAME} ./simulation -file ${FOLDER}/${BASENAME} -CFL ${CFL} -bpdx $[${BPD}/8] -bpdy ${BPD} -tend 25. -rhoS ${R} -ypos .85 -nu 0.0000044442057811955 -sim falling -shape disk -radius 0.0078125 -dlm ${D}
cd -
done
for BPD in 256
do
BASENAME=FallingCylinder_Namkoong_24012016_DLM${D}_RK2_CFL${CFL}_bpdy${BPD}_R${R}
FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Thesis/${BASENAME}
mkdir ${FOLDER}
cp ../makefiles/simulation ${FOLDER}
cp ThesisFallingNamkoong.sh ${FOLDER}
cd ${FOLDER}
export OMP_NUM_THREADS=48
bsub -n 48 -W 168:00 -o ${BASENAME} -J ${BASENAME} ./simulation -file ${FOLDER}/${BASENAME} -CFL ${CFL} -bpdx $[${BPD}/8] -bpdy ${BPD} -tend 25. -rhoS ${R} -ypos .85 -nu 0.0000044442057811955 -sim falling -shape disk -radius 0.0078125 -dlm ${D}
cd -
done
done
done
done