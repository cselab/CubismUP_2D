export OMP_NUM_THREADS=48

BASEPATH=/cluster/scratch_xp/public/cconti/CubismUP/Thesis/
BASENAME=OscillatingCylinder_24012016_CFL0.01

cd ../makefiles/
make clean
make config=production bc=box poisson=split-fftw precision=double particles=false dlm=true -j
cd ../launch/

for R in 0.05 0.005
do
for B in 16 32 64
do
	OPTIONS=" -bpdx ${B} -bpdy ${B} -CFL 0.01 -radius ${R} -rhoS 1. -dlm 1"

	NAME=_r${R}_$[${B}*32]

	FOLDER=${BASEPATH}${BASENAME}${NAME}
	mkdir ${FOLDER}
	cp ../makefiles/simulation ${FOLDER}
	cp oscillatingCylinder.sh ${FOLDER}
	cd ${FOLDER}
	bsub -n 48 -W 10:00 -o ${BASENAME}${NAME} -J ${BASENAME}${NAME} ./simulation -file ${BASENAME}${NAME} -sim oscillating -Re 100 -kc 5 -tend 10 ${OPTIONS} -fdump 100
	cd -
done
done