cd ../makefiles
make clean
make config=production poisson=split-fftw bc=mixed precision=single particles=false rk2=true dlm=true -j
cd ../launch

FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/Andersen_ParamStudy_DLM_11032016_J0_Stats
mkdir ${FOLDER}
cp ../makefiles/simulation ${FOLDER}/simulation
cp Andersen_ParameterStudy.sh ${FOLDER}

for BPD in 32
# 8 16
do
for CFL in 0.01 0.1
#0.01
do
for RHO in 2.7
#1.35
#1.25 1.35 1.45 2.6 2.7 2.8 4.0 5.0
do
for X in 0.5
#0.3
do
	BASENAME=Andersen_ParamStudy_11032016_CFL${CFL}_bpd${BPD}_rho${RHO}_xpos${X}_J0
	cd ${FOLDER}

	OPTIONS=" -CFL ${CFL} -bpdx ${BPD} -bpdy ${BPD} -tend 10. -rhoS ${RHO} -angle 0.2 -sim falling -shape ellipse -xpos ${X} -ypos .9"
#export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME}_L -J ${BASENAME}_L ./simulation -file ${FOLDER}/${BASENAME}_L ${OPTIONS} -nu 0.000152606310013717 -semiAxisX 0.0125 -semiAxisY 0.1
#export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME}_M -J ${BASENAME}_M ./simulation -file ${FOLDER}/${BASENAME}_M ${OPTIONS} -nu 0.0000539544783312781 -semiAxisX 0.00625 -semiAxisY 0.05
	export OMP_NUM_THREADS=48;bsub -n 48 -W 10:00 -o ${BASENAME}_S -J ${BASENAME}_S ./simulation -file ${FOLDER}/${BASENAME}_S ${OPTIONS} -nu 0.0000190757887517147 -semiAxisX 0.003125 -semiAxisY 0.025
	cd -
done
done
done
done