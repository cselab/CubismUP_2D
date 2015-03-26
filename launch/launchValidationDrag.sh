module load gcc
export OMP_NUM_THREADS=48

cd ../makefiles/;make clean;make config=production bc=periodic precision=single -j;cd ../launch/

for B in 32 64 128 256
do
	for R in 40 100 1000
	do
		FOLDER_FIXED=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastFixedCylinderRe${R}_$[${B}*32]
		mkdir ${FOLDER_FIXED}
		cp ../makefiles/simulation ${FOLDER_FIXED}
		bsub -n 48 -W 168:00 -o FlowPastFixedCylinderRe${R}_$[${B}*32] ${FOLDER_FIXED}/simulation -file ${FOLDER_FIXED}/FlowPastFixedCylinderRe${R}_$[${B}*32] -bpdx ${B} -bpdy ${B} -radius .025 -uinf .01 -Re ${R} -tend 200 -rhoS 1.00 -sim fixed


		FOLDER_MOVING=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastMovingCylinderRe${R}_$[${B}*32]
		mkdir ${FOLDER_MOVING}
		cp ../makefiles/simulation ${FOLDER_MOVING}
		bsub -n 48 -W 168:00 -o FlowPastMovingCylinderRe${R}_$[${B}*32] ${FOLDER_MOVING}/simulation -file ${FOLDER_MOVING}/FlowPastMovingCylinderRe${R}_$[${B}*32] -bpdx ${B} -bpdy ${B} -radius .025 -uinf .01 -Re ${R} -tend 200 -rhoS 1.00 -sim moving
	done
done