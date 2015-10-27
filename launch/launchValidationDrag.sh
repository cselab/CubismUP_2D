#module load gcc
export OMP_NUM_THREADS=48

cd ../makefiles/;make clean;make config=production bc=periodic precision=single particles=false -j;cd ../launch/

for B in 32
#128 256
#32 64 128 256
do
	# Re 40
	FOLDER_FIXED=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastFixedCylinderRe40_2110_lambda1e6_$[${B}*32]
	mkdir ${FOLDER_FIXED}
	cp ../makefiles/simulation ${FOLDER_FIXED}
	bsub -n 48 -W 10:00 -o FlowPastFixedCylinderRe40_$[${B}*32] ${FOLDER_FIXED}/simulation -file ${FOLDER_FIXED}/FlowPastFixedCylinderRe40_$[${B}*32] -bpdx ${B} -bpdy ${B} -CFL 0.01 -radius .025 -uinf .01 -Re 40 -tend 8 -rhoS 1.00 -sim fixed -lambda 1e6


	FOLDER_MOVING=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastMovingCylinderRe40_2110_lambda1e6_$[${B}*32]
	mkdir ${FOLDER_MOVING}
	cp ../makefiles/simulation ${FOLDER_MOVING}
	bsub -n 48 -W 10:00 -o FlowPastMovingCylinderRe40_$[${B}*32] ${FOLDER_MOVING}/simulation -file ${FOLDER_MOVING}/FlowPastMovingCylinderRe40_$[${B}*32] -bpdx ${B} -bpdy ${B} -CFL 0.01 -radius .025 -uinf .01 -Re 40 -tend 8 -rhoS 1.00 -sim moving -lambda 1e6

	# Re 100
	FOLDER_FIXED=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastFixedCylinderRe100_2110_lambda1e6_$[${B}*32]
	mkdir ${FOLDER_FIXED}
	cp ../makefiles/simulation ${FOLDER_FIXED}
	bsub -n 48 -W 10:00 -o FlowPastFixedCylinderRe100_$[${B}*32] ${FOLDER_FIXED}/simulation -file ${FOLDER_FIXED}/FlowPastFixedCylinderRe100_$[${B}*32] -bpdx ${B} -bpdy ${B} -CFL 0.01 -radius .025 -uinf .01 -Re 100 -tend 30 -rhoS 1.00 -sim fixed -lambda 1e6


	FOLDER_MOVING=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastMovingCylinderRe100_2110_lambda1e6_$[${B}*32]
	mkdir ${FOLDER_MOVING}
	cp ../makefiles/simulation ${FOLDER_MOVING}
	bsub -n 48 -W 10:00 -o FlowPastMovingCylinderRe100_$[${B}*32] ${FOLDER_MOVING}/simulation -file ${FOLDER_MOVING}/FlowPastMovingCylinderRe100_$[${B}*32] -bpdx ${B} -bpdy ${B} -CFL 0.01 -radius .025 -uinf .01 -Re 100 -tend 30 -rhoS 1.00 -sim moving -lambda 1e6

	# Re 1000
	FOLDER_FIXED=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastFixedCylinderRe1000_2110_lambda1e6_$[${B}*32]
	mkdir ${FOLDER_FIXED}
	cp ../makefiles/simulation ${FOLDER_FIXED}
	bsub -n 48 -W 10:00 -o FlowPastFixedCylinderRe1000_$[${B}*32] ${FOLDER_FIXED}/simulation -file ${FOLDER_FIXED}/FlowPastFixedCylinderRe1000_$[${B}*32] -bpdx ${B} -bpdy ${B} -CFL 0.01 -radius .025 -uinf .01 -Re 1000 -tend 5 -rhoS 1.00 -sim fixed -lambda 1e6


	FOLDER_MOVING=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastMovingCylinderRe1000_2110_lambda1e6_$[${B}*32]
	mkdir ${FOLDER_MOVING}
	cp ../makefiles/simulation ${FOLDER_MOVING}
	bsub -n 48 -W 10:00 -o FlowPastMovingCylinderRe1000_$[${B}*32] ${FOLDER_MOVING}/simulation -file ${FOLDER_MOVING}/FlowPastMovingCylinderRe1000_$[${B}*32] -bpdx ${B} -bpdy ${B} -CFL 0.01 -radius .025 -uinf .01 -Re 1000 -tend 5 -rhoS 1.00 -sim moving -lambda 1e6
done
