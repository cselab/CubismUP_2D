#module load gcc
export OMP_NUM_THREADS=48

cd ../makefiles/;make clean;make config=production bc=periodic precision=single particles=false -j;cd ../launch/

for B in 256
#32 64 128 256
do
	# Re 40
#	FOLDER_FIXED=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastFixedCylinderRe40_$[${B}*32]
#	mkdir ${FOLDER_FIXED}
#	cp ../makefiles/simulation ${FOLDER_FIXED}
#	bsub -n 48 -W 168:00 -o FlowPastFixedCylinderRe40_$[${B}*32] ${FOLDER_FIXED}/simulation -file ${FOLDER_FIXED}/FlowPastFixedCylinderRe40_$[${B}*32] -bpdx ${B} -bpdy ${B} -CFL 0.01 -radius .025 -uinf .01 -Re 40 -tend 8 -rhoS 1.00 -sim fixed -restart -serialization ${FOLDER_FIXED}/FlowPastFixedCylinderRe40_$[${B}*32]


#	FOLDER_MOVING=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastMovingCylinderRe40_$[${B}*32]
#	mkdir ${FOLDER_MOVING}
#	cp ../makefiles/simulation ${FOLDER_MOVING}
#	bsub -n 48 -W 168:00 -o FlowPastMovingCylinderRe40_$[${B}*32] ${FOLDER_MOVING}/simulation -file ${FOLDER_MOVING}/FlowPastMovingCylinderRe40_$[${B}*32] -bpdx ${B} -bpdy ${B} -CFL 0.01 -radius .025 -uinf .01 -Re 40 -tend 8 -rhoS 1.00 -sim moving -restart -serialization ${FOLDER_FIXED}/FlowPastMovingCylinderRe40_$[${B}*32]

	# Re 100
	FOLDER_FIXED=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastFixedCylinderRe100_$[${B}*32]
	mkdir ${FOLDER_FIXED}
	cp ../makefiles/simulation ${FOLDER_FIXED}
	bsub -n 48 -W 168:00 -o FlowPastFixedCylinderRe100_$[${B}*32] ${FOLDER_FIXED}/simulation -file ${FOLDER_FIXED}/FlowPastFixedCylinderRe100_$[${B}*32] -bpdx ${B} -bpdy ${B} -CFL 0.01 -radius .025 -uinf .01 -Re 100 -tend 30 -rhoS 1.00 -sim fixed -restart -serialization ${FOLDER_FIXED}/FlowPastFixedCylinderRe100_$[${B}*32]


	FOLDER_MOVING=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastMovingCylinderRe100_$[${B}*32]
	mkdir ${FOLDER_MOVING}
	cp ../makefiles/simulation ${FOLDER_MOVING}
	bsub -n 48 -W 168:00 -o FlowPastMovingCylinderRe100_$[${B}*32] ${FOLDER_MOVING}/simulation -file ${FOLDER_MOVING}/FlowPastMovingCylinderRe100_$[${B}*32] -bpdx ${B} -bpdy ${B} -CFL 0.01 -radius .025 -uinf .01 -Re 100 -tend 30 -rhoS 1.00 -sim moving -restart -serialization ${FOLDER_MOVING}/FlowPastMovingCylinderRe100_$[${B}*32]

	# Re 1000
#	FOLDER_FIXED=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastFixedCylinderRe1000_$[${B}*32]
#	mkdir ${FOLDER_FIXED}
#	cp ../makefiles/simulation ${FOLDER_FIXED}
#	bsub -n 48 -W 168:00 -o FlowPastFixedCylinderRe1000_$[${B}*32] ${FOLDER_FIXED}/simulation -file ${FOLDER_FIXED}/FlowPastFixedCylinderRe1000_$[${B}*32] -bpdx ${B} -bpdy ${B} -CFL 0.01 -radius .025 -uinf .01 -Re 1000 -tend 5 -rhoS 1.00 -sim fixed -restart -serialization ${FOLDER_FIXED}/FlowPastFixedCylinderRe1000_$[${B}*32]


#	FOLDER_MOVING=/cluster/scratch_xp/public/cconti/CubismUP/Drag_FlowPastMovingCylinderRe1000_lambda1e8_$[${B}*32]
#	mkdir ${FOLDER_MOVING}
#	cp ../makefiles/simulation ${FOLDER_MOVING}
#	bsub -n 48 -W 168:00 -o FlowPastMovingCylinderRe1000_$[${B}*32] ${FOLDER_MOVING}/simulation -file ${FOLDER_MOVING}/FlowPastMovingCylinderRe1000_$[${B}*32] -bpdx ${B} -bpdy ${B} -CFL 0.01 -radius .025 -uinf .01 -Re 1000 -tend 5 -rhoS 1.00 -sim moving -lambda 1e8
done
