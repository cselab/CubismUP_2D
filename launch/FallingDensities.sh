cd ../makefiles
make clean
make config=production poisson=split-fftw bc=mixed precision=single rk2=true particles=false -j
cd ../launch
for C in 0.05
do
    for BPD in 16
    do
#        for D in 1.1 2.0 10. 100. 1000.
#        do
#            NAME=Falling_2010_Lambda1e6_bpd${BPD}_CFL${C}_rhos${D}
#            FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${NAME}
#            mkdir ${FOLDER}
#            cp ../makefiles/simulation ${FOLDER}
#            export OMP_NUM_THREADS=48
#            bsub -n 48 -W 10:00 -o ${NAME} -J ${NAME} ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL ${C} -bpdx ${BPD} -bpdy $[${BPD}*4] -radius .025 -tend 10. -rhoS ${D} -ypos .9 -nu 0.000846515308813225 -sim falling -tdump 0.025 -lambda 1e6
#        done
        for D in 0.01 0.1
#0.5 0.9
        do
            NAME=Rising_2110_Lambda1e6_bpd${BPD}_CFL${C}_rhos${D}
            FOLDER=/cluster/scratch_xp/public/cconti/CubismUP/${NAME}
            mkdir ${FOLDER}
            cp ../makefiles/simulation ${FOLDER}
            export OMP_NUM_THREADS=48
            bsub -n 48 -W 168:00 -o ${NAME} -J ${NAME} ${FOLDER}/simulation -file /cluster/scratch_xp/public/cconti/CubismUP/${NAME}/${NAME} -CFL ${C} -bpdx ${BPD} -bpdy $[${BPD}*4] -radius .025 -tend 10. -rhoS ${D} -ypos .1 -nu 0.000846515308813225 -sim falling -tdump 0.025 -lambda 1e6
        done
	done
done
