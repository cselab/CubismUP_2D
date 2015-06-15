#make clean;make -j
#make clean;make bc=periodic -j
#make clean;make precision=double -j
#make clean;make poisson=split-fftw -j
#make clean;make poisson=hypre -j
#make clean;make multiphase=true

#make clean;make config=production -j
#make clean;make config=production bc=periodic -j
#make clean;make config=production precision=double -j
#make clean;make config=production poisson=split-fftw
#make clean;make config=production multiphase=true

for CONFIG in debug production
do
	for BC in periodic mixed pipe vortex
	do
		for PRECISION in single double
		do
			for POISSON in split-fftw hypre
			do
				for MULTIPHASE in true false
				do
					for PARTICLES in true false
					do
						for VERTEX in true false
						do
							make clean
							make config=${CONFIG} bc=${BC} precision=${PRECISION} multiphase=${MULTIPHASE} particles=${PARTICLES} vertexcentered=${VERTEX} -j
						done
					done
				done
			done
		done
	done
done