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
							for DLM in true false
							do
								make clean
								make config=${CONFIG} bc=${BC} precision=${PRECISION} multiphase=${MULTIPHASE} particles=${PARTICLES} vertexcentered=${VERTEX} dlm=${DLM} -j
							done
						done
					done
				done
			done
		done
	done
done