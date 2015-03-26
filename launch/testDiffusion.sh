for BS in 32 64 128 256 512 1024 2048
do
	cd ../makefiles
	make clean
	make config=production poisson=hypre bc=periodic precision=double bs=${BS} -j
	cd ../launch
	../makefiles/test -minBPD 1 -maxBPD 1 -test diffusion > diff_${BS}
done

for BS in 32 64 128 256 512 1024 2048
do
	cat diff_${BS}
done
