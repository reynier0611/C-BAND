#!/bin/sh

make clean
make

DATA_PATH="/work/clas12/rg-b/production/recon/pass0/v5/dst/filtered/00"

for run in 6215 6240 6289 6310 6333 6380 6433 6467 6501 6524 6559 6571 6595
do
	for reaction in 1 3 4
	do
		./reactions_BAND $run $reaction $DATA_PATH$run"/dst_inc_00"$run".hipo" "output/outRoot_"$run
	done
done
