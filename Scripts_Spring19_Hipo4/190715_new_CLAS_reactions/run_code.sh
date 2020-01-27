#!/bin/sh

make clean
make

#DATA_PATH="/cache/clas12/rg-b/production/recon/pass0/v16/dst/filtered/00"

for run in 6233 6433 6524 6546 6559 6571 6586 6595
do
	./reactions $run 4 "/cache/clas12/rg-b/production/recon/pass0/v16/dst/filtered/00"$run"/dst_inc_00"$run".hipo"
done

#6233 6302 6303 6305 6307 6310 6311 6313 6321 6326 6327 6328 6346 6347 6349 6420 6428 6433 6442 6450 6467 6474 6481 6492 6501 6502 6515 6522 6524 6546 6559 6571 6586 6595
