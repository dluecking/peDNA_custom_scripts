for i in top20/gtas/interpro/input_proteins/*
do
	f="$(basename -- ${i})"
	echo "working on: "
	echo $f
	python ~/bin/webservice-clients/python/iprscan5.py --sequence=${i}  \
	--outfile=top20/gtas/interpro/interpro_out/${f} --email=dlueckin@mpi-bremen.de
done
