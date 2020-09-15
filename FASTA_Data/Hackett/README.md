for x in *fasta; do run_seqtools.py -masksites 8 -infile $x -outfile $x.masked8; done
