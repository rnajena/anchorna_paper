echo Vanilla MAFFT
time mafft --auto pesti.fasta > pesti_mafft.aln
#anchorna view anchors50.gff --fname pesti_mafft.aln --mode nt
echo AnchoRNA+MAFFT
anchorna cutout anchors50.gff start "A14<" --fname pesti.fasta -o pesti4mafft_start_A14.fasta
anchorna cutout anchors50.gff "A14<" end --fname pesti.fasta -o pesti4mafft_A14_end.fasta
time mafft --auto pesti4mafft_start_A14.fasta > pesti4mafft_start_A14.aln
time mafft --auto pesti4mafft_A14_end.fasta > pesti4mafft_A14_end.aln
sugar merge pesti4mafft_start_A14.aln pesti4mafft_A14_end.aln -o pesti_anchorna+mafft.aln --fmtout fasta
#anchorna view anchors50.gff --fname pesti_anchorna+mafft.aln --mode nt
