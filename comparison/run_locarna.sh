echo Vanilla LocARNA
time mlocarna -v --keep-sequence-order --stockholm --threads=2 --tgtdir locarna_vanilla pesti.fasta
cp locarna_vanilla/results/result.stk pesti_locarna_raw.stk

echo AnchoRNA export+LocARNA
# export anchors to LocARNA format
anchorna export "anchors.gff||8" --fmt locarna -o anchors4locarna.bed
# run LocARCNA with anchor file
time mlocarna -v --keep-sequence-order --stockholm --anchor-constraints=anchors4locarna.bed --threads=2 --tgtdir locarna_anchorna_export pesti.fasta
cp locarna_anchorna_export/results/result.stk pesti_anchorna_export+locarna_raw.stk

echo AnchoRNA cutout+LocARNA
cd ..  # go to root directory
mkdir -p comparison/anchorna_cutout+locarna
time python 02_anchorna_cutout+locarna.py
