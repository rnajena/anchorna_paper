# Code repository for the AnchoRNA publication

This repository contains the source code for the reproduction of results and figures of the following publication:

Tom Eulenfeld, Sandra Triebel, Peter F. Stadler, Manja Marz (2025),
AnchoRNA: Full virus genome alignments through conserved anchor regions,
doi:[10.1101/2025.01.30.635689](https://doi.org/10.1101/2025.01.30.635689)
[[pdf](https://www.biorxiv.org/content/10.1101/2025.01.30.635689.full.pdf)]

### Create a dedicated conda environment

We also report the versions of the relevant packages.

```sh
conda create -n anchorna_paper -c conda-forge -c bioconda python==3.12 mafft==7.526 locarna=2.0
conda activate anchorna_paper
pip install rnajena-sugar==0.8 anchorna==2.0
```

After creating the environment, you can run all of the sections below, or just some of them.
Files in this repository will be recreated, except for these files:
The hand-curated alignment file `pesti_hand_curated_55seqs_raw.stk`,
the sequence IDs for the comparison `ids_50.txt`,
the IDs of the sequences shown in figure 5 `ids_24.txt`,
and the Python and BASH scripts.

### Search anchors

First, we copy the pesti data file with 55 representative sequences,
it is included in the AnchoRNA distribution.

```sh
cd comparison
anchorna create --tutorial
```

Now we  the `run 01_prepare.py` script to create the files with 50 out of 55 sequences.

```sh
cd ..; python 01_prepare.py; cd comparison
```

We are still in the comparison folder.
Search anchors for 55 sequences and for 50 sequences:

```sh
anchorna go anchors.gff
anchorna go --fname pesti.gff anchors50.gff
anchorna export anchors.gff --mode nt -o anchors_nt.gff  # convert aa anchors to nucleotides
```

Create the screen shots for figure 3:

```sh
anchorna view anchors.gff --align A1
anchorna cutout anchors.gff start "A1+6" | anchorna view anchors.gff --mode nt --align A0 --fname -
```

### Create MAFFT alignments with version v7.526

```sh
echo Vanilla MAFFT
time mafft --auto pesti.fasta > pesti_mafft.aln
anchorna view anchors50.gff --fname pesti_mafft.aln --mode nt

echo AnchoRNA+MAFFT
anchorna cutout anchors50.gff start "A14<" --fname pesti.fasta -o pesti4mafft_start_A14.fasta
anchorna cutout anchors50.gff "A14<" end --fname pesti.fasta -o pesti4mafft_A14_end.fasta
time mafft --auto pesti4mafft_start_A14.fasta > pesti4mafft_start_A14.aln
time mafft --auto pesti4mafft_A14_end.fasta > pesti4mafft_A14_end.aln
sugar merge pesti4mafft_start_A14.aln pesti4mafft_A14_end.aln -o pesti_anchorna+mafft.aln --fmtout fasta
anchorna view anchors50.gff --fname pesti_anchorna+mafft.aln --mode nt
```

To log the output you can also run the script `bash run_mafft.sh 2&> run_mafft_log.txt`.

### Create Dialign alignments using the server at https://dialign.gobics.de/

* Use Dialign to create pesti_dialign_uncorrected.aln
* Use CHAOS+Dialign to create pesti_chaos+dialign_uncorrected.aln
* Create anchor file for Dialign
  ```sh
  anchorna export "anchors50.gff|13,14,15,40,41,42" -m nt --fmt dialign --fname pesti.gff -o anchors4dialign.txt
  ```
* Use Anchored Dialign (AnchoRNA+Dialign) to create `pesti_anchorna+dialign_uncorrected.aln`

### Create Locarna alignments with LocARNA 2.0.1.

```sh
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
```

To log the output you can also run the script `bash run_locarna.sh 2&> run_locarna_log.txt`.

### Run the remaining analysis and plotting scripts

Be sure, that you are in the root directory.

```sh
python 03_postprocess_alis.py  # clean up Dialign files, postprocess LocARNA files
python 04_alignment_scores.py  # calculate scores for table 1
python 05_plots.py  # create plots for figures 4, 5, 6
```
