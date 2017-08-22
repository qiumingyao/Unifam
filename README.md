# Unifam

Metagenome or metaproteome functional annotation by full length gene or protein sequences
1. bin/*: third party tools, i.e. hmmsearch2 and prodigal. These are executable files that can be used directly by unifam.
2. unifam/src/* : the python code for annotating the whole length proteins. This can be run on a single machine.
3. unifam_mm/* : the C code wrapper, for distributing the hmmsearch on multiple computer nodes.
4. unifam_mma/* : the C code wrapper, for distributing the annotation on multiple computer nodes.

##install Unifam and its distributed components:

0. The distrubuted module needs MPI library. For example, module load openmpi/1.8.0
1. copy or clone the repository to your own folder
2. unzip the file folder using tar -xvzf
3. cd unifam_mm/debug; make;
4. cd unifam_mma/debug; make;

##how to use it

1. Run small jobs on single node

python UniFam.py -c configFile -i inputfile

2. Run big jobs

use aprun or mpirun

e.g. aprun -n12800 -d1 -N16 -m 2G command

command:

unifam_mm -ih hmm_model_folder -ip protein_folder -fh hmm_model_summary_file -fp protein_list_summary_file -od mm_output_dir

unifam_mma -hd mm_output_dir -fd protein_folder -data Annotation_database -dtype all/prok -od mma_output_dir

3. Concatinate result:

cat mma_output_dir/*.annot > allresult.annot

cat mma_output_dir/*_annot.faa > allannotated.faa
