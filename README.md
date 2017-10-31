#### Running-MAKER-with-little-to-no-evidence
#Construct ab initio gene prediction using only BUSCO augustus models. Add in transcriptome for extra support

my_genome=Sipuncula_Muscle.fasta
###optional transcriptome
my_transcriptome=good.all.orthomerged.fasta
protein_evidence=/mnt/lustre/hcgs/joseph7e/databases/swiss_prot/uniprot_sprot.fasta

### Step 0: run busco with the --long option, this creates a species training model
sbatch quality_check_genome_busco $my_genome
path_to_species_model= 
path_to_proteins=


### Step 0: run Repeat Library Construction http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction--Basic
http://www.repeatmasker.org/RepeatModeler/


### Step 1: Generate Generic control files, http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained
maker -OPTS && maker -BOPTS && maker -EXE

### Step 2: Edit OPTS files
##### Line you will likely change
## genome=
## est= OR altest= (#EST/cDNA sequence file in fasta format from an alternate organism)
## protein= (use swissprot or a set from a closely related species fi yoou have them)
# add genome variable THESE COMMANDS DO NOT WORK WITH ABSOLUE PATHS, JUST MANUALLY EDIT THEM
sed -i 's/'^genome='/'genome="$my_genome"'/g' maker_opts.ctl
# add transcriptome to est line
sed -i 's/'^est='/'est="$my_transcriptome"'/g' maker_opts.ctl
#add protein evidence
### protein=/mnt/lustre/hcgs/joseph7e/databases/swiss_prot/uniprot_sprot.fasta

#add augustus species model




