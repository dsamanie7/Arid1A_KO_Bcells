##### MAP with STAR


##### Do the mapping in one for loop instead of submitting several jobs at the time.

#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=12:00:00
#SBATCH --job-name=star_ajay11
#SBATCH --output=/projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/scripts/star_50_vM23_mm10_star_loop_PBS.out
#SBATCH --error=/projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/scripts/star_50_vM23_mm10_star_loop_PBS.err

cd /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro

module load bowtie
module load samtools/1.10.1
module load java/jdk1.8.0_191
module load picard/2.21.4
module load deeptools

for i in $(cat /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/index/map.txt)
do
echo ${i} > TEMP
NAME=$(cut -f1 -d: TEMP)


STAR --runThreadN 24 -- genomeDir  /projects/p31767/Bioinformatics/mouse/star/star_50_vM23/ --readFilesIn /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/fastq/${NAME}_R1_001.fastq /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/fastq/${NAME}_R2_001.fastq --genomeLoad LoadAndRemove --outFilterMismatchNmax 4 --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100 --genomeFileSizes /projects/p31767/Bioinformatics/mouse/star/star_50/chrNameLength.txt --outFileNamePrefix /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/star/${NAME}- 

done


##### SORT AND REMOVE PCR DUPLICATES

for i in $(cat /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/index/map.txt)
do
echo ${i} > TEMP
NAME=$(cut -f1 -d: TEMP)

cat <<EOT >> /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/scripts/star_50_vM23mm10_map_${NAME}_PBS.sh
#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=12:00:00
#SBATCH --job-name=star_Ajay_11
#SBATCH --output=/projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/scripts/star_50_vM23mm10_map_${NAME}_PBS.out
#SBATCH --error=/projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/scripts/star_50_vM23mm10_map_${NAME}_PBS.err

cd /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro

module load bowtie
module load samtools/1.10.1
module load java/jdk1.8.0_191
module load picard/2.21.4
module load deeptools


samtools view -bS /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/star/${NAME}-Aligned.out.sam > /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/star/${NAME}-Aligned.out.bam
samtools sort /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/star/${NAME}-Aligned.out.bam -o /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/star/${NAME}-Aligned.out.sorted.bam
picard MarkDuplicates I=/projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/star/${NAME}-Aligned.out.sorted.bam O=/projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/star/${NAME}-Aligned.out.sorted.PCR.bam M=metrics.txt REMOVE_DUPLICATES=true
samtools index /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/star/${NAME}-Aligned.out.sorted.PCR.bam
bamCoverage --samFlagInclude 64 --filterRNAstrand forward -b /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/star/${NAME}-Aligned.out.sorted.PCR.bam --normalizeUsing CPM -o /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/star/${NAME}-Aligned.out.sorted.PCR.forward.bw
bamCoverage --samFlagInclude 64 --filterRNAstrand reverse -b /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/star/${NAME}-Aligned.out.sorted.PCR.bam --normalizeUsing CPM -o /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/star/${NAME}-Aligned.out.sorted.PCR.reverse.bw

EOT

sbatch /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/scripts/star_50_vM23mm10_map_${NAME}_PBS.sh

rm TEMP

done


 
##### Use featurecounts to get the raw counts

featureCounts -p -B -g gene_name -a /projects/p31767/Bioinformatics/mouse/mm10/annotation_gencode_vM23/gencode.vM23.annotation.gtf -o /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/featurecounts/gene_name_counts_s2_rf.txt -s 2 -S rf /projects/b1042/Shukla_lab/Daniela/projects/11.RNA-seq_Cg1_invitro/star/*-Aligned.out.sorted.PCR.bam

