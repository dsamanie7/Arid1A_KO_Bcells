
### Map

for i in $(cat /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/index/map.txt)
do
echo ${i} > TEMP
OLD=$(cut -f1 -d: TEMP)

cat <<EOT >> /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/scripts/bowtie_mm10_map_${OLD}_PBS.sh
#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=12:00:00
#SBATCH --job-name=bowtie_only
#SBATCH --output=/projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/scripts/bowtie_mm10_map_${OLD}_PBS.out
#SBATCH --error=/projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/scripts/bowtie_mm10_map_${OLD}_PBS.err

cd /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay

module load bowtie
module load samtools/1.10.1
module load java/jdk1.8.0_191
module load picard/2.21.4
module load deeptools



bowtie -p 12 -m 1 --best --strata -X 2000 -S --fr --chunkmbs 1024 /projects/p31767/Bioinformatics/mouse/mm10/mm10_bowtie_index/mm10 -1 /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/fastq/${OLD}_R1_001.fastq -2 /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/fastq/${OLD}_R2_001.fastq -S  /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.sam 2> /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.err

samtools view -bS /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.sam > /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.bam

samtools sort /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.bam -o /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.sorted.bam 

samtools index  /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.sorted.bam 

bamCoverage -b  /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.sorted.bam  --normalizeUsing CPM -o  /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.sorted.bw

picard MarkDuplicates I=/projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.sorted.bam O=/projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.sorted_PCR.bam M=metrics.txt REMOVE_DUPLICATES=true

samtools index /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.sorted_PCR.bam

bamCoverage --samFlagInclude 64 -b /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.sorted_PCR.bam --normalizeUsing CPM -o /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.sorted_PCR_paired.bw


EOT

sbatch /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/scripts/bowtie_mm10_map_${OLD}_PBS.sh

rm TEMP
done




########### call peaks

for i in $(cat /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/index/map.txt)
do
echo ${i} > TEMP
OLD=$(cut -f1 -d: TEMP)

cat <<EOT >> /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/scripts/tagdir_${OLD}_PBS.sh
#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=4:00:00
#SBATCH --job-name=bowtie
#SBATCH --output=/projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/scripts/tagdir_${OLD}_PBS.out
#SBATCH --error=/projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/scripts/tagdir_${OLD}_PBS.err

cd /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay

module load homer
module load samtools/1.10.1
module load java/jdk1.8.0_191
module load picard/2.21.4
module load deeptools


makeTagDirectory /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/tag_directories/${OLD} -genome mm10 /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/bowtie/${OLD}_accepted_hits_mm10.sorted_PCR.bam 

EOT

sbatch /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/scripts/tagdir_${OLD}_PBS.sh

rm TEMP
done



for i in $(cat /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/index/map.txt)
do
echo ${i} > TEMP
OLD=$(cut -f1 -d: TEMP)

cat <<EOT >> /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/scripts/peaks_${OLD}_PBS.sh
#!/bin/bash

#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=4:00:00
#SBATCH --job-name=bowtie
#SBATCH --output=/projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/scripts/peaks_${OLD}_PBS.out
#SBATCH --error=/projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/scripts/peaks_${OLD}_PBS.err

cd /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay

module load homer
module load samtools/1.10.1
module load java/jdk1.8.0_191
module load picard/2.21.4
module load deeptools


findPeaks /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/tag_directories/${OLD} -center -style dnase -o /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/peaks/${OLD} 

EOT

sbatch /projects/b1042/Shukla_lab/Daniela/projects/19.ATACseq_Ajay/scripts/peaks_${OLD}_PBS.sh

rm TEMP
done


#### merge the replicates

mergePeaks *cre* -venn venn_mergepeaks_cre.txt >  mergepeaks_cre.txt

mergePeaks *KO* -venn venn_mergepeaks_KO.txt >  mergepeaks_KO.txt


##### select the intersect  from replicates

grep 'R_4218_cre-_S1|R_4239_cre-_S2' mergepeaks_cre.txt > intersection_mergepeaks_control.txt
grep 'R_4216_KO_S3|R_4240_KO_S4' mergepeaks_KO.txt > intersection_mergepeaks_KO.txt


### get the union of each intersected peaks set

mergePeaks intersection_mergepeaks_* -venn venn_union_mergepeaks_controlKO.txt > union_mergepeaks_controlKO.txt

