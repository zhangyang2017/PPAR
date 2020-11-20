rsync -avhW --no-compress --progress /nfs2/hts/illumina/190910_J00107_0215_AHF3JTBBXY_1475/L78/Unde* /nfs0/BB/Gombart_Lab/1_rawDATA/RNAseq/conv/IM
rsync -avhW --no-compress --progress /nfs2/hts/illumina/190910_J00107_0215_AHF3JTBBXY_1475/L78/*GF* /nfs0/BB/Gombart_Lab/1_rawDATA/RNAseq/
rsync -avhW --no-compress --progress /nfs0/BB/Gombart_Lab/1_rawDATA/RNAseq/conv2/EWAT/*.gz /nfs0/BB/Gombart_Lab/1_rawDATA/RNAseq/conv2/EWAT/raw/

#cat /nfs0/BB/Gombart_Lab/1_rawDATA/RNAseq/conv/ewat/combined/gfewat10_raw.fastq /nfs0/BB/Gombart_Lab/1_rawDATA/RNAseq/conv2/EWAT/combined/gfewat10_raw.fastq > ewat70_cat.fastq


cat lane7*LFD_12_* lane8*LFD_12_* > combined/ewat12_raw.fastq
cat lane7*GF*_10_* lane8*GF*_10_* > combined/gfewat10_raw.fastq


bbduk.sh -Xmx1g in=../combined/IM01_cat.fastq out=IM01_clean2.fastq ref=/nfs0/BB/Gombart_Lab/5_RNA_stuff/2019txn/polyA.fa,/nfs0/BB/Gombart_Lab/5_RNA_stuff/2019txn/truseq_rna.fa k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=20 minlength=20


bbduk.sh in1=/nfs0/BB/Gombart_Lab/5_RNA_stuff/2019txn/LV/all/cleaned/LV56_clean.fastq out1=clean.fastq ref=/nfs0/BB/Gombart_Lab/5_RNA_stuff/2019txn/ribo_mm_hs.fa k=31 mm=f


#https://www.ncbi.nlm.nih.gov/nuccore/BK000964.3?report=fasta
#https://www.ncbi.nlm.nih.gov/nuccore/U13369.1?report=fasta

## Version 37.95

STAR --runThreadN 10 --genomeDir /nfs0/BB/Gombart_Lab/5_RNA_stuff/RNAseq/mmSTAR --sjdbGTFfile /nfs0/BB/Gombart_Lab/5_RNA_stuff/RNAseq/GRCm38/gencode.vM22.primary_assembly.annotation.gtf --sjdbOverhang 80 --readFilesIn clean.fastq --outFileNamePrefix mapped_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard


