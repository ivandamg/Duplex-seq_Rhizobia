# Du novo analysis


# 1. add Barcode to header of each file
Barcode is in file R2.

We first add the barcode to the header to the foward read (R1)

                  for FILE in $(ls *_L*_R1_fastq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,3)fastp --time=0-01:00:00 --mem-per-cpu=64G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,3)_fastp.out --error=$(echo $FILE | cut -d'_' -f1,3)_fastp.error --mail-type=END,FAIL --wrap " ~/00_Software/fastp -i $FILE -I $(echo $FILE | cut -d'_' -f1,2)_R2_fastq.gz -o /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/$(echo $FILE | cut -d'_' -f1,2,3)_UMI_fastq.gz -O /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/$(echo $FILE | cut -d'_' -f1,2)_tmp.gz --umi --umi_loc=read2 --umi_len=11 -Q -A -L -w 1 -u 100 -n 11 -Y 100 -G "; sleep  1; done

We then add the barcode to the header of the reverse read (R3)

           for FILE in $(ls *_L*_R3_fastq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,3)fastp --time=0-01:00:00 --mem-per-cpu=64G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,3)_fastp.out --error=$(echo $FILE | cut -d'_' -f1,3)_fastp.error --mail-type=END,FAIL --wrap " ~/00_Software/fastp -i $FILE -I $(echo $FILE | cut -d'_' -f1,2)_R2_fastq.gz -o /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/$(echo $FILE | cut -d'_' -f1,2,3)_UMI_fastq.gz -O /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/$(echo $FILE | cut -d'_' -f1,2)_tmp2.gz --umi --umi_loc=read2 --umi_len=11 -Q -A -L -w 1 -u 100 -n 11 -Y 100 -G "; sleep  1; done


We can now delete the tmp files.

            rm *tmp*


# 2. Trimming and quality filter

            for FILE in $(ls *_R1_UMI_fastq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,3,4)_bbduk --time=01:00:00 --mem-per-cpu=64G --ntasks=2 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,3,4)_bbduk.out  --mail-type=END,FAIL --wrap "module load UHTS/Analysis/BBMap/38.91; cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp ;  bbduk.sh in=$FILE in2=$(echo $FILE | cut -d'_' -f1,2)_R3_UMI_fastq.gz out=/data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/01_trimmedfiles/$(echo $FILE | cut -d'_' -f1,2,3,4)_clean.1.fq.gz out2=/data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/01_trimmedfiles/$(echo $FILE | cut -d'_' -f1,2,3,4)_clean.2.fq.gz ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=24"   ; sleep 1; done

     
# 3. Mapping to Rhizobia

Index reference 

              sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,2)Pi --time=0-10:30:00 --mem-per-cpu=12G --ntasks=2 --cpus-per-task=1 --output=index.out --error=index.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/01_trimmedfiles/ ; /home/imateusgonzalez/00_Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index  p_ctg_oric.fasta "




Map to reference genome

      for FILE in $(ls Argon*_L*_R1_UMI_clean.1.fq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,2)_bw --time=0-05:00:00 --mem-per-cpu=64G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_bwa-mem2.out --error=$(echo $FILE | cut -d'_' -f1,2)_bwa-mem2.error --mail-type=END,FAIL --wrap "module load UHTS/Analysis/samtools/1.10; cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/01_trimmedfiles; /home/imateusgonzalez/00_Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 8 p_ctg_oric.fasta $FILE $(echo $FILE | cut -d'_' -f1,2,3,4)_clean.2.fq.gz > $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Rhizobia.sam; samtools view -bS $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Rhizobia.sam > $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Rhizobia.bam; mv $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Rhizobia.bam /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/02_MapRhizobiaDirect/ ;  rm $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Rhizobia.sam "; sleep  1; done

Sorting


Indexing



# 4. Deduplication 

Use UMI-tools to deduplicate the bam file. https://umi-tools.readthedocs.io/en/latest/QUICK_START.html#step-3--extract-the-umis

umi_tools dedup -I example.bam --paired --output-stats=deduplicated -S deduplicated.bam
