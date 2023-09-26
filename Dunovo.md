# Du novo analysis


# 1. add Barcode to header of each file
Barcode is in file R2.

We first add the barcode to the header to the foward read (R1)

                  for FILE in $(ls *_L1_R1_fastq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,3)fastp --time=0-01:00:00 --mem-per-cpu=64G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,3)_fastp.out --error=$(echo $FILE | cut -d'_' -f1,3)_fastp.error --mail-type=END,FAIL --wrap " ~/00_Software/fastp -i $FILE -I $(echo $FILE | cut -d'_' -f1,2)_R2_fastq.gz -o /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/$(echo $FILE | cut -d'_' -f1,2,3)_UMI_fastq.gz -O /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/$(echo $FILE | cut -d'_' -f1,2)_tmp.gz --umi --umi_loc=read2 --umi_len=11 -Q -A -L -w 1 -u 100 -n 11 -Y 100 -G "; sleep  1; done

We then add the barcode to the header of the reverse read (R3)

 for FILE in $(ls *_L1_R3_fastq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,3)fastp --time=0-01:00:00 --mem-per-cpu=64G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,3)_fastp.out --error=$(echo $FILE | cut -d'_' -f1,3)_fastp.error --mail-type=END,FAIL --wrap " ~/00_Software/fastp -i $FILE -I $(echo $FILE | cut -d'_' -f1,2)_R2_fastq.gz -o /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/$(echo $FILE | cut -d'_' -f1,2,3)_UMI_fastq.gz -O /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/$(echo $FILE | cut -d'_' -f1,2)_tmp2.gz --umi --umi_loc=read2 --umi_len=11 -Q -A -L -w 1 -u 100 -n 11 -Y 100 -G "; sleep  1; done


We can now delete the tmp files.

            rm *tmp*


# 2. demultiplex with fastp

     
