# Duplex-seq_Rhizobia
Duplex-seq analysis of rhizobial sequencing


# Sequencing 

The sequencing output consist of 6 files per samples:

L1.... R1
L1.... R2
L1 .... R3
L2.....R1
L2.....R2
L2.....R3

The sequencing was paired-end with UMI barcodes and on two different lanes per sample.

L1: Lane number
R1: Reads Forward
R3: Reads Reverse
R2: UMI Barcodes.

fastp https://github.com/OpenGene/fastp/issues/23
UMI_tools https://umi-tools.readthedocs.io


# 1. change names to standard.
cp WT3_NEB_L2_R2_001_6FQSNPjRSwkh.fastq.gz WT3_L2_R2_fastq.gz


# 2. add Barcode to header of each file
Barcode is in file R2.

We first add the barcode to the header to the foward read (R1)

                  for FILE in $(ls *_L*_R1_fastq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,3)fastp --time=0-01:00:00 --mem-per-cpu=64G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,3)_fastp.out --error=$(echo $FILE | cut -d'_' -f1,3)_fastp.error --mail-type=END,FAIL --wrap " cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/01_beforeTrim; ~/00_Software/fastp -i $FILE -I $(echo $FILE | cut -d'_' -f1,2)_R2_fastq.gz -o /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/02_WithUMI/$(echo $FILE | cut -d'_' -f1,2,3)_UMI_fastq.gz -O /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/02_WithUMI/$(echo $FILE | cut -d'_' -f1,2)_tmp.gz --umi --umi_loc=read2 --umi_len=11 -Q -A -L -w 1"; sleep  1; done

We then add the barcode to the header of the reverse read (R3)

           for FILE in $(ls *_L*_R3_fastq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,3)fastp --time=0-01:00:00 --mem-per-cpu=64G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,3)_fastp.out --error=$(echo $FILE | cut -d'_' -f1,3)_fastp.error --mail-type=END,FAIL --wrap " cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/01_beforeTrim; ~/00_Software/fastp -i $FILE -I $(echo $FILE | cut -d'_' -f1,2)_R2_fastq.gz -o /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/02_WithUMI/$(echo $FILE | cut -d'_' -f1,2,3)_UMI_fastq.gz -O /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/02_WithUMI/$(echo $FILE | cut -d'_' -f1,2)_tmp2.gz --umi --umi_loc=read2 --umi_len=11 -Q -A -L -w 1"; sleep  1; done


We can now delete the tmp files.

            rm *tmp*


Modify read name 
|sed -E 's/^(@[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+):/\1_/g'|

          for FILE in $(ls *_L*_R*_fastq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,3)fastp --time=0-01:00:00 --mem-per-cpu=64G --ntasks=2 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,3)_fastp.out --error=$(echo $FILE | cut -d'_' -f1,3)_fastp.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/02_WithUMI/; zcat $FILE |sed -E 's/^(@[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+:[^:]+):/\1_/g'| gzip -f > $(echo $FILE | cut -d'_' -f1,2,3,4)_NewName_fastq.gz "; sleep  1; done



# 3. Mapping to Rhizobia

Index reference 

              sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,2)Pi --time=0-10:30:00 --mem-per-cpu=12G --ntasks=2 --cpus-per-task=1 --output=index.out --error=index.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/01_trimmedfiles/ ; /home/imateusgonzalez/00_Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2 index  p_ctg_oric.fasta "




Map to reference genome

          for FILE in $(ls Argon1_L*_R1_UMI_NewName_fastq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,2)_bw --time=0-05:00:00 --mem-per-cpu=64G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_bwa-mem2.out --error=$(echo $FILE | cut -d'_' -f1,2)_bwa-mem2.error --mail-type=END,FAIL --wrap "module load UHTS/Analysis/samtools/1.10; cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/04*/; /home/imateusgonzalez/00_Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2 mem -t 8 p_ctg_oric.fasta $FILE $(echo $FILE | cut -d'_' -f1,2)_R3_UMI_NewName_fastq.gz > $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Rhizobia.sam; samtools view -bS $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Rhizobia.sam > $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Rhizobia.bam; mv $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Rhizobia.bam /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/02_MapRhizobiaDirect/ ;  rm $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Rhizobia.sam "; sleep  1; done

Sorting and indexing

            for FILE in $(ls Argon1_L1_bwa-mem2_Rhizobia.bam); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,2)ST --time=0-03:00:00 --mem-per-cpu=64G --ntasks=2 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_ST.out --error=$(echo $FILE | cut -d'_' -f1,2)_ST.error --mail-type=END,FAIL --wrap "module load UHTS/Analysis/samtools/1.10; cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/04*/; samtools sort $FILE -o $(echo $FILE | cut -d'.' -f1)_Sorted.bam; samtools index $(echo $FILE | cut -d'.' -f1)_Sorted.bam; "; sleep 1; done



# 5. Deduplication 

Use UMI-tools to deduplicate the bam file. https://umi-tools.readthedocs.io/en/latest/QUICK_START.html#step-3--extract-the-umis

        for FILE in $(ls *Sorted.bam); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,2)ST --time=0-03:00:00 --mem-per-cpu=64G --ntasks=2 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_ST.out --error=$(echo $FILE | cut -d'_' -f1,2)_ST.error --mail-type=END,FAIL --wrap "conda activate rhizo; module load UHTS/Analysis/samtools/1.10; cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/; samtools index $FILE; umi_tools dedup -I $FILE --paired --output-stats=$(echo $FILE | cut -d'_' -f1,2)_deduplicated -S $(echo $FILE | cut -d'_' -f1,2)_deduplicated.bam "; sleep 1; done


# 6. Group

ex. https://www.biostars.org/p/279951/


            for FILE in $(ls Argon1_L1*Sorted.bam); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,2)ST --time=0-03:00:00 --mem-per-cpu=64G --ntasks=2 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_ST.out --error=$(echo $FILE | cut -d'_' -f1,2)_ST.error --mail-type=END,FAIL --wrap "conda activate rhizo; module load UHTS/Analysis/samtools/1.10; cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/04*/; samtools index $FILE; umi_tools group -I $FILE --paired --group-out=$(echo $FILE | cut -d'_' -f1,2)_grouped.tsv --output-bam -S $(echo $FILE | cut -d'_' -f1,2)_deduplicated_TAGGED.bam --read-length -L $(echo $FILE | cut -d'_' -f1,2)_stats.txt"; sleep 1; done

# 6. Variant calling

        for FILE in $(ls *deduplicated.bam); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,2)ST2 --time=0-03:00:00 --mem-per-cpu=64G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_ST.out --error=$(echo $FILE | cut -d'_' -f1,2)_FB.error --mail-type=END,FAIL --wrap "cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/10_fastp/; module add UHTS/Analysis/samtools/1.10; module load UHTS/Analysis/EPACTS/3.2.6; bcftools mpileup --threads 8 -a AD,DP,SP -f p_ctg_oric.fasta $FILE | bcftools call --threads 8 -mv -Ov -o $(echo $FILE | cut -d'_' -f1,2).vcf; bcftools view --threads 8 --exclude 'QUAL <= 30 ' $(echo $FILE | cut -d'_' -f1,2).vcf -Oz -o $(echo $FILE | cut -d'_' -f1,2)_bcftoolsV1_Q30.vcf.gz"   ; sleep 1; done
























################### -obsolete



1. Clean reads with bbmap. ( only R1 and R3) forward and reverse.

      for FILE in $(ls *_R1_fastq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,3,4)_bbduk --time=01:00:00 --mem-per-cpu=64G --ntasks=2 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,3,4)_bbduk.out  --mail-type=END,FAIL --wrap "module load UHTS/Analysis/BBMap/38.91; cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/01_beforeTrim;  bbduk.sh in=$FILE in2=$(echo $FILE | cut -d'_' -f1,2)_R3_fastq.gz out=$(echo $FILE | cut -d'_' -f1,2)_clean.1.fq.gz out2=$(echo $FILE | cut -d'_' -f1,2)_clean.2.fq.gz ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=24"   ; sleep 1; done

2. Make genome ref indexes



3. map reads to medicago

    for FILE in $(ls Free*_clean.1.fq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,2)_bw --time=0-10:00:00 --mem-per-cpu=64G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_bwa-mem2.out --error=$(echo $FILE | cut -d'_' -f1,2)_bwa-mem2.error --mail-type=END,FAIL --wrap "module load UHTS/Analysis/samtools/1.10; cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/01_beforeTrim; /home/imateusgonzalez/00_Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.sse41 mem -t 8 Medicago_truncatula.fa $FILE $(echo $FILE | cut -d'_' -f1,2)_clean.2.fq.gz > $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Mtruncatula.sam; samtools view -bS $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Mtruncatula.sam > $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Mtruncatula.bam; mv $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Mtruncatula.bam /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/02_MapToMedicago/ ;  rm $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Mtruncatula.sam "; sleep  1; done



1.  Pre-processing UMI outputs and adding into read headers

          ~/00_Software/fastp -i /data/projects/p495_SinorhizobiumMeliloti/TEST/Argon1_NEB_L1_R1_001_R2S3cplEhFN6.fastq.gz -I /data/projects/p495_SinorhizobiumMeliloti/TEST/Argon1_NEB_L1_R2_001_xijpGOJF0jXL.fastq.gz -o /data/projects/p495_SinorhizobiumMeliloti/TEST/R1.out.fq -O /data/projects/p495_SinorhizobiumMeliloti/TEST/UMI.out.fq --umi --umi_loc=read2 --umi_len=11 -Q -A -L -w 1 -u 100 -n 11 -Y 100 -G

          ~/00_Software/fastp -i /data/projects/p495_SinorhizobiumMeliloti/TEST/Argon1_NEB_L1_R3_001_QN2XT8O0reaH.fastq.gz -I /data/projects/p495_SinorhizobiumMeliloti/TEST/Argon1_NEB_L1_R2_001_xijpGOJF0jXL.fastq.gz -o /data/projects/p495_SinorhizobiumMeliloti/TEST/R2.out.fq -O /data/projects/p495_SinorhizobiumMeliloti/TEST/R3.out.fq --umi --umi_loc=read2 --umi_len=11 -Q -A -L -w 1 -u 100 -n 11 -Y 100 -G


Protocol by [UMIERRORCONRRECT](https://github.com/stahlberggroup/umierrorcorrect/)

2. Mapping

        module load UHTS/Aligner/bwa-mem2/2.2.0
        module load UHTS/Aligner/bwa/0.7.17

        run_mapping.py -o /data/projects/p495_SinorhizobiumMeliloti/TEST/01_Mapping -r1 /data/projects/p495_SinorhizobiumMeliloti/TEST/R1.out.fq -r2 /data/projects/p495_SinorhizobiumMeliloti/TEST/R2.out.fq -r /data/projects/p495_SinorhizobiumMeliloti/01_REF/Sinorhizobium_meliloti_2011_ASM34606v1_genomic.fna -s TESTE


3. UMI clustering
    
 
                    #in cluster      
           sbatch --partition=pall --job-name=UMIERCOR --time=24:00:00 --mem-per-cpu=64G --cpus-per-task=18 --output=UMIERCOR.out --error=UMIERCOR.error --mail-type=END,FAIL --wrap "umi_error_correct.py -o /data/projects/p495_SinorhizobiumMeliloti/TEST/02_UMIClustering -b /data/projects/p495_SinorhizobiumMeliloti/TEST/01_Mapping/TESTE.sorted.bam -r /data/projects/p495_SinorhizobiumMeliloti/01_REF/Sinorhizobium_meliloti_2011_ASM34606v1_genomic.fna -s TESTE2"  
            
            module load UHTS/Analysis/salmon/0.11.2; salmon quant -i /data/projects/p782_RNA_seq_Argania_spinosa/02_SalmonQuantification/argane_index -l A -1 $FILE -2 $(echo $FILE | cut -d'_' -f1)_2_clean.fastq.gz -p 8 --validateMappings -o /data/projects/p782_RNA_seq_Argania_spinosa/02_SalmonQuantification/$(echo $FILE | cut -d'_' -f1)_quant 
        
        umi_error_correct.py -o /data/projects/p495_SinorhizobiumMeliloti/TEST/02_UMIClustering -b /data/projects/p495_SinorhizobiumMeliloti/TEST/01_Mapping/TESTE.sorted.bam -r /data/projects/p495_SinorhizobiumMeliloti/01_REF/Sinorhizobium_meliloti_2011_ASM34606v1_genomic.fna -s TESTE2


4. Create concensus reads

        get_consensus_statistics.py -o /data/projects/p495_SinorhizobiumMeliloti/TEST/03_Concensus -c /data/projects/p495_SinorhizobiumMeliloti/TEST/02_UMIClustering -hist /data/projects/p495_SinorhizobiumMeliloti/TEST/03_Concensus --output_raw s TESTE3


5. Create concensus output

6. PErform variant calling
