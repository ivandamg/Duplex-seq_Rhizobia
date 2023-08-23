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

# fastp analysis
https://github.com/OpenGene/fastp/issues/23

# 1. change names to standard.
cp WT3_NEB_L2_R2_001_6FQSNPjRSwkh.fastq.gz WT3_L2_R2_fastq.gz



# 1. Clean reads with bbmap. ( only R1 and R3) forward and reverse.

      for FILE in $(ls *_R1_fastq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,3,4)_bbduk --time=01:00:00 --mem-per-cpu=64G --ntasks=2 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,3,4)_bbduk.out  --mail-type=END,FAIL --wrap "module load UHTS/Analysis/BBMap/38.91; cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/01_beforeTrim;  bbduk.sh in=$FILE in2=$(echo $FILE | cut -d'_' -f1,2)_R3_fastq.gz out=$(echo $FILE | cut -d'_' -f1,2)_clean.1.fq.gz out2=$(echo $FILE | cut -d'_' -f1,2)_clean.2.fq.gz ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=24"   ; sleep 1; done

# 2. Make genome ref indexes



#3. map reads to medicago

    for FILE in $(ls Free*_clean.1.fq.gz); do echo $FILE; sbatch --partition=pall --job-name=$(echo $FILE | cut -d'_' -f1,2)_bw --time=0-10:00:00 --mem-per-cpu=64G --ntasks=8 --cpus-per-task=1 --output=$(echo $FILE | cut -d'_' -f1,2)_bwa-mem2.out --error=$(echo $FILE | cut -d'_' -f1,2)_bwa-mem2.error --mail-type=END,FAIL --wrap "module load UHTS/Analysis/samtools/1.10; cd /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/01_beforeTrim; /home/imateusgonzalez/00_Software/bwa-mem2-2.2.1_x64-linux/bwa-mem2.sse41 mem -t 8 Medicago_truncatula.fa $FILE $(echo $FILE | cut -d'_' -f1,2)_clean.2.fq.gz > $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Mtruncatula.sam; samtools view -bS $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Mtruncatula.sam > $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Mtruncatula.bam; mv $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Mtruncatula.bam /data/projects/p495_SinorhizobiumMeliloti/02_DuplexSeq/02_MapToMedicago/ ;  rm $(echo $FILE | cut -d'_' -f1,2)_bwa-mem2_Mtruncatula.sam "; sleep  1; done




# 1.  Pre-processing UMI outputs and adding into read headers

          ~/00_Software/fastp -i /data/projects/p495_SinorhizobiumMeliloti/TEST/Argon1_NEB_L1_R1_001_R2S3cplEhFN6.fastq.gz -I /data/projects/p495_SinorhizobiumMeliloti/TEST/Argon1_NEB_L1_R2_001_xijpGOJF0jXL.fastq.gz -o /data/projects/p495_SinorhizobiumMeliloti/TEST/R1.out.fq -O /data/projects/p495_SinorhizobiumMeliloti/TEST/UMI.out.fq --umi --umi_loc=read2 --umi_len=11 -Q -A -L -w 1 -u 100 -n 11 -Y 100 -G

          ~/00_Software/fastp -i /data/projects/p495_SinorhizobiumMeliloti/TEST/Argon1_NEB_L1_R3_001_QN2XT8O0reaH.fastq.gz -I /data/projects/p495_SinorhizobiumMeliloti/TEST/Argon1_NEB_L1_R2_001_xijpGOJF0jXL.fastq.gz -o /data/projects/p495_SinorhizobiumMeliloti/TEST/R2.out.fq -O /data/projects/p495_SinorhizobiumMeliloti/TEST/R3.out.fq --umi --umi_loc=read2 --umi_len=11 -Q -A -L -w 1 -u 100 -n 11 -Y 100 -G


Protocol by [UMIERRORCONRRECT](https://github.com/stahlberggroup/umierrorcorrect/)

#2. Mapping

        module load UHTS/Aligner/bwa-mem2/2.2.0
        module load UHTS/Aligner/bwa/0.7.17

        run_mapping.py -o /data/projects/p495_SinorhizobiumMeliloti/TEST/01_Mapping -r1 /data/projects/p495_SinorhizobiumMeliloti/TEST/R1.out.fq -r2 /data/projects/p495_SinorhizobiumMeliloti/TEST/R2.out.fq -r /data/projects/p495_SinorhizobiumMeliloti/01_REF/Sinorhizobium_meliloti_2011_ASM34606v1_genomic.fna -s TESTE


#3. UMI clustering
    
 
                    #in cluster      
           sbatch --partition=pall --job-name=UMIERCOR --time=24:00:00 --mem-per-cpu=64G --cpus-per-task=18 --output=UMIERCOR.out --error=UMIERCOR.error --mail-type=END,FAIL --wrap "umi_error_correct.py -o /data/projects/p495_SinorhizobiumMeliloti/TEST/02_UMIClustering -b /data/projects/p495_SinorhizobiumMeliloti/TEST/01_Mapping/TESTE.sorted.bam -r /data/projects/p495_SinorhizobiumMeliloti/01_REF/Sinorhizobium_meliloti_2011_ASM34606v1_genomic.fna -s TESTE2"  
            
            module load UHTS/Analysis/salmon/0.11.2; salmon quant -i /data/projects/p782_RNA_seq_Argania_spinosa/02_SalmonQuantification/argane_index -l A -1 $FILE -2 $(echo $FILE | cut -d'_' -f1)_2_clean.fastq.gz -p 8 --validateMappings -o /data/projects/p782_RNA_seq_Argania_spinosa/02_SalmonQuantification/$(echo $FILE | cut -d'_' -f1)_quant 
        
        umi_error_correct.py -o /data/projects/p495_SinorhizobiumMeliloti/TEST/02_UMIClustering -b /data/projects/p495_SinorhizobiumMeliloti/TEST/01_Mapping/TESTE.sorted.bam -r /data/projects/p495_SinorhizobiumMeliloti/01_REF/Sinorhizobium_meliloti_2011_ASM34606v1_genomic.fna -s TESTE2


#4. Create concensus reads

        get_consensus_statistics.py -o /data/projects/p495_SinorhizobiumMeliloti/TEST/03_Concensus -c /data/projects/p495_SinorhizobiumMeliloti/TEST/02_UMIClustering -hist /data/projects/p495_SinorhizobiumMeliloti/TEST/03_Concensus --output_raw s TESTE3


#5. Create concensus output

#6. PErform variant calling
