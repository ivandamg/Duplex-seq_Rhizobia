# Du novo analysis


# 1. add Barcode to header of each file




# 2. Add Illumina adaptor to begining of each file


      cat R2.out.fq | head | sed '2~4s/^/@@@@@@@/'


# 3. Add UMI barcode to each fastq read


      awk '(FNR) % 4 == 1 { -F; seq=$8; next }
     (FNR) % 4 == 2 { line[FNR]=$0; print $0 seq}' R1test.fq
