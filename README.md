# DNA Compressor
## Reference Book
+ GRCH37 / HG19
  1. Download GRCH37
     
     https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
  2. Delete 'N'
     ```
     sed "s/N//g" Homo_sapiens.GRCh37.dna.primary_assembly.fa >> grch37_tmp.fasta
     ```
  3. Delete white spaces
     ```
     sed '/^$/d' grch37_tmp.fasta >> grch37.fasta
     ```
  4. Run our code
     ```
     ./obj/main
     ```
  5. Run KMC
     ```
     kmc -k32 -m64 -fa -cs1000000000 -b hg19.fasta 32mers kmc_tmp
     ```
     ```
     kmc_tools transform 32mers dump 32mers.txt
     ```
  6. Sort the file in descending order
     ```
     sort -k 2 -n -r 32mers.txt > hg19refbook.txt
     ```
+ GRCH38 / HG38
  1. Download GRCH38
     
     https://ftp.ensembl.org/pub/release-112/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
  2. Delete 'N'
     ```
     sed "s/N//g" Homo_sapiens.GRCh38.dna.primary_assembly.fa >> grch38_tmp.fasta
     ```
  3. Delete white spaces
     ```
     sed '/^$/d' grch38_tmp.fasta >> grch38.fasta
     ```
  4. Run our code
     ```
     ./obj/main
     ```
  5. Run KMC
     ```
     kmc -k32 -m64 -fa -cs1000000000 -b hg38.fasta 32mers kmc_tmp
     ```
     ```
     kmc_tools transform 32mers dump 32mers.txt
     ```
  6. Sort the file in descending order
      ```
      sort -k 2 -n -r 32mers.txt > hg38refbook.txt
      ```
## Compression
