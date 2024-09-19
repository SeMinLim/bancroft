# DNA Compressor
## Reference Book
### Preprocessing for Assembled References
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
### Make a 32-Mers Reference Book
#### HG19
  1. Run KMC
     ```
     kmc -k32 -m64 -fa -cs1000000000 -b hg19.fasta 32mers kmc_tmp
     ```
     ```
     kmc_tools transform 32mers dump 32mers.txt
     ```
  2. Sort the file in descending order
     ```
     sort -k 2 -n -r 32mers.txt > hg19refbook.txt
     ```
#### HG38
  1. Run KMC
     ```
     kmc -k32 -m64 -fa -cs1000000000 -b hg38.fasta 32mers kmc_tmp
     ```
     ```
     kmc_tools transform 32mers dump 32mers.txt
     ```
  2. Sort the file in descending order
      ```
      sort -k 2 -n -r 32mers.txt > hg38refbook.txt
      ```
### Make a 256-Mers Reference Book
#### HG19
#### HG38
#### HG19 + HG38
  1. Run KMC (To get the 256-Mers that occurred 3 times or more)
     ```
     kmc -k256 -m64 -fa -ci3 -cs1000000000 -b @files.lst 256mers_1 kmc_tmp
     ```
     ```
     kmc_tools transfrom 256mers_1 dump 256mers_1.txt
     ```
  2. Sort the file in decending order
     ```
     sort -k 2 -n -r 256mers_1.txt > hg19hg38RefBook256Mers_1.txt
     ```
  3. Run KMC (To get the 256-Mers that occurred 2 times) 
     ```
     kmc -k256 -m64 -fa -cs1000000000 -cx2 -b @files.lst 256mers_2 kmc_tmp
     ```
     ```
     kmc_tools transform 256mers_2 dump 256mers_2.txt
     ```
  4. 2-Bit Encoding
     ```
     ./2_RefBookMaker/2BitEncoding/hg19hg38/obj/main
     ```   
## Compression
