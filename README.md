# :dna: Microbiome analysis and transcriptional profile of sugarcane varieties in resistance to the fungus _Sporisorium scitamineum_

This repository will be used to describe, in detail, the data processing and results related to the _Sporisorium scitamineum_ microbiome project. The experiment consists of two varieties of sugar cane in low soils
and high fertility. The seedlings were planted in December 2022 at Usina Iracema, in São Martinho (SP), and the plants were  kept in the field for seven months. Experiments were arranged in random blocks with six replicates.
Each experimental plot consists of four rows measuring four meters in length and planting spacing of 0.5 m, totaling 32 plants per plot. As varieties chosen are IACSP01-5503 (adapted to low fertility soils and
resistant to _S. scitamineum_) and IACSP-6007 (not adapted to low fertility soils and with intermediate-resistance to _S. scitamineum_).

## :package: Raw data storage: Where is our data located?

1. :open_file_folder: Raw files from the table are located in ```fisher.esalq.usp.br```, specifically at:
   
```/media/hdzin/pedro/NGS726/data/raw_data```

> It is important to note that the external HDD that we are using is usually buggy, so if you try to access the folder with ```cd``` and it doesn't work, it probably means the HDD is not mounted. So use the following codes: ```sudo umount /media/hdzin/``` and ```sudo mount -a```. In case you don't possess sudo permissions, please contact Renato Bombardelli.

In the below table, we describe the file correspondences to their respective treatment and condition:

| File | Biological replicate | Plant Genotype | Condition | Treatment |
|------------- | ------------- | ------------- | ------------- | ------------- |
| NGS726_1_S1_L001_R1_001.fastq  | 1 | IACSP-5503 | Sandy soil | Inoculated |
| NGS726_1_S1_L001_R2_001.fastq  | 1 | IACSP-5503 | Sandy soil | Inoculated |
| NGS726_2_S2_L001_R1_001.fastq  | 2 | IACSP-5503 | Sandy soil | Inoculated |
| NGS726_2_S2_L001_R2_001.fastq  | 2 | IACSP-5503 | Sandy soil | Inoculated |
| NGS726_3_S3_L001_R1_001.fastq  | 3 | IACSP-5503 | Sandy soil | Inoculated |
| NGS726_3_S3_L001_R2_001.fastq  | 3 | IACSP-5503 | Sandy soil | Inoculated |
| NGS726_4_S4_L001_R1_001.fastq  | 1 | IACSP-6007 | Sandy soil | Inoculated |
| NGS726_4_S4_L001_R2_001.fastq  | 1 | IACSP-6007 | Sandy soil | Inoculated |
| NGS726_5_S5_L001_R1_001.fastq  | 2 | IACSP-6007 | Sandy soil | Inoculated |
| NGS726_5_S5_L001_R2_001.fastq  | 2 | IACSP-6007 | Sandy soil | Inoculated |
| NGS726_6_S6_L001_R1_001.fastq  | 3 | IACSP-6007 | Sandy soil | Inoculated |
| NGS726_6_S6_L001_R2_001.fastq  | 3 | IACSP-6007 | Sandy soil | Inoculated |
| NGS726_7_S7_L001_R1_001.fastq  | 1 | IACSP-5503 | Sandy soil | Control |
| NGS726_7_S7_L001_R2_001.fastq  | 1 | IACSP-5503 | Sandy soil | Control |
| NGS726_8_S8_L001_R1_001.fastq  | 2 | IACSP-5503 | Sandy soil | Control |
| NGS726_8_S8_L001_R2_001.fastq  | 2 | IACSP-5503 | Sandy soil | Control |
| NGS726_9_S9_L001_R1_001.fastq  | 3 | IACSP-5503 | Sandy soil | Control |
| NGS726_9_S9_L001_R2_001.fastq  | 3 | IACSP-5503 | Sandy soil | Control |
| NGS726_10_S10_L001_R1_001.fastq  | 1 | IACSP-6007 | Sandy soil | Control |
| NGS726_10_S10_L001_R2_001.fastq  | 1 | IACSP-6007 | Sandy soil | Control |
| NGS726_11_S11_L001_R1_001.fastq  | 2 | IACSP-6007 | Sandy soil | Control |
| NGS726_11_S11_L001_R2_001.fastq  | 2 | IACSP-6007 | Sandy soil | Control |
| NGS726_12_S12_L001_R1_001.fastq  | 3 | IACSP-6007 | Sandy soil | Control |
| NGS726_12_S12_L001_R2_001.fastq  | 3 | IACSP-6007 | Sandy soil | Control |
| NGS726_13_S13_L001_R1_001.fastq  | 1 | IACSP-5503 | Clay soil | Inoculated |
| NGS726_13_S13_L001_R2_001.fastq  | 1 | IACSP-5503 | Clay soil | Inoculated |
| NGS726_14_S14_L001_R1_001.fastq  | 2 | IACSP-5503 | Clay soil | Inoculated |
| NGS726_14_S14_L001_R2_001.fastq  | 2 | IACSP-5503 | Clay soil | Inoculated |
| NGS726_15_S15_L001_R1_001.fastq  | 3 | IACSP-5503 | Clay soil | Inoculated |
| NGS726_15_S15_L001_R2_001.fastq  | 3 | IACSP-5503 | Clay soil | Inoculated |
| NGS726_16_S16_L001_R1_001.fastq  | 1 | IACSP-6007 | Clay soil | Inoculated |
| NGS726_16_S16_L001_R2_001.fastq  | 1 | IACSP-6007 | Clay soil | Inoculated |
| NGS726_17_S17_L001_R1_001.fastq  | 2 | IACSP-6007 | Clay soil | Inoculated |
| NGS726_17_S17_L001_R2_001.fastq  | 2 | IACSP-6007 | Clay soil | Inoculated |
| NGS726_18_S18_L001_R1_001.fastq  | 3 | IACSP-6007 | Clay soil | Inoculated |
| NGS726_18_S18_L001_R2_001.fastq  | 3 | IACSP-6007 | Clay soil | Inoculated |
| NGS726_19_S19_L001_R1_001.fastq  | 1 | IACSP-5503 | Clay soil | Control |
| NGS726_19_S19_L001_R2_001.fastq  | 1 | IACSP-5503 | Clay soil | Control |
| NGS726_20_S20_L001_R1_001.fastq  | 2 | IACSP-5503 | Clay soil | Control |
| NGS726_20_S20_L001_R2_001.fastq  | 2 | IACSP-5503 | Clay soil | Control |
| NGS726_21_S21_L001_R1_001.fastq  | 3 | IACSP-5503 | Clay soil | Control |
| NGS726_21_S21_L001_R2_001.fastq  | 3 | IACSP-5503 | Clay soil | Control |
| NGS726_22_S22_L001_R1_001.fastq  | 1 | IACSP-6007 | Clay soil | Control |
| NGS726_22_S22_L001_R2_001.fastq  | 1 | IACSP-6007 | Clay soil | Control |
| NGS726_23_S23_L001_R1_001.fastq  | 2 | IACSP-6007 | Clay soil | Control |
| NGS726_23_S23_L001_R2_001.fastq  | 2 | IACSP-6007 | Clay soil | Control |
| NGS726_24_S24_L001_R1_001.fastq  | 3 | IACSP-6007 | Clay soil | Control |
| NGS726_24_S24_L001_R2_001.fastq  | 3 | IACSP-6007 | Clay soil | Control |

## :hammer_and_wrench: Data processing: What did we do?

First of all, we needed to assess the read quality from the sequencing data. How confident we are in the sequenced bases? For that, we used **FastQC v.0.11.5** (Andrews, 2010). We used this quick code below:
```
#!/bin/bash
mkdir -p ./FastQC_raw
for i in *fastq.gz
do
fastqc "$i" -o ./FastQC_raw -t 8
done
```

The output results are located in :open_file_folder:  ```/media/hdzin/pedro/NGS726/data/FastQC_raw```

Next, read filtering step consisted in removing library adapters, low quality reads (Q < 20) and ambiguous bases (N) using **Cutadapt v.1.18** (Martin, 2011):
```
#!/bin/bash

for i in ../raw_data/*_R1_001.fastq.gz
do
    file2=`echo "$i" | sed "s/_R1/_R2/g"`
    output1=$(basename "$i" | sed "s/_R1_001.fastq.gz//g")

    cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
    -m 50 --max-n 2 -q 20 -j 6 \
    -o "$output1"_cut_PE1.fastq.gz \
    -p "$output1"_cut_PE2.fastq.gz \
    "./$i"  "./$file2" > "$output1"_cutadapt_summary.txt

done
```
>Note in the code above, we specifically gave as input the adapter sequences we want the software to look for. **Always** check the correct adapters used during the sequencing library to give as input to the program. Because they might be different for each one of your data sets.

The output results are located in :open_file_folder:  ```/media/hdzin/pedro/NGS726/data/cutadapt/```

We then ran FastQC again on trimmed reads using the same code as above. To concatenate the huge amount of results, we used ```MultiQC``` software.
1. :link: Raw read quality assessment can be assessed downloading the HTML file [here](https://github.com/pedrofvilanova/smut_microbiome/blob/main/multiqc_report_raw.html).
2. :link: Trimmed read quality assessment can be assessed downloading the HTML file [here](https://github.com/pedrofvilanova/smut_microbiome/blob/main/multiqc_report_trimmed.html).

3. The output of FastQC on trimmed reads can be accessed: ```/media/hdzin/pedro/NGS726/data/cutadapt/FastQC_trimmed/```

You can also access a PDF file containing tables summarizing the total number of raw reads and the total number of reads after trimming, clicking [here](https://github.com/pedrofvilanova/smut_microbiome/blob/main/smut_microbiome_rnaseq_reads_report.pdf)



## :dna: Transcript alignment:
### FISHER SERVER -  BY PEDRO VILANOVA
The analysis executed in this step are inside the folder :open_file_folder:  ```/media/hdzin/pedro/NGS726/data/cutadapt/rRNA_clean```

>Inside this folder you will find the statistics for the rRNA cleaning step and folders referring to the two sugarcane genotypes IACSP-5503 and IACSP-6007

The analyses were done in the following order:
1. Get rRNA cleaned reads and map them first to _Sporisorium scitamineum_ reference genome located in  :open_file_folder: ```/media/hdzin/pedro/NGS726/data/cutadapt/rRNA_clean/references/NCBI_Ssci_genome.fa```

2. Get unmapped reads to _S. scitamineum_ genome (reads that are not fungus) and map them to COMPGG transcripts located in :open_file_folder: ```/media/hdzin/pedro/NGS726/data/cutadapt/rRNA_clean/references/comps_plus_GGs.fas.gz```

3. Get unclassified reads that are not mapped neither to _S. scitamineum_ nor sugarcane COMPGG transcripts.

The directory structures is as follows:

```
- rRNA_clean
  - IACSP-5503
    - clay
      - control
      - inoc
    - sandy
      - control
      - inoc
  - IACSP-6007
    - clay
      - control
      - inoc
    - sandy
      - control
      - inoc
```
   
>Inside each one of the folders ```control``` or ```inoc``` there is a folder called ```clean_reads``` which contains ```.fastq``` trimmed and rRNA clean files of the respective treatment and condition.

For transcript alignment, we used ```Hisat2 v.2.1.0```. First of all, we align trimmed and rRNA cleaned reads to Sporisorium scitamineum genome, no matter if the treatment is control or inoculated. We hope that mock-inoculated plants will not have any or a very small number (due to spurious alignments) of alignments to the fungus genome. After, aligning reads to the fungus genome, we're gonna retrieve the reads that were *unmapped*, because if they do not map to the fungus, they belong to some other organisms, likely sugarcane. 

Mapping to _S. scitamineum_ will always be stored in a folder called ```1_mapping2_ssci``` for each one of the treatments and conditions. Inside a folder like that you will find three codes: 

1. ```1_run_mapping.sh```: responsible for running the map of .fastq reads to S. scitamineum genome
2. ```2_get_unmapped.sh```: responsible for retrieving only unmapped reads, which are likely anything other than the fungus.
3. ```3_get_mapped.sh```: we use this code to actually separate the reads that belong to the fungus, which might be of interest later.

Let's break down these codes!

Above is the content of ```1_run_mapping.sh```:

```
#!/bin/bash

ref=/media/hdzin/pedro/NGS726/data/cutadapt/rRNA_clean/references/NCBI_Ssci_genome.fa
index=/media/hdzin/pedro/NGS726/data/cutadapt/rRNA_clean/references/NCBI_Ssci_genome.fa.index

# Index reference genome
#hisat2-build "$ref" "$index"

# Map reads against genome
for i in ../clean_reads/*PE1.fastq.gz
do
        file2=`echo "$i" | sed "s/PE1/PE2/g"`
        rep=$(basename "$i" | sed -E "s/NGS726_([0-9]+).*/\1/g")

        echo "Analyzing $i and $file2"


        hisat2 -p 8 --rg-id IACSP-5503_clay_control --rg SM:IACSP-5503_clay_control_rep"$rep" \
        --summary-file ./summary_IACSP-5503_clay_control_rep"$rep"_mapped2_ssci.txt \
        -x "$index" -1 "$i" -2 "$file2" \
        -S ./IACSP-5503_clay_control_rep"$rep"_mapped2_ssci.sam

        # Map convert to bam and index mapping
        samtools sort ./IACSP-5503_clay_control_rep"$rep"_mapped2_ssci.sam > ./IACSP-5503_clay_control_rep"$rep"_mapped2_ssci.bam
        samtools index ./IACSP-5503_clay_control_rep"$rep"_mapped2_ssci.bam
        rm ./IACSP-5503_clay_control_rep"$rep"_mapped2_ssci.sam

done
```

The first lines of the code indicate where is my reference sequence to be mapped against, expressed in ```ref``` and also where the index is gonna be and how it is going to be called ```index```.
Then, the first command of ```hisat2``` will actually build the index ```#hisat2-build "$ref" "$index"```
Next, we're actually gonna map the reads against the fungus genome, as you can see I'm mapping all files inside ```clean_reads```, using both ```PE1``` and ```PE2``` reads.
The final step uses ```samtools v. 1.10```to convert the result ```.sam``` files into ```.bam``` files. Then, we index the ```.bam``` files.

Now, we have ```.bam``` files containing mapped and unmapped reads to S. scitamineum genome. We're now gonna get everything that has not been mapped to the fungus with ```2_get_unmapped.sh``` code:

```
#!/bin/bash

mkdir -p ssci_free

for i in ./*.bam
do

        # get only the unmapped reads - "4" indicates the paired end unmapped reads according to samtools
        samtools view --threads 10 -b -f 4 "$i" > "$i"_unmapped.bam

        # save fastq reads in separate R1 and R2 files
        samtools fastq -@ 10 "$i"_unmapped.bam\
        -1 ./ssci_free/"$i"_ssci_free_PE1.fastq.gz \
        -2 ./ssci_free/"$i"_ssci_free_PE2.fastq.gz \
        -0 /dev/null -s /dev/null -n

        rm "$i"_unmapped.bam
done
```

Firstly, we create a folder called ```ssci_free``` where reads free of the fungus transcripts will be stored. Then we'll use ```samtools view``` command to get the unmapped reads using the ```4``` flag, which indicates that we only want reads that are unmapped to our reference. After that, we need to convert the resulting ```.bam``` files into ```.fastq``` files again, and we use ```samtools fastq``` to do the conversion. 

Now, we want to retrieve reads that actually mapped to the fungus genome, so we're gonna use ```3_get_mapped.sh```

```
#!/bin/bash

mkdir -p mapped2_ssci

for i in ./*.bam
do

        # get only the mapped reads - "2" indicates the paired end mapped reads according to samtools
        samtools view --threads 10 -b -f 2 "$i" > "$i"_mapped.bam

        # save fastq reads in separate R1 and R2 files
        samtools fastq -@ 10 "$i"_mapped.bam\
        -1 ./mapped2_ssci/"$i"_ssci_PE1.fastq.gz \
        -2 ./mapped2_ssci/"$i"_ssci_PE2.fastq.gz \
        -0 /dev/null -s /dev/null -n

        rm "$i"_mapped.bam
done
```

We create, this time, a folder called ```mapped2_ssci``` to actually store reads that mapped to the fungus, and we're going to repeat the same step as code 2, but instead, we're gonna use flag ```2``` to get the paired-end mapped reads.

Now, inside the ```1_mapping2_ssci``` folder we're left with fungal reads in ```mapped2_ssci``` and unmapped reads that contain sugarcane reads of our interest. 


## NIOO - SERVER - The analyses were done in the following order:
1. Get "Ssci_free" reads (from fisher server - as describe before) and map them to R570_poliploid_version_2024/SofficinarumxspontaneumR570_771_v2.0.fa. reference genome located in  :open_file_folder: ```joycef@nioo0003.nioo.int:/home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data/References/R570_poliploid_version_2024```

2. Get mapped and unmapped reads to R570_poliploid_version_2024 genome and salve as
mapped in :open_file_folder: ```../Data/2Alignments_R570/R570_mapped```
unmapped in :open_file_folder: ```../Data/2Alignments_R570/R570_unmapped```

### Installing packages through the mamba environment
Packages versions : 
**HiSat2: 2.2.1**; 
**SAMtools: 1.20**; 
**featureCounts** (parte do pacote subread): **2.0.6**; 

### Executing the Script with Screen ###
```
screen -S hisat2_alignment
```
Within the screen session, run the script:
```
bash _scriptX
```
Detach from the screen session: Press Ctrl + A followed by D.

Reattach to the screen session (if needed):
```
screen -r hisat2_alignment
```

### 1 Step: Reference Index

```
#!/bin/bash

if [ ! -d ../../Data/References/R570_poliploid_version_2024/hisat2_index ]
then
    mkdir -p ../../Data/References/R570_poliploid_version_2024/hisat2_index
    hisat2-build -p 10 ../../Data/References/R570_poliploid_version_2024/SofficinarumxspontaneumR570_771_v2.0.fa ../../Data/References/R570_poliploid_version_2024/hisat2_index/SofficinarumxspontaneumR570_771_v2.0
fi

```

### 2 Step: Reference Alignment

```
#!/bin/bash

#Verifica se o diretório de mapeamento já existe, se não, cria
if [ ! -d ../../Data/2Alignments_R570 ]; then
  mkdir -p ../../Data/2Alignments_R570
fi

#Cria diretórios para leituras mapeadas e não mapeadas
mapped_dir="../../Data/2Alignments_R570/R570_mapped"
unmapped_dir="../../Data/2Alignments_R570/R570_unmapped"
mkdir -p "$mapped_dir"
mkdir -p "$unmapped_dir"

#Define o diretório de trabalho
cd ../../Data/2Alignments_R570

#Caminho para o índice do HISAT2
index="../References/R570_poliploid_version_2024/hisat2_index/SofficinarumxspontaneumR570_771_v2.0"

#Nome do arquivo de controle
checkpoint_file="./processed_files.txt"

#Verifica se o arquivo de controle existe, se não existir, cria
if [ ! -e "$checkpoint_file" ]; then
  touch "$checkpoint_file"
fi

#Função para lidar com o sinal de interrupção (Ctrl+C)
function cleanup {
  echo "Script interrupted. Exiting..."
  exit 1
}

#Captura o sinal de interrupção (Ctrl+C)
trap cleanup SIGINT

#Loop para processar arquivos
for PE1 in /home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data/rRNA_clean/Ssci_free/*PE1.fastq.gz; do
  PE2="${PE1/PE1/PE2}"
  BASE_PE1=$(basename "$PE1")

#Verifica se o par já foi processado antes de executar o HISAT2
  if ! grep -q "$BASE_PE1" "$checkpoint_File"; then
    echo "Analyzing $PE1 and $PE2"

   hisat2 -p 20 --rg-id "$PE1" --rg SM:"$PE1" \
    --summary-file ./summary_"$BASE_PE1"_Ssci_free.txt \
    -x "$index" -1 "$PE1" -2 "$PE2" \
    -S ./"$BASE_PE1"_Ssci_free.sam

    #Ordena o arquivo SAM e converte para BAM
    samtools sort -@ 20 ./"$BASE_PE1"_Ssci_free.sam > ./"$BASE_PE1"_Ssci_free.bam

    #Devido ao tamanho de cromossomo maior que 512Mb, não é possível criar o bam index do tipo .bai; usamos então o parâmetro "-c"
    samtools index -@ 20 -c ./"$BASE_PE1"_Ssci_free.bam
    rm ./"$BASE_PE1"_Ssci_free.sam

   #Filtra as leituras mapeadas e as move para o diretório de leituras mapeadas
 samtools view --threads 20 -b -F 4 ./"$BASE_PE1"_Ssci_free.bam > "$mapped_dir/${BASE_PE1}_R570mapped.bam"

    #Filtra as leituras não mapeadas e as move para o diretório de leituras não mapeadas
    samtools view --threads 20 -b -f 4 ./"$BASE_PE1"_Ssci_free.bam > "$unmapped_dir/${BASE_PE1}_R570unmapped.bam"

    #Adiciona o nome do arquivo ao arquivo de controle
    echo "$BASE_PE1" >> "$checkpoint_file"
  else
    echo "File $BASE_PE1 has been processed previously. Ignoring."
  fi
done

```

Separation of Mapped and Unmapped Reads:
After the HISAT2 classification step, the script filters the mapped and unmapped reads using samtools view.
The mapped reads are moved to the directory R570_mapped, while the unmapped reads are moved to the directory R570_unmapped.

## Reads Counts:
**featureCounts** (parte do pacote subread): **2.0.6**;

### 1 Step: Convert reference annotation - gff3 to gft

```                                                                  
#!/bin/bash

#Define os caminhos dos arquivos
input_gff3="../../Data/References/R570_poliploid_version_2024/annotation/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gff3.gz"
output_gtf="../../Data/References/R570_poliploid_version_2024/annotation/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf"

#Verifica se o arquivo de entrada está compactado
if [[ "$input_gff3" == *.gz ]]; then
    echo "Descompactando $input_gff3..."
    gunzip -c "$input_gff3" > "${input_gff3%.gz}"
    input_gff3="${input_gff3%.gz}"
fi

#Converte o GFF3 para GTF
gffread "$input_gff3" -T -o "$output_gtf"

#Verifica se a conversão foi bem-sucedida
if [ $? -eq 0 ]; then
    echo "Conversion successful: $output_gtf"
else
    echo "Conversion failed."
fi

```

### 2 Step: Reads Counts

These options ensure that the script counts paired-end read pairs, includes multi-mapping reads, uses the gene_id for grouping, and processes the data efficiently with 20 threads.

```
#!/bin/bash

#Verifica se o diretório de saída existe, caso contrário, cria
if [ ! -d ../../Data/3featureCounts_R570__ShSHN ]; then
    mkdir -p ../../Data/3featureCounts_R570__ShSHN
fi

#Muda para o diretório de saída
cd ../../Data/3featureCounts_R570__ShSHN/

#Executa o featureCounts com os parâmetros fornecidos
featureCounts -s 0 -p --countReadPairs -T 20 -t CDS -g gene_id -M \
-a /home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data/References/R570_poliploid_version_2024/annotation/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf \
-o ./Quantified_all_R570_genome_ShSHN_ref.tsv \
$(ls /home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data/2Alignments_R570/R570_mapped/*.bam) \
2> ./featureCounts_R570_genome_ShSHN_log.txt


```

## EdgeR - RStudio
Until now, we have carried out all analyzes on the server. From now on we will work through RStudio for the following analyses.
The script performs differential gene expression analysis using RNA-Seq data with the edgeR library. It includes steps to load and prepare the data, normalize the counts, perform statistical tests to identify differentially expressed genes, and export the results to CSV files.

```
library(dplyr)
library(edgeR)

# Load count table
count_table <- read.table(
  "C:/Users/joyde/OneDrive/Área de Trabalho/Desktop NIOO/RNAseq/EdgeR/Quantified_all_R570_genome_ShSHN_ref.tsv",
  header = TRUE,
  row.names = "Geneid",
  check.names = FALSE
)

# Store gene length and remove annotation columns
gene_length <- count_table[, 5]
names(gene_length) <- rownames(count_table)

count_table <- count_table[, -(1:5)]

# Load correspondence file
correspondence <- read.table(
  "C:/Users/joyde/OneDrive/Área de Trabalho/Desktop NIOO/RNAseq/EdgeR/samples_sscifree.txt",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(correspondence) <- c("old_name", "new_name")

# Create name mapping
mapping <- setNames(correspondence$new_name, correspondence$old_name)

# Fix count table column names so they match the correspondence file
old_colnames <- basename(colnames(count_table))
old_colnames <- sub("_R570mapped\\.bam$", "", old_colnames)

# Rename columns using correspondence file
new_colnames <- mapping[old_colnames]

# Stop if some columns were not found in the correspondence file
if (any(is.na(new_colnames))) {
  cat("Columns not found in correspondence file:\n")
  print(old_colnames[is.na(new_colnames)])
  stop("Some count table columns could not be renamed.")
}

colnames(count_table) <- unname(new_colnames)

# Pairwise comparison
comparison <- c("IACSP-6007_clay_inoc", "IACSP-6007_clay_control")

# Select only samples from this comparison
selected_cols <- grep(
  paste0("^(", paste(comparison, collapse = "|"), ")[0-9]+$"),
  colnames(count_table),
  value = TRUE
)

print(selected_cols)

count_table <- count_table[, selected_cols]

# Check selected samples
print(ncol(count_table))
print(colnames(count_table))

# Stop if expected number of samples was not selected
if (ncol(count_table) != 6) {
  stop("Expected 6 samples: 3 inoc and 3 control. Check selected column names.")
}

# Create sample information automatically from column names
sample_info <- data.frame(
  condition = case_when(
    grepl("_control[0-9]+$", colnames(count_table)) ~ "control",
    grepl("_inoc[0-9]+$", colnames(count_table)) ~ "inoc",
    TRUE ~ NA_character_
  )
)

rownames(sample_info) <- colnames(count_table)

# Convert condition to factor, with control as reference
sample_info$condition <- factor(
  sample_info$condition,
  levels = c("control", "inoc")
)

# Check final metadata
print(sample_info)

##EdgeR
group <- sample_info$condition

y <- DGEList(counts = count_table, group = group)
cpms <- cpm(y)

# Filter genes expressed with CPM >= 1 in all 3 replicates of at least one group
filter_criteria <- function(row, group) {
  any(tapply(row, group, function(x) all(x >= 1)))
}

keep <- apply(cpms, 1, filter_criteria, group = group)

y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# MDS plot
plotMDS(
  y,
  col = ifelse(group == "control", "blue", "red")
)

# Experimental design
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

y <- estimateDisp(y, design, robust = TRUE)

plotBCV(y)

# Export normalized counts
cpm_tmm <- cpm(y, normalized.lib.sizes = TRUE)

cpms_with_gene_id <- cbind(
  gene_id = rownames(cpm_tmm),
  cpm_tmm
)

filename <- "CPM_TMM_IACSP-6007_clay_inoc-vs-control.tsv"

write.table(
  cpms_with_gene_id,
  file = filename,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Differential expression test
fit <- glmFit(y, design)

my.contrast <- makeContrasts(inoc - control, levels = design)

lrt <- glmLRT(fit, contrast = my.contrast)

topTags(lrt)

# DEG decision with FDR <= 0.05 and absolute logFC >= 1
DE <- decideTests(
  lrt,
  adjust.method = "BH",
  p.value = 0.05,
  lfc = 1
)

summary(DE)

DE_tags <- rownames(y)[as.logical(DE)]

# MA plot
plotSmear(lrt, de.tags = DE_tags, cex = 0.3)
abline(h = c(-1, 1), col = "blue")

# Export all expressed genes with statistics
top_tags <- topTags(lrt, n = Inf)

results_with_gene_id <- cbind(
  gene_id = rownames(top_tags$table),
  top_tags$table
)

filename <- "Expressed_IACSP-6007_clay_inoc-vs-control.tsv"

write.table(
  results_with_gene_id,
  file = filename,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Export only DEGs
threshold <- 0.05
logFC_threshold <- 1

DEGs_LogFC <- results_with_gene_id[
  results_with_gene_id$FDR < threshold &
    abs(results_with_gene_id$logFC) >= logFC_threshold,
]

filename <- "DEGs_IACSP-6007_clay_inoc-vs-control.tsv"

write.table(
  DEGs_LogFC,
  file = filename,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


```
## FOR GO TERMS ANALYSIS - WITH TOP2GO
# Load packages

```
library("tidyr")
library("dplyr")
library("stringr")

#seting working directory
setwd("working/directory/path")

```

### Preparing Gene2GO file

```
#Reading the annotation file from Blast2GO containing the gene IDs and GO terms
annot2 <- read.csv("annotation file/path/", header = T, sep = '\t')
annot2$SeqName <- sub("\\.1\\.p$", "", annot2$SeqName) 

#Selecting only gene name and correspondent GO term
gene2GO <- annot2[,c(3,10)]
colnames(gene2GO) <- c("Gene","GO")
head(gene2GO)

#Removing genes that doesn't contain a GO term
gene2GO <- gene2GO %>% filter(!is.na(GO) & GO != "")
head(gene2GO, 10)

#Separating the GO term one per row
long_gene2GO <- separate_rows(gene2GO, GO, sep = ";")
head(long_gene2GO, 10)

#remove white spaces
long_gene2GO[] <- lapply(long_gene2GO, str_trim)
head(long_gene2GO, 10)

#remove P: F: C: .1.p
long_gene2GO$GO <- gsub("P:", "", long_gene2GO$GO)
long_gene2GO$GO <- gsub("F:", "", long_gene2GO$GO)
long_gene2GO$GO <- gsub("C:", "", long_gene2GO$GO)
long_gene2GO$Gene <- sub("\\.1\\.p$", "", long_gene2GO$Gene) 
head(long_gene2GO, 10)

#Generating a list with each Gene containing a list of GO terms
gene2GO <- tapply(long_gene2GO$GO, long_gene2GO$Gene, function(x)x)
head(gene2GO, 10)

```

## Prepating DEGs and TopGO analysis

```
Load packages
library("topGO")  
library("tidyr")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("svDialogs")


#seting working directory
setwd("working/directory/path")

``` 
Choose input 

```
#Choose data to analyse
#data name
name = 'DEGs_IACSP-6007_sandy_inoc-control'

#Load csv to be analysed
DESeq2_results <- read.csv(paste0(name,".csv"), sep = "\t")

```

### Preparing DEGs 

```
head(DESeq2_results, 10)


#Removing the suffix .v2.1 from the first column gene_id
DESeq2_results$gene_id <- sub("\\.v2\\.1$", "", DESeq2_results$gene_id)
head(DESeq2_results, 10)

#separete in down and downregulated
DESeq2_results_6007_sandy_inoc_control <- DESeq2_results[order(DESeq2_results$logFC),]
DESeq2_results_up = filter(DESeq2_results, logFC > 0)
DESeq2_results_down = filter(DESeq2_results, logFC < 0)


################################## UP ########################################

#making geneList file
up_a = DESeq2_results_up %>% dplyr::select(1, 5)
names(up_a)[1] <- "SeqName"
head(up_a)

gene_list_making = annot2[3:3]
head(gene_list_making)

gene_list_making2 = full_join(up_a,gene_list_making, by = "SeqName")

up_a = gene_list_making2
up_a = up_a %>% replace(is.na(.), 1)
head(up_a)

#order
up_a <- up_a[order(up_a$PValue),]

#Creating the gene universe list containing the gene name and value 1 if significative or 0 if not
pcutoff <- 0.05 # cutoff for defining significant genes
tmp <- ifelse(up_a$PValue < pcutoff, 1, 0)

geneList <- tmp
head(geneList)

names(geneList) <- up_a$SeqName
head(geneList)

#verify the number of 0 and 1 (0 = not DEGs; 1 = DEGs)
table(geneList)

#Saving result to a csv file
write.csv(DESeq2_results_up, "DEGs_SANDY_5503-6007_CONTROL_up.csv", row.names = T)


########################## topGO analyses UP ################################


#Generating the topGO object for the Biological Process category
GOdata_up_6007_sandy <- new("topGOdata",
                            ontology = "BP",
                            description = "topGO",
                            allGenes = geneList,
                            geneSelectionFun = function(x)(x == 1),
                            annot = annFUN.gene2GO, 
                            gene2GO = gene2GO_bianca2)


#Running the topGO test weight01 fisher
resultFisher_up_w <- runTest(GOdata_up_6007_sandy, algorithm = "weight01", statistic = "fisher")


allRes_up <- GenTable(GOdata_up_6007_sandy, p.value = resultFisher_up_w,
                  orderBy = "weight01Fisher", topNodes = length(resultFisher_up_w@score), numChar = 120)


#Saving result to a csv file
write.csv(allRes, paste0("topGO_enrichment_",name,"-up_unfiltered",".csv"), row.names = F)


############################# DOWN ##########################################
#making geneList file
down_a = DESeq2_results_down %>% dplyr::select(1, 5)
names(down_a)[1] <- "SeqName"
head(down_a)

gene_list_making = annot2[3:3]
head(gene_list_making)

gene_list_making2 = full_join(down_a,gene_list_making, by = "SeqName")
head(gene_list_making2)

down_a = gene_list_making2
down_a = down_a %>% replace(is.na(.), 1)

#order
down_a <- down_a[order(down_a$PValue),]

#Creating the gene universe list containing the gene name and value 1 if significative or 0 if not
pcutoff <- 0.05 # cutoff for defining significant genes
tmp <- ifelse(down_a$PValue < pcutoff, 1, 0)

geneList <- tmp
head(geneList)

names(geneList) <- down_a$SeqName
head(geneList)

#verify the number of 0 and 1 (0 = not DEGs; 1 = DEGs)
table(tmp)
table(geneList)


#Saving result to a csv file
write.csv(DESeq2_results_down, "DEGs_SANDY_5503-6007_CONTROL_down.csv", row.names = T)


############################# topGO analyses DOWN #####################################


#Generating the topGO object for the Biological Process category
GOdata_down_6007_sandy <- new("topGOdata",
                              ontology = "BP",
                              description = "topGO",
                              allGenes = geneList,
                              geneSelectionFun = function(x)(x == 1),
                              annot = annFUN.gene2GO, 
                              gene2GO = gene2GO_bianca2)


#Running the topGO test weight01 fisher
resultFisher_down_w <- runTest(GOdata_down_6007_sandy, algorithm = "weight01", statistic = "fisher")


allRes_down <- GenTable(GOdata_down_6007_sandy, p.value = resultFisher_down_w,
                        topNodes = length(resultFisher_down_w@score), numChar = 120)

#Saving result to a csv file
write.csv(allRes_down, paste0("topGO_enrichment_",name,"-down_unfiltered",".csv"), row.names = F)

```

## Bubble plot

```
library(tidyverse)
library(dplyr)
library("viridis")      
library("patchwork")

#Load the TopGO enrichment results
data = read.csv("topGO_enrichment_DEGs_IACSP-5503_clay_inoc-control-up_unfiltered.csv",header = TRUE, sep = ,)

#################### Result visualization ########################


#Selecting the term GOs to be shown in bubble plot.
ggdata_up <- data[c(1,3,10,16,26,28,40,42,132,136,137,301,376),] #5503_clay_up

ggdata_up

colnames(ggdata_up)[6] <- "raw.p.value"
ggdata_up

#Aplicar a função à coluna
ggdata_up$raw.p.value <- gsub('<', '', ggdata_up$raw.p.value) #remove the < character
ggdata_up$raw.p.value <- as.numeric(as.character(ggdata_up$raw.p.value))
ggdata_up$Term <- factor(ggdata_up$Term, levels = rev(ggdata_up$Term)) # fixes order

#chart
bubblue_plot1 = ggplot(ggdata_up,aes(x=-log10(raw.p.value), y=Term, size=Significant, color=-log10(raw.p.value))) +
  geom_point(alpha=0.7) +
  scale_size(range = c(2, 12), name="Number of genes") +
  xlim(0,35) + #changes the X axis scale
  scale_color_viridis(option = "D") +
  ylab('') + xlab('Enrichment score') +
  #labs(title = 'GO Biological processes - 6007 down regulated DEGs',
  #subtitle = 'All terms ordered by Fisher p-value',
  #caption = '') +
  theme_bw(base_size = 12) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 12, face = 'bold', hjust = 0.9), # Define the appearance of the plot's main title
    plot.subtitle = element_text(angle = 0, size = 7, face = 'bold', vjust = 1), # Define the appearance of the plot's subtitle
    plot.caption = element_text(angle = 0, size = 6, face = 'bold', vjust = 1), # Define the appearance of the plot's caption
    
    axis.text.x = element_text(angle = 0, size = 10, face = 'bold', hjust = 1.10,colour = 'black'), # Define the appearance of the x-axis labels
    axis.text.y = element_text(angle = 0, size = 10, face = 'bold', vjust = 0.5,colour = 'black'), # Define the appearance of the y-axis labels
    #axis.title = element_text(size = 6, face = 'bold'), # General appearance for both x and y axis titles
    axis.title.x = element_text(size = 10, face = 'bold'), # Appearance specifically for x-axis title
    axis.title.y = element_text(size = 10, face = 'bold'), # Appearance specifically for y-axis title
    axis.line = element_line(colour = 'black', linewidth = 0.2), # Define the appearance of the axis lines (black color)
    axis.ticks=element_line(size=0.8),
    
    #Legend
    #legend.key = element_blank(), # removes the border
    legend.key.size = unit(0.5, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 8, face = "bold"), # Text size
    title = element_text(size = 10, face = "bold"))

bubblue_plot1


ggsave(filename = "gos_bubble.png", device = "png", plot = bubblue_plot1, width = 18, height = 18)

```

## Heatmap plot

```
library(tidyverse)
library(dplyr)
library("RColorBrewer")
library(pheatmap)

#Select GO term to build the Heatmap
go_term = "GO:0006979"

#Before building the heatmap it id needed to have all TopGO objects that will be analysed

######################## recover GO genes 5503 clay #######################


#Recover gene which has this go term from the TopGO object
go_term_up <- genesInTerm(GOdata_up_5503_clay, whichGO=go_term)
go_term_down <- genesInTerm(GOdata_down_5503_clay, whichGO=go_term)

#Append UP and DOWN in one list
if (length(go_term_up) == 1 && length(go_term_down) == 1) {
  #Combine both lists if each has exactly one element
  go_term_gene_list <- as.character(c(go_term_up[[1]], go_term_down[[1]]))
} else if (length(go_term_up) == 1) {
  #Use go_term_up if it has one element
  go_term_gene_list <- as.character(go_term_up[[1]])
} else if (length(go_term_down) == 1) {
  #Use go_term_down if it has one element
  go_term_gene_list <- as.character(go_term_down[[1]])
} else if (length(go_term_up) == 0 && length(go_term_down) == 0) {
  message <- sprintf("O termo %s não está na lista de termos GOs", go_term)
  dlgMessage(message, "ok")
}


#Recover the logFC from each gene

#load data
name = 'DEGs_IACSP-5503_clay_inoc-control'
DESeq2_results <- read.csv(
  paste0(name,".csv"), sep = "\t")

#Removing the suffix .v2.1 from the first column gene_id
DESeq2_results$gene_id <- sub("\\.v2\\.1$", "", DESeq2_results$gene_id)

DESeq2_results_5503_clay_inoc_control <- DESeq2_results[order(DESeq2_results$logFC),]

go_term_gene_list_5503_clay <- DESeq2_results_5503_clay_inoc_control %>%
  filter(gene_id %in% go_term_gene_list)


######################## recover GO genes 5503 sandy #######################


#Recover gene which has this go term from the TopGO object
go_term_up_2 <- genesInTerm(GOdata_up_5503_sandy, whichGO=go_term)
go_term_down_2 <- genesInTerm(GOdata_down_5503_sandy, whichGO=go_term)

#Append UP and DOWN in one list
if (length(go_term_up_2) == 1 && length(go_term_down_2) == 1) {
  #Combine both lists if each has exactly one element
  go_term_gene_list <- as.character(c(go_term_up_2[[1]], go_term_down_2[[1]]))
} else if (length(go_term_up_2) == 1) {
  #Use go_term_up if it has one element
  go_term_gene_list <- as.character(go_term_up_2[[1]])
} else if (length(go_term_down_2) == 1) {
  #Use go_term_down if it has one element
  go_term_gene_list <- as.character(go_term_down_2[[1]])
} else if (length(go_term_up_2) == 0 && length(go_term_down_2) == 0) {
  message <- sprintf("O termo %s não está na lista de termos GOs", go_term)
  dlgMessage(message, "ok")
}


#Recover the logFC from each gene

#load data
name = 'DEGs_IACSP-5503_sandy_inoc-control'
DESeq2_results <- read.csv(
  paste0(name,".csv"), sep = "\t")

#Removing the suffix .v2.1 from the first column gene_id
DESeq2_results$gene_id <- sub("\\.v2\\.1$", "", DESeq2_results$gene_id)

DESeq2_results_5503_sandy_inoc_control <- DESeq2_results[order(DESeq2_results$logFC),]

go_term_gene_list_5503_sandy <- DESeq2_results_5503_sandy_inoc_control %>%
  filter(gene_id %in% go_term_gene_list)

######################## recover GO genes 6007 clay #######################


#Recover gene which has this go term from the TopGO object
go_term_up <- genesInTerm(GOdata_up_6007_clay, whichGO=go_term)
go_term_down <- genesInTerm(GOdata_down_6007_clay, whichGO=go_term)

#Append UP and DOWN in one list
if (length(go_term_up) == 1 && length(go_term_down) == 1) {
  #Combine both lists if each has exactly one element
  go_term_gene_list <- as.character(c(go_term_up[[1]], go_term_down[[1]]))
} else if (length(go_term_up) == 1) {
  #Use go_term_up if it has one element
  go_term_gene_list <- as.character(go_term_up[[1]])
} else if (length(go_term_down) == 1) {
  #Use go_term_down if it has one element
  go_term_gene_list <- as.character(go_term_down[[1]])
} else if (length(go_term_up) == 0 && length(go_term_down) == 0) {
  message <- sprintf("O termo %s não está na lista de termos GOs", go_term)
  dlgMessage(message, "ok")
}


#Recover the logFC from each gene

#load data
name = 'DEGs_IACSP-6007_clay_inoc-control'
DESeq2_results <- read.csv(
  paste0(name,".csv"), sep = "\t")

#Removing the suffix .v2.1 from the first column gene_id
DESeq2_results$gene_id <- sub("\\.v2\\.1$", "", DESeq2_results$gene_id)

DESeq2_results_6007_clay_inoc_control <- DESeq2_results[order(DESeq2_results$logFC),]


go_term_gene_list_6007_clay <- DESeq2_results_6007_clay_inoc_control %>%
  filter(gene_id %in% go_term_gene_list)


######################## recover GO genes 6007 sandy #######################


#Recover gene which has this go term from the TopGO object
go_term_up_2 <- genesInTerm(GOdata_up_6007_sandy, whichGO=go_term)
go_term_down_2 <- genesInTerm(GOdata_down_6007_sandy, whichGO=go_term)

#Append UP and DOWN in one list
if (length(go_term_up_2) == 1 && length(go_term_down_2) == 1) {
  #Combine both lists if each has exactly one element
  go_term_gene_list <- as.character(c(go_term_up_2[[1]], go_term_down_2[[1]]))
} else if (length(go_term_up_2) == 1) {
  #Use go_term_up if it has one element
  go_term_gene_list <- as.character(go_term_up_2[[1]])
} else if (length(go_term_down_2) == 1) {
  #Use go_term_down if it has one element
  go_term_gene_list <- as.character(go_term_down_2[[1]])
} else if (length(go_term_up_2) == 0 && length(go_term_down_2) == 0) {
  message <- sprintf("O termo %s não está na lista de termos GOs", go_term)
  dlgMessage(message, "ok")
}


#Recover the logFC from each gene

#load data
name = 'DEGs_IACSP-6007_sandy_inoc-control'
DESeq2_results <- read.csv(
  paste0(name,".csv"), sep = "\t")

#Removing the suffix .v2.1 from the first column gene_id
DESeq2_results$gene_id <- sub("\\.v2\\.1$", "", DESeq2_results$gene_id)

DESeq2_results_6007_sandy_inoc_control <- DESeq2_results[order(DESeq2_results$logFC),]

go_term_gene_list_6007_sandy <- DESeq2_results_6007_sandy_inoc_control %>%
  filter(gene_id %in% go_term_gene_list)

######################## Building the heatmap matrix #####################################

go_term_gene_list_merge = rbind(go_term_gene_list_5503_clay, 
                                go_term_gene_list_5503_sandy, 
                                go_term_gene_list_6007_clay, 
                                go_term_gene_list_6007_sandy)

go_term_gene_list_unique = unique(go_term_gene_list_merge$gene_id)



#Build a matrix to storage the logFC
go_matrix = matrix(NA, nrow = length(unique(go_term_gene_list_merge$gene_id)), ncol = 4,
                   dimnames = list(unique(go_term_gene_list_merge$gene_id), c("5503_clay", "5503_sandy",
                                                                              "6007_clay", "6007_sandy")))

#FIll the matrix with the logFC valeus
go_matrix[go_term_gene_list_5503_clay$gene_id, "5503_clay"] = (go_term_gene_list_5503_clay$logFC)
go_matrix[go_term_gene_list_5503_sandy$gene_id, "5503_sandy"] = (go_term_gene_list_5503_sandy$logFC)
go_matrix[go_term_gene_list_6007_clay$gene_id, "6007_clay"] = (go_term_gene_list_6007_clay$logFC)
go_matrix[go_term_gene_list_6007_sandy$gene_id, "6007_sandy"] = (go_term_gene_list_6007_sandy$logFC)

#Replace NA for 0
go_matrix[is.na(go_matrix)] = 0

#Build the heatmap colour palette
myColor <- colorRampPalette(c("blue", "white", "red"))(50)

myBreaks <- c(seq(min(go_matrix), 0, length.out=ceiling(50/2) + 1), 
              seq(max(go_matrix)/50, max(go_matrix), length.out=floor(50/2)))

#Recobe the gene annotation from the annotation file
gene_name_annot = annot2 %>% #annot2 is the annotation file
  distinct(SeqName, .keep_all = T) %>% 
  filter(SeqName %in% go_term_gene_list_unique)

gene_annot_df = data.frame(SeqName = go_term_gene_list_unique)

gene_annot_df <- gene_annot_df %>%
  left_join(gene_name_annot %>% dplyr::select(SeqName, Description), by = "SeqName")


row.names(gene_annot_df) <- gene_annot_df$SeqName


gene_annot_df <- gene_annot_df %>% mutate_all(na_if,"")


#Build the heatmap column annotation
ann_df <- data.frame( Group = c("5503_clay","5503_sandy","6007_clay","6007_sandy"))
row.names(ann_df) <- colnames(go_matrix)

ann_colors <- list(
  Group = c("5503_clay" = "darkolivegreen",
            "5503_sandy" = "darkolivegreen3",
            "6007_clay" = "violetred3",
            "6007_sandy" = "violetred1"))

main =  paste0("Innoc-Control ",go_term)

#Summarize by gene annotation. 
#Gene with the same annotation are treat as one gene and the logFC is the average of the genes with the same annotation
go_matrix2 = as.data.frame(go_matrix)
go_matrix2$SeqName <- rownames(go_matrix2)


go_matrix2 <- go_matrix2 %>%
  left_join(gene_name_annot %>% dplyr::select(SeqName, Description), by = "SeqName")

colnames(go_matrix2) <- c("cinco_clay","cinco_sandy","seis_clay","seis_sandy","genes","Description")

go_matrix3 = go_matrix2 %>%
  group_by(Description) %>% summarize("5503_clay" = mean(cinco_clay), 
                                      "5503_sandy" = mean(cinco_sandy),
                                      "6007_clay" = mean(seis_clay), 
                                      "6007_sandy" = mean(seis_sandy))

go_matrix3 = go_matrix3 %>% remove_rownames %>% column_to_rownames(var="Description")

go_matrix3 = as.matrix(go_matrix3)
go_matrix3

#Build heatmap for all genes
heatmap1 = pheatmap::pheatmap(go_matrix3, 
                              clustering_distance_cols = 'euclidean',
                              clustering_distance_rows = "euclidean",
                              clustering_method_rows = "euclidean",
                              cluster_rows = T,               # Clusterizar as linhas
                              cluster_cols = FALSE,              # Desativar clusterização das colunas
                              show_rownames = TRUE,              # Mostrar nomes das linhas (termos GO)
                              show_colnames = F,              # Mostrar nomes das colunas (Up/Down)
                              color = myColor,  # Paleta de cores
                              border_color = "gray",            # Cor da borda das células
                              fontsize_row = 10,                 # Tamanho da fonte das linhas
                              fontsize_col = 10,                 # Tamanho da fonte das colunas
                              angle_col = "0",
                              main = main,
                              breaks = myBreaks,
                              annotation_col = ann_df,  # Adicionar anotações de coluna
                              annotation_colors = ann_colors,   # Definir cores das anotações
)

heatmap1

######################## Building matrix for selected genes ######################

#Build matrix
selected_genes_df = data.frame(go_matrix3)
selected_genes_df

selected_genes_matrix <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(selected_genes_matrix) <- c("","IAC-5503_clay","IAC-5503_sandy","IAC-6007_clay","IAC-6007_sandy")
selected_genes_matrix

########################## Rbinding the selected genes ##########################
selected_genes_df = data.frame(go_matrix3)

#Visualise the row index of the genes
v <- rownames(go_matrix3)
v

#A list with the selected row index number of the selected genes
a = c(53:61,39,20,1:3,14,13,46,52,73)
a

#'for' loop in the list of number of the selected genes
for (i in a){
  selected_genes_matrix <- rbind(selected_genes_matrix,selected_genes_df[i,])
}
selected_genes_matrix    # selected genes matrix


#Saving result to a csv file
write.csv(selected_genes_matrix, "gene_selected_inoc-control-0006979.csv", row.names = T)



##########################################################

selected_genes_matrix <- as.matrix(selected_genes_matrix)
selected_genes_matrix

go_matrix3 = selected_genes_matrix

##########################################################

rg <- max(abs(matrix))

#Build the heatmap colour palette
myColor <- colorRampPalette(c("blue", "white", "red"))(100)


myBreaks <- c(seq(min(selected_genes_matrix), 0, length.out=ceiling(100/2)+1),
              seq(max(selected_genes_matrix)/100, max(selected_genes_matrix), length.out=floor(100/2)))


#Build the heatmap column annotation
ann_df <- data.frame( Group = c("IAC-5503_clay","IAC-5503_sandy","IAC-6007_clay","IAC-6007_sandy"))
row.names(ann_df) <- colnames(selected_genes_matrix)

ann_colors <- list(
  Group = c("IAC-5503_clay" = "darkolivegreen",
            "IAC-5503_sandy" = "darkolivegreen3",
            "IAC-6007_clay" = "violetred3",
            "IAC-6007_sandy" = "violetred1"))    


main =  paste0("Inoc-Control selected ", go_term )


#Build heatmap for selected genes                                                                                                                    
heatmap2 = pheatmap::pheatmap(go_matrix3, 
                              clustering_distance_cols = 'euclidean',
                              clustering_distance_rows = "euclidean",
                              clustering_method_rows = "euclidean",,
                              cluster_rows = T,               # Clusterizar as linhas
                              cluster_cols = FALSE,              # Desativar clusterização das colunas
                              show_rownames = TRUE,              # Mostrar nomes das linhas (termos GO)
                              show_colnames = F,              # Mostrar nomes das colunas (Up/Down)
                              color = myColor,  # Paleta de cores
                              border_color = "gray",            # Cor da borda das células
                              fontsize_row = 10,                 # Tamanho da fonte das linhas
                              fontsize_col = 10,                 # Tamanho da fonte das colunas
                              angle_col = "0",
                              main = main,
                              breaks = myBreaks,
                              annotation_col = ann_df,  # Adicionar anotações de coluna
                              annotation_colors = ann_colors,   # Definir cores das anotações
                             
)

heatmap2

```

### Plots: Heatmap With Expressed Genes

```
rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(pheatmap)
  library(grid)
})

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_dir)

# -----------------------------------------------------------------------------
# PURPOSE
# -----------------------------------------------------------------------------
# Rebuild the selected heatmaps in the same pheatmap style used before, but now
# filling the matrix with logFC values from ALL treatments for the same selected
# genes/descriptions.
#
# Selection logic:
# - the old gene_selected_inoc-control-*.csv files define which rows appear
#   and in which order
# - the values are reloaded from the full Expressed_* tables
#
# Differential-expression direction:
#   logFC = inoculated - control
#   logFC > 0  -> higher in inoculated
#   logFC < 0  -> higher in control

# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------
selected_files_pattern <- "^gene_selected_inoc-control-[0-9]+\\.csv$"

contrast_metadata <- tibble::tribble(
  ~contrast_id,   ~file_name,                                      ~soil,   ~genotype,
  "5503_clay",    "C:/Users/joyde/OneDrive/Área de Trabalho/Desktop NIOO/RNAseq/4EdgeR/INOC-CONTROL/IACSP-5503_CLAY/Expressed_IACSP-5503_clay_inoc-control.csv", "clay",  "IACSP-5503",
  "5503_sandy",   "C:/Users/joyde/OneDrive/Área de Trabalho/Desktop NIOO/RNAseq/4EdgeR/INOC-CONTROL/IACSP-5503_SANDY/Expressed_IACSP-5503_inoc-control.csv", "sandy", "IACSP-5503",
  "6007_clay",    "C:/Users/joyde/OneDrive/Área de Trabalho/Desktop NIOO/RNAseq/4EdgeR/INOC-CONTROL/IACSP-6007_CLAY/Expressed_IACSP-6007_clay_inoc-vs-control.tsv", "clay",  "IACSP-6007",
  "6007_sandy",   "C:/Users/joyde/OneDrive/Área de Trabalho/Desktop NIOO/RNAseq/4EdgeR/INOC-CONTROL/IACSP-6007_SANDY/Expressed_IACSP-6007_sandy_inoc-control.csv", "sandy", "IACSP-6007"
)

column_order <- c("5503_clay", "5503_sandy", "6007_clay", "6007_sandy")

soil_colors <- c(clay = "#22C7D5", sandy = "#9AD500")
genotype_colors <- c("IACSP-5503" = "#B77AF5", "IACSP-6007" = "#FF8D8A")

cluster_rows <- TRUE
cluster_cols <- TRUE
show_colnames <- FALSE
show_rownames <- TRUE
fontsize_row <- 9
fontsize_col <- 9
cellwidth <- 26
cellheight <- 16
border_color <- "grey88"
na_color <- "grey95"

color_limits <- c(-1.5, 1.5)
palette_n <- 100

output_dir <- file.path(project_dir, "selected_terms_all_treatments_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# HELPERS
# -----------------------------------------------------------------------------
guess_delim <- function(path) {
  first_line <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
  counts <- c(
    "," = stringr::str_count(first_line, ","),
    ";" = stringr::str_count(first_line, ";"),
    "\t" = stringr::str_count(first_line, "\t")
  )
  names(which.max(counts))
}

read_table_auto <- function(path) {
  readr::read_delim(
    path,
    delim = guess_delim(path),
    show_col_types = FALSE,
    trim_ws = TRUE
  )
}

sanitize_go_id <- function(go_id) {
  stringr::str_replace_all(go_id, "[: ]", "")
}

normalize_description <- function(x) {
  x %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[;,]", " ") %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_trim()
}

discover_selected_files <- function() {
  files <- list.files(
    path = project_dir,
    pattern = selected_files_pattern,
    full.names = TRUE
  )

  if (length(files) == 0) {
    stop("No selected-gene files found with pattern: ", selected_files_pattern)
  }

  files
}

parse_go_dictionary <- function(annotation_file) {
  ann <- read_table_auto(annotation_file)

  ann %>%
    transmute(
      go_ids = as.character(`GO IDs`),
      go_names = as.character(`GO Names`)
    ) %>%
    filter(!is.na(go_ids), go_ids != "", !is.na(go_names), go_names != "") %>%
    separate_rows(go_ids, go_names, sep = ";") %>%
    mutate(
      go_ids = stringr::str_trim(go_ids),
      go_names = stringr::str_trim(go_names),
      go_ids = stringr::str_remove(go_ids, "^[PCF]:"),
      go_names = stringr::str_remove(go_names, "^[PCF]:")
    ) %>%
    filter(go_ids != "", go_names != "") %>%
    group_by(go_ids, go_names) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    arrange(go_ids, desc(n), go_names) %>%
    group_by(go_ids) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    transmute(go_id = go_ids, go_name = go_names)
}

load_annotation_table <- function(annotation_file) {
  ann <- read_table_auto(annotation_file)

  required_cols <- c("SeqName", "Description", "GO IDs")
  missing_cols <- setdiff(required_cols, colnames(ann))
  if (length(missing_cols) > 0) {
    stop("Missing annotation columns: ", paste(missing_cols, collapse = ", "))
  }

  ann %>%
    transmute(
      gene_id = sub("\\.1\\.p$", "", SeqName),
      Description = as.character(Description),
      GO_IDs = as.character(`GO IDs`)
    ) %>%
    distinct(gene_id, .keep_all = TRUE)
}

read_expressed_table <- function(path, contrast_id) {
  df <- read_table_auto(path)

  required_cols <- c("gene_id", "logFC", "FDR")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns in ", basename(path), ": ", paste(missing_cols, collapse = ", "))
  }

   df <- df %>%
    mutate(
      logFC = as.numeric(logFC),
      FDR = as.numeric(FDR)
    )

  # Sanity checks to catch corrupted exports before plotting.
  if (any(!is.na(df$FDR) & (df$FDR < 0 | df$FDR > 1))) {
    stop("Invalid FDR values detected in ", basename(path), ". The file may be corrupted.")
  }

  if ("PValue" %in% colnames(df)) {
    df <- df %>% mutate(PValue = as.numeric(PValue))
    if (any(!is.na(df$PValue) & (df$PValue < 0 | df$PValue > 1))) {
      stop("Invalid PValue values detected in ", basename(path), ". The file may be corrupted.")
    }
  }

  if (any(!is.na(df$logFC) & abs(df$logFC) > 50)) {
    stop("Unusually large |logFC| values detected in ", basename(path), ". The file may be corrupted.")
  }

  df %>%
    transmute(
      gene_id = sub("\\.v2\\.1$", "", gene_id),
      logFC = logFC,
      FDR = FDR,
      is_sig = !is.na(FDR) & FDR <= 0.05,
      contrast_id = contrast_id
    )
}

load_all_expressed_long <- function(annotation_df) {
  purrr::map2_dfr(
    contrast_metadata$file_name,
    contrast_metadata$contrast_id,
    ~ read_expressed_table(
      if (grepl("^[A-Za-z]:[/\\\\]", .x)) .x else file.path(project_dir, .x),
      .y
    )
  ) %>%
    left_join(annotation_df, by = "gene_id") %>%
    left_join(contrast_metadata, by = "contrast_id")
}

read_selected_row_order <- function(selected_file) {
  delim <- guess_delim(selected_file)

  if (delim == ";") {
    lines <- readLines(selected_file, warn = FALSE, encoding = "UTF-8")
    if (length(lines) <= 1) {
      return(character(0))
    }

    row_lab <- vapply(lines[-1], function(line) {
      parts <- strsplit(line, ";", fixed = TRUE)[[1]]

      while (length(parts) > 0 && parts[length(parts)] == "") {
        parts <- parts[-length(parts)]
      }

      if (length(parts) <= 4) {
        return(NA_character_)
      }

      desc <- paste(parts[seq_len(length(parts) - 4)], collapse = ";")
      desc <- stringr::str_trim(desc)
      desc <- stringr::str_remove(desc, '^"')
      desc <- stringr::str_remove(desc, '"$')
      desc
    }, character(1))

    row_lab <- row_lab[!is.na(row_lab) & row_lab != ""]
    return(unique(as.character(row_lab)))
  }

  df <- read.csv(
    selected_file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  row_lab <- df[[1]]
  row_lab <- row_lab[!is.na(row_lab) & row_lab != ""]
  as.character(row_lab)
}

extract_go_gene_ids <- function(annotation_df, go_id_current) {
  annotation_df %>%
    filter(!is.na(GO_IDs), GO_IDs != "") %>%
    tidyr::separate_rows(GO_IDs, sep = ";") %>%
    mutate(
      GO_IDs = stringr::str_trim(GO_IDs),
      GO_IDs = stringr::str_remove(GO_IDs, "^[PCF]:")
    ) %>%
    filter(GO_IDs == go_id_current) %>%
    distinct(gene_id) %>%
    pull(gene_id)
}

rebuild_matrices_for_selected_rows <- function(selected_rows, expressed_long) {
  selected_df <- tibble(
    Description = selected_rows,
    norm_description = normalize_description(selected_rows),
    row_id = make.unique(selected_rows, sep = "___dup")
  )

  summarized <- expressed_long %>%
    filter(
      !is.na(Description),
      normalize_description(Description) %in% selected_df$norm_description
    ) %>%
    mutate(norm_description = normalize_description(Description)) %>%
    filter(norm_description %in% selected_df$norm_description) %>%
    group_by(norm_description, contrast_id) %>%
    summarise(
      logFC = mean(logFC, na.rm = TRUE),
      is_sig = any(is_sig, na.rm = TRUE),
      n_gene_ids = dplyr::n_distinct(gene_id),
      .groups = "drop"
    ) %>%
    mutate(logFC = ifelse(is.nan(logFC), NA_real_, logFC))

  full_df <- tidyr::expand_grid(
    row_id = selected_df$row_id,
    contrast_id = column_order
  ) %>%
    left_join(selected_df, by = "row_id") %>%
    left_join(summarized, by = c("norm_description", "contrast_id"))

  mat_df <- full_df %>%
    select(row_id, contrast_id, logFC) %>%
    pivot_wider(names_from = contrast_id, values_from = logFC) %>%
    as.data.frame(stringsAsFactors = FALSE)

  rownames(mat_df) <- mat_df$row_id
  mat_df$row_id <- NULL
  mat_df <- mat_df[, column_order, drop = FALSE]
  mat_df[] <- lapply(mat_df, as.numeric)

  mat <- as.matrix(mat_df)
  storage.mode(mat) <- "numeric"

  sig_df <- full_df %>%
    select(row_id, contrast_id, is_sig) %>%
    mutate(is_sig = ifelse(is.na(is_sig), FALSE, is_sig)) %>%
    pivot_wider(names_from = contrast_id, values_from = is_sig) %>%
    as.data.frame(stringsAsFactors = FALSE)

  rownames(sig_df) <- sig_df$row_id
  sig_df$row_id <- NULL
  sig_df <- sig_df[, column_order, drop = FALSE]
  sig_df[] <- lapply(sig_df, as.logical)

  sig_mat <- as.matrix(sig_df)
  storage.mode(sig_mat) <- "logical"

  list(
    logfc_mat = mat[selected_df$row_id, column_order, drop = FALSE],
    sig_mat = sig_mat[selected_df$row_id, column_order, drop = FALSE],
    row_labels = selected_df$Description
  )
}

build_annotation_col <- function() {
  ann_col <- contrast_metadata %>%
    select(contrast_id, soil, genotype) %>%
    as.data.frame(stringsAsFactors = FALSE)

  rownames(ann_col) <- ann_col$contrast_id
  ann_col$contrast_id <- NULL
  ann_col[column_order, , drop = FALSE]
}

build_breaks <- function() {
  c(
    seq(color_limits[1], 0, length.out = ceiling(palette_n / 2) + 1),
    seq(color_limits[2] / palette_n, color_limits[2], length.out = floor(palette_n / 2))
  )
}

build_safe_hclust <- function(mat, margin = c("row", "col"), method = "complete") {
  margin <- match.arg(margin)
  mat2 <- mat

  if (!is.matrix(mat2)) {
    mat2 <- as.matrix(mat2)
  }

  storage.mode(mat2) <- "numeric"

  if (margin == "row") {
    if (nrow(mat2) < 2) {
      return(FALSE)
    }

    row_means <- apply(mat2, 1, function(x) {
      if (all(is.na(x))) 0 else mean(x, na.rm = TRUE)
    })

    for (i in seq_len(nrow(mat2))) {
      idx <- is.na(mat2[i, ])
      if (any(idx)) {
        mat2[i, idx] <- row_means[i]
      }
    }

    return(hclust(dist(mat2), method = method))
  }

  if (ncol(mat2) < 2) {
    return(FALSE)
  }

  col_means <- apply(mat2, 2, function(x) {
    if (all(is.na(x))) 0 else mean(x, na.rm = TRUE)
  })

  for (j in seq_len(ncol(mat2))) {
    idx <- is.na(mat2[, j])
    if (any(idx)) {
      mat2[idx, j] <- col_means[j]
    }
  }

  hclust(dist(t(mat2)), method = method)
}

draw_one_heatmap <- function(logfc_mat, sig_mat, row_labels, go_id, go_name, output_png, output_pdf) {
  if (!is.matrix(logfc_mat)) {
    logfc_mat <- as.matrix(logfc_mat)
  }
  storage.mode(logfc_mat) <- "numeric"

  display_numbers <- matrix("", nrow = nrow(sig_mat), ncol = ncol(sig_mat))
  display_numbers[sig_mat] <- "*"
  rownames(display_numbers) <- rownames(sig_mat)
  colnames(display_numbers) <- colnames(sig_mat)

  ann_col <- build_annotation_col()
  ann_colors <- list(
    soil = soil_colors,
    genotype = genotype_colors
  )

  myColor <- colorRampPalette(c("blue", "white", "red"))(palette_n)
  myBreaks <- build_breaks()

  row_cluster_obj <- if (isTRUE(cluster_rows)) {
    build_safe_hclust(logfc_mat, margin = "row", method = "complete")
  } else {
    FALSE
  }

  col_cluster_obj <- if (isTRUE(cluster_cols)) {
    build_safe_hclust(logfc_mat, margin = "col", method = "complete")
  } else {
    FALSE
  }

  title_text <- paste0(go_id, ' - "', go_name, '"')

  png(
    filename = output_png,
    width = 3200,
    height = max(2200, 1200 + nrow(logfc_mat) * 45),
    res = 300,
    bg = "white"
  )
  pheatmap::pheatmap(
    mat = logfc_mat,
    cluster_rows = row_cluster_obj,
    cluster_cols = col_cluster_obj,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    color = myColor,
    border_color = border_color,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    angle_col = 0,
    main = title_text,
    breaks = myBreaks,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    na_col = na_color,
    labels_row = row_labels,
    display_numbers = display_numbers,
    number_color = "black",
    fontsize_number = 11,
    cellwidth = cellwidth,
    cellheight = cellheight,
    silent = FALSE
  )
  grid.text(
    "* = significant DEG | logFC = inoculated - control",
    x = 0.98, y = 0.02, just = c("right", "bottom"),
    gp = gpar(fontsize = 10)
  )
  dev.off()

  pdf(
    file = output_pdf,
    width = 10.5,
    height = max(7, 3 + nrow(logfc_mat) * 0.18),
    bg = "white",
    useDingbats = FALSE
  )
  pheatmap::pheatmap(
    mat = logfc_mat,
    cluster_rows = row_cluster_obj,
    cluster_cols = col_cluster_obj,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    color = myColor,
    border_color = border_color,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    angle_col = 0,
    main = title_text,
    breaks = myBreaks,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    na_col = na_color,
    labels_row = row_labels,
    display_numbers = display_numbers,
    number_color = "black",
    fontsize_number = 11,
    cellwidth = cellwidth,
    cellheight = cellheight,
    silent = FALSE
  )
  grid.text(
    "* = significant DEG | logFC = inoculated - control",
    x = 0.98, y = 0.02, just = c("right", "bottom"),
    gp = gpar(fontsize = 10)
  )
  dev.off()
}

# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------
annotation_file <- file.path(project_dir, "omicsbox_table_all_proteins.txt")
annotation_df <- load_annotation_table(annotation_file)
expressed_long <- load_all_expressed_long(annotation_df)
selected_files <- discover_selected_files()
go_dictionary <- parse_go_dictionary(annotation_file)

summary_rows <- lapply(selected_files, function(selected_file) {
  go_num <- stringr::str_extract(basename(selected_file), "[0-9]+")
  go_id_current <- paste0("GO:", go_num)

  go_name <- go_dictionary %>%
    filter(go_id == go_id_current) %>%
    pull(go_name)

  if (length(go_name) == 0) {
    go_name <- go_id_current
  } else {
    go_name <- go_name[1]
  }

  selected_rows <- read_selected_row_order(selected_file)
  mats <- rebuild_matrices_for_selected_rows(selected_rows, expressed_long)

  safe_go <- sanitize_go_id(go_id_current)

  write.csv(
    data.frame(row_id = rownames(mats$logfc_mat), Description = mats$row_labels, mats$logfc_mat, check.names = FALSE),
    file = file.path(output_dir, "tables", paste0("heatmap_matrix_", safe_go, ".csv")),
    row.names = FALSE
  )
  write.csv(
    data.frame(row_id = rownames(mats$sig_mat), Description = mats$row_labels, mats$sig_mat, check.names = FALSE),
    file = file.path(output_dir, "tables", paste0("significance_matrix_", safe_go, ".csv")),
    row.names = FALSE
  )

  draw_one_heatmap(
    logfc_mat = mats$logfc_mat,
    sig_mat = mats$sig_mat,
    row_labels = mats$row_labels,
    go_id = go_id_current,
    go_name = go_name,
    output_png = file.path(output_dir, "png", paste0("heatmap_", safe_go, ".png")),
    output_pdf = file.path(output_dir, "pdf", paste0("heatmap_", safe_go, ".pdf"))
  )

    tibble(
      go_id = go_id_current,
      go_name = go_name,
      n_rows = nrow(mats$logfc_mat),
      source_file = basename(selected_file),
      status = "ok"
    )
})

summary_df <- bind_rows(summary_rows)
write_csv(summary_df, file.path(output_dir, "heatmap_summary.csv"))

writeLines(
  c(
    "Methodology",
    "1. Existing gene_selected_inoc-control-*.csv files were used only to define row selection and row order.",
    "2. logFC values were reloaded from the full Expressed_* treatment tables.",
    "3. The same selected genes/descriptions are shown, but all treatments are filled when values exist.",
    "4. logFC direction is preserved as inoculated - control.",
    "5. Positive values indicate higher expression in inoculated plants.",
    "6. Negative values indicate higher expression in control plants."
  ),
  con = file.path(output_dir, "methodology.txt")
)

message("Selected-term heatmaps with all treatments were generated.")

```

