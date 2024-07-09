# :dna: Microbiome analysis and transcriptional profile of sugarcane varieties in resistance to the fungus _Sporisorium scitamineum_

This repository will be used to describe, in detail, the data processing and results related to the _Sporisorium scitamineum_ microbiome project. The experiment consists of two varieties of sugar cane in low soils
and high fertility. The seedlings were planted in December 2022 at Usina Iracema, in São Martinho (SP), and the plants were  kept in the field for seven months. Experiments were arranged in random blocks with six replicates.
Each experimental plot consists of four rows measuring four meters in length and planting spacing of 0.5 m, totaling 32 plants per plot. As varieties chosen are IACSP01-5503 (adapted to low fertility soils and
resistant to _S. scitamineum_) and IACSP-6007 (not adapted to low fertility soils and resistant to _S. scitamineum_).

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
**Kraken2: 2.1.3**

```
mamba install -c bioconda hisat2 samtools subread kraken2
```

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

# Verifica se o diretório de mapeamento já existe, se não, cria
if [ ! -d ../../Data/2Alignments_R570 ]; then
  mkdir -p ../../Data/2Alignments_R570
fi

# Cria diretórios para leituras mapeadas e não mapeadas
mapped_dir="../../Data/2Alignments_R570/R570_mapped"
unmapped_dir="../../Data/2Alignments_R570/R570_unmapped"
mkdir -p "$mapped_dir"
mkdir -p "$unmapped_dir"

# Define o diretório de trabalho
cd ../../Data/2Alignments_R570

# Caminho para o índice do HISAT2
index="../References/R570_poliploid_version_2024/hisat2_index/SofficinarumxspontaneumR570_771_v2.0"

# Nome do arquivo de controle
checkpoint_file="./processed_files.txt"

# Verifica se o arquivo de controle existe, se não existir, cria
if [ ! -e "$checkpoint_file" ]; then
  touch "$checkpoint_file"
fi

# Função para lidar com o sinal de interrupção (Ctrl+C)
function cleanup {
  echo "Script interrupted. Exiting..."
  exit 1
}

# Captura o sinal de interrupção (Ctrl+C)
trap cleanup SIGINT

# Loop para processar arquivos
for PE1 in /home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data/rRNA_clean/Ssci_free/*PE1.fastq.gz; do
  PE2="${PE1/PE1/PE2}"
  BASE_PE1=$(basename "$PE1")

# Verifica se o par já foi processado antes de executar o HISAT2
  if ! grep -q "$BASE_PE1" "$checkpoint_File"; then
    echo "Analyzing $PE1 and $PE2"

   hisat2 -p 20 --rg-id "$PE1" --rg SM:"$PE1" \
    --summary-file ./summary_"$BASE_PE1"_Ssci_free.txt \
    -x "$index" -1 "$PE1" -2 "$PE2" \
    -S ./"$BASE_PE1"_Ssci_free.sam

    # Ordena o arquivo SAM e converte para BAM
    samtools sort -@ 20 ./"$BASE_PE1"_Ssci_free.sam > ./"$BASE_PE1"_Ssci_free.bam

    # Devido ao tamanho de cromossomo maior que 512Mb, não é possível criar o bam index do tipo .bai; usamos então o parâmetro "-c"
    samtools index -@ 20 -c ./"$BASE_PE1"_Ssci_free.bam
    rm ./"$BASE_PE1"_Ssci_free.sam

   # Filtra as leituras mapeadas e as move para o diretório de leituras mapeadas
 samtools view --threads 20 -b -F 4 ./"$BASE_PE1"_Ssci_free.bam > "$mapped_dir/${BASE_PE1}_R570mapped.bam"

    # Filtra as leituras não mapeadas e as move para o diretório de leituras não mapeadas
    samtools view --threads 20 -b -f 4 ./"$BASE_PE1"_Ssci_free.bam > "$unmapped_dir/${BASE_PE1}_R570unmapped.bam"

    # Adiciona o nome do arquivo ao arquivo de controle
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

# Define os caminhos dos arquivos
input_gff3="../../Data/References/R570_poliploid_version_2024/annotation/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gff3.gz"
output_gtf="../../Data/References/R570_poliploid_version_2024/annotation/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf"

# Verifica se o arquivo de entrada está compactado
if [[ "$input_gff3" == *.gz ]]; then
    echo "Descompactando $input_gff3..."
    gunzip -c "$input_gff3" > "${input_gff3%.gz}"
    input_gff3="${input_gff3%.gz}"
fi

# Converte o GFF3 para GTF
gffread "$input_gff3" -T -o "$output_gtf"

# Verifica se a conversão foi bem-sucedida
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

# Verifica se o diretório de saída existe, caso contrário, cria
if [ ! -d ../../Data/3featureCounts_R570__ShSHN ]; then
    mkdir -p ../../Data/3featureCounts_R570__ShSHN
fi

# Muda para o diretório de saída
cd ../../Data/3featureCounts_R570__ShSHN/

# Executa o featureCounts com os parâmetros fornecidos
featureCounts -s 0 -p --countReadPairs -T 20 -t CDS -g gene_id -M \
-a /home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data/References/R570_poliploid_version_2024/annotation/SofficinarumxspontaneumR570_771_v2.1.gene_exons.gtf \
-o ./Quantified_all_R570_genome_ShSHN_ref.tsv \
$(ls /home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data/2Alignments_R570/R570_mapped/*.bam) \
2> ./featureCounts_R570_genome_ShSHN_log.txt


```

## EdgeR - RStudio
Until now, we have carried out all analyzes on the server. From now on we will work through RStudio for the following analyses.


```
library(dplyr)
library(edgeR)

#Load count table
count_table <- read.table("/home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data/Cutadapt_alignment_input/EdgeR/Joyce_edgeR/Quantified_all_R570_genome_ShSHN_ref.tsv", header = TRUE, row.names = "Geneid")
gene_length <- count_table[,5]
names(gene_length) <- rownames(count_table)
count_table <- count_table[, -(1:5)]

#Load correspondance file
correspondence <- read.table("/home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data/Cutadapt_alignment_input/EdgeR/Joyce_edgeR/samples.txt", header = FALSE)

# Load correspondence table
### This way its ensured that the names of the samples and sequencing filenames matches, and its not just a mistake regarding the column order
colnames(correspondence) <- c("old_name", "new_name")

# Create a mapping between old and new column names
mapping <- setNames(correspondence$new_name, correspondence$old_name)

# Process column names
new_names <- sapply(colnames(count_table), function(column) {
  sub(".*\\.([^_]+_\\d+_S\\d+_L\\d+)_cut_PE1\\.fastq\\.gz.*", "\\1", column)
})
colnames(count_table) <- new_names

# Rename columns based on the mapping
colnames(count_table) <- sapply(colnames(count_table), function(old_name) {
  if (old_name %in% names(mapping)) {
    return(mapping[old_name])
  } else {
    return(old_name)
  }
})

# Sort columns
count_table <- count_table %>%
  select(sort(colnames(count_table)))

#### Para comparação par a par, selecionando apenas o inoculado e controle para a comparação
comparison <- c("IACSP_5503_sandy_inoc", "IACSP_5503_sandy_control")
# Select columns containing the base names
count_table <- count_table %>%
  select(contains(comparison))


# Repeat each base name three times for te strain object
strain <- c("inoc", "inoc", "inoc", "control", "control", "control")
sample_info <- data.frame(condition = factor(strain))
rownames(sample_info) <- colnames(count_table)
group <- factor(paste(strain))


### Follow the EdgeR analysis
library("edgeR")

y <- DGEList(counts = count_table, group = group)
cpms <- cpm(y)

# chatGPT automatization expression filtering
# filter if any samples triplicate has more than three cpms in the row sum
filter_criteria <- function(row, step = 3) {
  num_cols <- length(row)
  indices <- seq(1, num_cols, by = step)
  
  any(sapply(indices, function(i) all(row[i:(i + step - 1)] >= 1)))
}

# Apply the function to check each gene.
keep <- sapply(1:nrow(cpms), function(i) {
  filter_criteria(cpms[i, ])
})

#keep <- rowSums(cpms > 1) >= 3 ### Filter to be considered expressed.
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# MDS is only based on the top 500 genes... maybe using more it can look really different
plotMDS(y, col = c("blue", "blue", "blue",
                   "red", "red", "red"))

### Experimental Design ###
design <- model.matrix(~ 0 + group, data = y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y, design, robust = TRUE)
plotBCV(y)
y$design

### Export the table with normalized counts (cpm_tmm)
filename <- "IACSP_5503_sandy_inoc-control.csv"
cpms_with_gene_id <- cbind(gene_id = rownames(cpm(y)), cpm(y))
write.table(file=filename, cpms_with_gene_id, sep='\t', row.names=FALSE)

### Até aqui, temos apenas os count normalizados (cpm) sem atribuição de diferencialmente expressos


### Aplicar teste para a expressão diferencial
### Likelihood Ratio Test (LRT) for Differential Expression Testing ###
fit <- glmFit(y, design, robust = TRUE)
# We can change the contrast to be tested
my.contrast <- makeContrasts("inoc - control", levels = design)
lrt <- glmLRT(fit, contrast = my.contrast)
topTags(lrt)
#### With FDR
DE <- decideTests(lrt, adjust.method = "BH", p.value = 0.05, lfc = 1)
summary(DE)

# ### Without FDR
# DE <- decideTestsDGE(lrt, adjust.method = "none", p.value = 0.05)
# summary(DE)
# #plotSmear(lrt)

DE_tags <- rownames(y)[as.logical(DE)]
#MAplot
plotSmear(lrt, de.tags=DE_tags, cex = 0.3)
abline(h=c(-1, 1), col="blue")

### Export table with genes considered expressed (undergo filtering by minimum CPM in a certain number of samples)

### there must be a better way to export the tables!

# Apply the topTags() function to obtain the results.
top_tags <- topTags(lrt, n = (nrow(y$counts) + 1))
# Add the gene IDs as the first column.
results_with_gene_id <- cbind(gene_id = rownames(top_tags$table), top_tags$table)
# Export the results to a CSV file.
filename <- "Expressed_IACSP-5503_inoc-control.csv"
write.table(results_with_gene_id, file = filename, sep = '\t', row.names = FALSE) #maybe just set to TRUE, and don't need to cbind the gene_names


## Export the table with the DEGs (Differentially Expressed Genes).
threshold <- 0.05
### column 6:FDR, column 5:Pvalue
DEGs <- results_with_gene_id[results_with_gene_id[, 6] < threshold, ]
# Set a cutoff value for the absolute LogFC (Logarithm of the Fold Change).
logFC_threshold <- 1
# Calculate the absolute value of the LogFC (Logarithm of the Fold Change).
DEGs$absLogFC <- abs(DEGs[, 2])
# Filter the DEGs based on the absolute value of the LogFC (Logarithm of the Fold Change).
DEGs_LogFC <- DEGs[DEGs$absLogFC >= logFC_threshold, ]
# Remove the temporary column 'absLogFC'.
DEGs_LogFC <- DEGs_LogFC[, -ncol(DEGs_LogFC)]

filename <- "DEGs_IACSP-5503_inoc-control.csv"
write.table(DEGs_LogFC, file = filename, row.names=FALSE, sep = "\t")


```

## Kraken2 - Contaminats Analysis: Taxonomic Annotation (This part is not been using anymore)
Kraken2 informations: https://github.com/DerrickWood/kraken2/wiki/Manual

### 1 Step. Standard Kraken 2 Database

```
# Definir o caminho e nome do banco de dados
DBNAME="/home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data/Software/Kraken2/Kraken2_database"

# Criar o diretório do banco de dados
mkdir -p "$DBNAME"

# Baixar e construir o banco de dados padrão do Kraken2
kraken2-build --standard --db "$DBNAME" --threads 15

```
### 2 Step. Kraken2 Classification

```
#!/bin/bash

# Diretório de trabalho e do banco de dados
work_dir="/home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data"
database_dir="/home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data/Software/Kraken2/Kraken2_database"
input_dir="/home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data/Cutadapt"

# Mudança para o diretório de trabalho
cd "$work_dir"

# Criação de diretórios para saída
mkdir -p kraken2
mkdir -p kraken2/classified
mkdir -p kraken2/unclassified

# Loop para processar arquivos trimados
for input in "$input_dir"/*PE1.fastq.gz; do
    input2=$(echo "$input" | sed "s/PE1/PE2/g")
    output1=$(basename "$input" | sed 's/_PE1.fastq.gz//g')

    kraken2 --use-names --threads 20 --gzip-compressed \
        --db "$database_dir" \
        --report ./kraken2/"${output1}"_report.txt \
        --paired "$input" "$input2" \
        --report-zero-counts \
        --use-mpa-style \
        --classified-out ./kraken2/classified/"${output1}"_contamination#.fastq \
        --unclassified-out ./kraken2/unclassified/"${output1}"_filtered#.fastq \
        > ./kraken2/"${output1}"_kraken.txt

    # Compacta os arquivos classificados e não classificados
    gzip ./kraken2/classified/"${output1}"_contamination_1.fastq
    gzip ./kraken2/classified/"${output1}"_contamination_2.fastq
    gzip ./kraken2/unclassified/"${output1}"_filtered_1.fastq
    gzip ./kraken2/unclassified/"${output1}"_filtered_2.fastq
done

```
