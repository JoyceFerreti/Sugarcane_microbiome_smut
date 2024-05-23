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


### Joyce Dellavechia Ferreti - Start at may/2024
# Pipeline 
![Untitled Diagram drawio (2)](https://github.com/JoyceFerreti/Sugarcane_microbiome_smut/assets/163196002/36cb5efb-1cea-49fd-8aa4-5f76edfaf2f0)



## :dna: Transcript alignment:


The analyses were done in the following order:
1. Get cutadapt reads and map them to R570_poliploid_version_2024/SofficinarumxspontaneumR570_771_v2.0.fa. reference genome located in  :open_file_folder: ```joycef@nioo0003.nioo.int:/home/nioo/joycef/Sugarcane_microbiome_rnaSeq/Data/References/R570_poliploid_version_2024```

2. Get mapped and unmapped reads to R570_poliploid_version_2024 genome and salve as
mapped in :open_file_folder: ```../Data/Alignments_sugarcane_microbiome_R570_genome_ShSHN/R570_mapped_filtered```
unmapped in :open_file_folder: ```../Data/Alignments_sugarcane_microbiome_R570_genome_ShSHN/R570_unmapped```

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
screen -S hisat2_alignment
```
Detach from the screen session: Press Ctrl + A followed by D.
Reattach to the screen session (if needed):
```
screen -r hisat2_alignment
```

### 1step: Reference Index

```
#!/bin/bash

if [ ! -d ../../Data/References/R570_poliploid_version_2024/hisat2_index ]
then
    mkdir -p ../../Data/References/R570_poliploid_version_2024/hisat2_index
    hisat2-build -p 10 ../../Data/References/R570_poliploid_version_2024/SofficinarumxspontaneumR570_771_v2.0.fa ../../Data/References/R570_poliploid_version_2024/hisat2_index/SofficinarumxspontaneumR570_771_v2.0
fi

```

### 2step: Reference Alignment

```
#!/bin/bash

# Verifica se o diretório de mapeamento já existe, se não, cria
if [ ! -d ../../Data/Alignments_sugarcane_microbiome_R570_genome_ShSHN ]; then
  mkdir -p ../../Data/Alignments_sugarcane_microbiome_R570_genome_ShSHN
fi

# Cria diretórios para leituras mapeadas e não mapeadas
mapped_dir="../../Data/Alignments_sugarcane_microbiome_R570_genome_ShSHN/R570_mapped"
unmapped_dir="../../Data/Alignments_sugarcane_microbiome_R570_genome_ShSHN/R570_unmapped"
mkdir -p "$mapped_dir"
mkdir -p "$unmapped_dir"

# Define o diretório de trabalho
cd ../../Data/Alignments_sugarcane_microbiome_R570_genome_ShSHN

# Caminho para o índice do HISAT2
index=../References/R570_poliploid_version_2024/hisat2_index/SofficinarumxspontaneumR570_771_v2.0

# Nome do arquivo de controle
checkpoint_File="./processed_files.txt"

# Verifica se o arquivo de controle existe, se não existir, cria
if [ ! -e "$checkpoint_File" ]; then
  touch "$checkpoint_File"
fi

# Função para lidar com o sinal de interrupção (Ctrl+C)
function cleanup {
  echo "Script interrupted. Exiting..."
  exit 1
}

# Captura o sinal de interrupção (Ctrl+C)
trap cleanup SIGINT

for PE1 in ../Cutadapt/*PE1.fastq.gz; do
  PE2="${PE1/PE1/PE2}"
  BASE_PE1=$(basename "$PE1")

  # Verifica se o par já foi processado antes de executar o HISAT2
  if ! grep -q "$BASE_PE1" "$checkpoint_File"; then
    echo "Analyzing $PE1 and $PE2"

    hisat2 -p 10 --rg-id "$PE1" --rg SM:"$PE1" \
    --summary-file ./summary_"$BASE_PE1"_cutadapt.txt \
    -x "$index" -1 "$PE1" -2 "$PE2" \
    -S ./"$BASE_PE1"_cutadapt.sam

    samtools sort -@ 10 ./"$BASE_PE1"_cutadapt.sam > \
    ./"$BASE_PE1"_cutadapt.bam
    # Devido ao tamanho de cromossomo maior que 512Mb, não é possível criar o bam index do tipo .bai; usamos então o parâmetro "-c" para criar o index do tipo .csi
    samtools index -@ 10 -c ./"$BASE_PE1"_cutadapt.bam
    rm ./"$BASE_PE1"_cutadapt.sam

    # Filtra as leituras mapeadas e as move para o diretório de leituras mapeadas
    samtools view --threads 10 -b -F 4 ./"$BASE_PE1"_cutadapt.bam > "$mapped_dir/${BASE_PE1}_mapped.bam"

    # Filtra as leituras não mapeadas e as move para o diretório de leituras não mapeadas
    samtools view --threads 10 -b -f 4 ./"$BASE_PE1"_cutadapt.bam > "$unmapped_dir/${BASE_PE1}_unmapped.bam"

    # Adiciona o nome do arquivo ao arquivo de controle
    echo "$BASE_PE1" >> "$checkpoint_File"
  else
    echo "Pair $BASE_PE1 has been processed previously. Ignoring."
  fi
done

```

Separation of Mapped and Unmapped Reads:
After the HISAT2 classification step, the script filters the mapped and unmapped reads using samtools view.
The mapped reads are moved to the directory R570_mapped_filtered, while the unmapped reads are moved to the directory R570_unmapped.

## Unmapped anotation - Kraken2
Kraken2 informations: https://github.com/DerrickWood/kraken/blob/master/docs/MANUAL.markdown

### 1 Step. 
