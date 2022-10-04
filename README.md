# APP : Alienness-by-Phyletic-Pattern

APP, is a phyletic pattern analysis tool that uses comparative genomics to detect horizontally acquired genes. 

## Methodology

[APP architecture](https://github.com/sohamsg90/APP-Alienness-by-Phyletic-Pattern/blob/main/APP_architecture.png)
![APP_architecture](https://user-images.githubusercontent.com/43525619/187593303-09026817-1875-4e5e-acd8-449185bb14f3.png)

## Download files

```
git clone https://github.com/sohamsg90/APP-Alienness-by-Phyletic-Pattern
```

## Using the APP docker image
Pull the docker image

```
cd APP-Alienness-by-Phyletic-Pattern/example
docker pull sohamsg90/image_app_v1
```
## Input
We are using a random genome *Yersinia pestis KIM10+* to perform our sample run . The NCBI accession number is  [NC_004088](https://www.ncbi.nlm.nih.gov/nuccore/NC_004088.1/).

APP requires two input files - 
1. A file with taxonomic IDs of the associated genome
2. A multiFASTA protein sequence file

### 1. Constructing taxonomic ID file

Run the docker image and use the [input_file_prepare.pl](https://github.com/sohamsg90/APP-Alienness-by-Phyletic-Pattern/blob/main/scripts/taxonomy/input_file_prepare.pl) to create a taxonomic ID file `3.query_speciesID_taxID.txt`. 

An input file with the genome accession number (NCBI), [accession_number.txt](https://github.com/sohamsg90/APP-Alienness-by-Phyletic-Pattern/blob/main/example/accession_number.txt)  is to be provided.

```
docker run --rm -v "$(pwd)":/dir -w /dir sohamsg90/image_app_v1 /usr/local/bin/input_file_prepare.pl accession_number.txt
```
This generates a file `3.query_speciesID_taxID.txt` with the species ID and taxonomic ID of the genome (from NCBI taxonomy database).

The output of the file looks like this:
| Acc       | speciesid | taxid  |
| --------- | --------- | ------ |
| NC_004088 | 632       | 187410 |

This file has to be imported into the program. A new file is to be generated whenever analyzing a new genome.

### 2. MultiFASTA protein sequence file

APP has the capability to take input either a multiFASTA file or download the entire proteome of the genome of choice. 

**Option 1:**
In the former case, a file must be provided to the program specifying the type of input.

```
 docker run --rm -v "$(pwd)":/dir -w /dir sohamsg90/image_app_v1 /usr/local/bin/APP.pl -q NC_004088.faa -t 3.query_speciesID_taxID.txt -o NC_004088 -f multifasta
 ```

**Option 2:**
Alternatively, the user can supply the genome accession number (e.g. NC_004088) as a command-line input.

```
docker run --rm -v "$(pwd)":/dir -w /dir sohamsg90/image_app_v1 /usr/local/bin/APP.pl NC_004088 -t 3.query_speciesID_taxID.txt -o NC_004088 -f accession
```

Note that the input type has changed to `accession`. With the help of pre-installed NCBI eutilities, the program  will download the complete proteome accordingly.

**Option 3:**
If marker gene enrichment is to be performed along with generation of circular genome maps.

```
 docker run --rm -v "$(pwd)":/dir -w /dir sohamsg90/image_app_v1 /usr/local/bin/APP.pl -q NC_004088.faa -t 3.query_speciesID_taxID.txt -o NC_004088 -f multifasta -m 1 -g 1
 ```

***Tips:***

**Option 1** is directed when the user has a limited set of sequences to analyze, belonging to the same genome.

**Option 2** is when the user wants to perform phyletic pattern analysis on an entire genome (proteome) of interest.

## Execution
To run the program, simply type in:

```
docker run --rm -v "$(pwd)":/dir -w /dir sohamsg90/image_app_v1 /usr/local/bin/APP.pl -q <query fileName> -t <query taxonomy file> -o <Output fileName> -f <fileType> [Options]
```
***Note:*** The main script file, query sequence file and taxonomic ID file must be placed in the same working/current directory.

### Options

Following set of options are available with the program.

#### required arguments:
-q - multifasta amino acid file or genome acession number (NCBI only).

-t - query specific taxonomic id file.

-o - Output file name to create.

-f - File type as 'multifasta' or 'accession'.

#### optional arguments:

-n - No. of CPU cores to use for performing blast. By default, uses all available cores.

-m - Marker gene enrichment. Default set to 0. Use 1 only when whole genome is analyzed.

-g - Generate circular gene map with alien genes. Default set to 0. Use 1 only when whole genome is analyzed.

-e - turn on Expert option (1; keep temporary and intermediate files). Default set to 0. 

-v - Provide basic progress messages. Default set to 1.

-d - Provide detailed progress messages. Default set to 0.

## Expected space requirements

APP downloads genomes from NCBI ftp server in real-time. Henceforth, depending on the availability of completely sequenced genomes, and the abundance of genomes of the associated taxonomic ranks, the size of the downloaded database will vary, from genome to genome.


## Output analysis
In the default mode, the program provides a genome-wide list of genes (accession numbers) which are deemed horizontally acquired. Additionally, the program also predicts whether the horizontal gene transfer event was recent or ancient.


## Running in HPC environment
In high-performance cluster environemnt, user might not have sudo access to install docker daemon. 
HPC enviroments have the facilty to run singularity images without sudo access. 
### Pull docker image to build local singularity image
```
singularity pull APP.sif docker://sohamsg90/image_app_v1
```
### Execute using following command
```
singularity exec "$(pwd)"/APP.sif APP.pl -q "$(pwd)"/NC_004088.faa -t "$(pwd)"/3.query_speciesID_taxID.txt -o NC_004088 -f multifasta
```

