
# Application specification: SubspeciesClassification

This is the application specification for service with identifier SubspeciesClassification.

The backend script implementing the application is [App-SubspeciesClassification.pl](../service-scripts/App-SubspeciesClassification.pl).

The raw JSON file for this specification is [SubspeciesClassification.json](SubspeciesClassification.json).

This service performs the following task:   Subspecies Classification

It takes the following parameters:

| id | label | type | required | default value |
| -- | ----- | ---- | :------: | ------------ |
| input_source | Source of input (id_list, fasta_data, fasta_file, genome_group) | enum  | :heavy_check_mark: |  |
| input_fasta_data | Input sequence in fasta formats | string  |  |  |
| input_fasta_file | Input sequence as a workspace file of fasta data | wsid  |  |  |
| input_genome_group | Input sequence as a workspace genome group | wsid  |  |  |
| ref_msa_fasta | Reference multiple sequence alignment (Fasta-formatted) | wsid  |  |  |
| virus_type | Virus type | string  | :heavy_check_mark: |  |
| output_path | Output Folder | folder  | :heavy_check_mark: |  |
| output_file | File Basename | wsid  | :heavy_check_mark: |  |

