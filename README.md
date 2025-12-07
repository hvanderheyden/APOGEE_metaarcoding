
### DESCRIPTION
```
```
### Dependencies
```
Required software (install via conda):
  - porechop: Adapter removal
  - NanoFilt: Length/quality filtering
  - vsearch: Sequence clustering
  - minimap2: Read mapping and alignment
  - yacrd: Chimera detection
  - Python 3 with pandas: Data processing
```
## Installation
```
Install all dependencies:
  conda env create -f APOGEE.yml
  conda activate Apogee-pipeline

chmod +x */APOGEE.sh
```
### Basic usage
```
conda activate Apogee-pipeline

/minimap2_v02.sh \
  -i <input_dir> \
  -o <output_dir> \
  -r <ref_db> \
  -t <threads> \
  -T <taxonomy_db> \
  -F <filter_script> \
  -S <taxonomy_script> \
  [-c <enable_clustering>] \
  [-x <identity>] \
  [-w <wordlength>] \
  [-C <confidence>]
```
### REQUIRED ARGUMENTS
```
-i <input_dir>
     Path to directory containing input FASTQ files (gzipped or uncompressed)

  -o <output_dir>
     Path to directory for output files (will be created if it doesn't exist)

  -r <ref_db>
     Path to minimap2 index file (.mmi) for reference database
     Note: Create with: minimap2 -d output.mmi sequences.fasta

  -s <script_dir>
     Path to directory containing taxonomyTable.py script

  -t <threads>
     Number of threads for parallel processing (positive integer)

  -T <taxonomy_db>
     Path to taxonomy database file (TSV format: accession<tab>taxonomy)

  -F <filter_script>
     Path to filter_with_confidence.py script for PAF filtering
```
### OPTIONAL ARGUMENTS
```
  -c <enable_clustering>
     Enable sequence clustering with VSEARCH (true/false, default: false)

  -x <identity>
     Identity threshold for clustering (0.0-1.0, default: 0.97)

  -w <wordlength>
     Word length for clustering k-mer matching (default: 10)

  -C <confidence>
     Minimum confidence threshold for filtering (default: 1)

```
### OUTPUT FILES
```
  otu_table.csv
    Format: Accession,sample1,sample2,...,TotalCount,MappingConfidence
    Contains: OTU IDs, read counts per sample, total counts, confidence scores

  taxonomy_table.csv
    Format: #OTU ID,Domain,Phylum,Class,Order,Family,Genus,Species
    Contains: OTU IDs and full taxonomic classification only
```







