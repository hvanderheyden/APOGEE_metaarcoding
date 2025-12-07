#!/bin/bash
set -euo pipefail

# Enable error tracing and cleanup on exit
trap 'if [[ $? -ne 0 ]]; then echo "ERROR: Script failed at line $LINENO" >&2; fi; cleanup' EXIT
trap 'echo "INFO: Script interrupted" >&2; cleanup; exit 130' INT TERM

cleanup() {
  # Add cleanup logic here if needed (e.g., remove temp files)
  :
}

#############
## Usage ####
#############

#chmod +x /media/herve/10TB/Apogee/6_mock/6_minimap2/minimap2_v02.sh

# /media/herve/10TB/Apogee/6_mock/6_minimap2/minimap2_v02.sh \
# -i /media/herve/10TB/Apogee/6_mock/6_minimap2/reads \
# -o /media/herve/10TB/Apogee/6_mock/11_minimap2_clustered \
# -r /media/herve/10TB/Apogee/6_mock/6_minimap2/ITS-RefDB_V02.mmi \
# -s /media/herve/10TB/Apogee/6_mock/6_minimap2 \
# -t 32 \
# -c true \
# -x 0.98 \
# -T /media/herve/10TB/Apogee/6_mock/6_minimap2/taxonomy.tsv \
# -F /media/herve/10TB/Apogee/5_Scripts/filter_with_confidence.py \
# -C 1

#############
# ARGUMENTS #
#############

# Parse command-line arguments
while getopts ":i:o:r:s:t:c:x:w:T:F:C:" opt; do
  case ${opt} in
    i ) input_dir="$OPTARG" ;;
    o ) output_dir="$OPTARG" ;;
    r ) ref_db="$OPTARG" ;;
    s ) script_dir="$OPTARG" ;;
    t ) threads="$OPTARG" ;;
    c ) enable_clustering="$OPTARG" ;;
    x ) identity="$OPTARG" ;;
    w ) wordlength="$OPTARG" ;;
    T ) taxonomy_db="$OPTARG" ;;
    F ) filter_script="$OPTARG" ;;
    C ) confidence="$OPTARG" ;;
    \? ) echo "ERROR: Invalid option: -$OPTARG" 1>&2; exit 1 ;;
    : ) echo "ERROR: Invalid option: -$OPTARG requires an argument" 1>&2; exit 1 ;;
  esac
done
shift $((OPTIND -1))

########################
# Validate Arguments   #
########################

if [[ -z "${input_dir:-}" ]] || [[ -z "${output_dir:-}" ]] || [[ -z "${ref_db:-}" ]] || [[ -z "${script_dir:-}" ]] || [[ -z "${threads:-}" ]] || [[ -z "${taxonomy_db:-}" ]] || [[ -z "${filter_script:-}" ]]; then
  echo "ERROR: Missing required arguments" >&2
  echo "Usage: $0 -i <input_dir> -o <output_dir> -r <ref_db> -s <script_dir> -t <threads> -T <taxonomy_db> -F <filter_script> [-c <enable_clustering>] [-x <identity>] [-w <wordlength>] [-C <confidence>]" >&2
  exit 1
fi

# Validate that input directory and reference DB exist
if [[ ! -d "${input_dir}" ]]; then
  echo "ERROR: Input directory does not exist: ${input_dir}" >&2
  exit 1
fi

if [[ ! -f "${ref_db}" ]]; then
  echo "ERROR: Reference database does not exist: ${ref_db}" >&2
  exit 1
fi

if [[ ! -f "${taxonomy_db}" ]]; then
  echo "ERROR: Taxonomy database does not exist: ${taxonomy_db}" >&2
  exit 1
fi

if [[ ! -f "${filter_script}" ]]; then
  echo "ERROR: Filter script does not exist: ${filter_script}" >&2
  exit 1
fi

# Set defaults if not provided
identity="${identity:-0.97}"
wordlength="${wordlength:-10}"
enable_clustering="${enable_clustering:-false}"
confidence="${confidence:-1}"

# Validate numeric parameters
if ! [[ "${threads}" =~ ^[0-9]+$ ]] || [[ "${threads}" -lt 1 ]]; then
  echo "ERROR: threads must be a positive integer" >&2
  exit 1
fi

if ! [[ "${identity}" =~ ^0\.[0-9]+$ ]] || (( $(echo "${identity} < 0 || ${identity} > 1" | bc -l) )); then
  echo "ERROR: identity must be between 0 and 1" >&2
  exit 1
fi

echo "==============================================="
echo "Metabarcoding Pipeline v02"
echo "==============================================="
echo "Input dir:      ${input_dir}"
echo "Output dir:     ${output_dir}"
echo "Reference DB:   ${ref_db}"
echo "Threads:        ${threads}"
echo "Clustering:     ${enable_clustering}"
echo "Identity:       ${identity}"
echo "Word length:    ${wordlength}"
echo "Taxonomy DB:    ${taxonomy_db}"
echo "Filter script:  ${filter_script}"
echo "Confidence:     ${confidence}"
echo "==============================================="
echo ""

########################
# Create Output Directories #
########################

mkdir -p "${output_dir}/porechop"
mkdir -p "${output_dir}/nanofilt"
mkdir -p "${output_dir}/clusters"
mkdir -p "${output_dir}/chimera"
mkdir -p "${output_dir}/refilt"
mkdir -p "${output_dir}/mapped"
mkdir -p "${output_dir}/filteredPAFs"

echo "Output directories created."
echo ""

########################
# Read Pre-processing  #
########################

echo "Starting analysis"
echo "Pre-processing the reads"
echo "------------------------"
echo ""

###################
# Adapter Removal #
###################

if [[ -n "$(find "${output_dir}/porechop" -maxdepth 1 -type f -name '*.fastq' 2>/dev/null)" ]]; then
    echo "✓ Step 1: Skipping adapter removal (output already exists)"
else
    echo "► Step 1: Adapter removal"
    echo "  Processing files from: ${input_dir}"
    input_count=$(find "${input_dir}" -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fastq.gz' \) 2>/dev/null | wc -l)
    if [[ ${input_count} -eq 0 ]]; then
      echo "ERROR: No FASTQ files found in ${input_dir}" >&2
      exit 1
    fi
    echo "  Found ${input_count} input file(s)"
    
    time find "${input_dir}" -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fastq.gz' \) -print0 | \
      xargs -0 -I {} -P "${threads}" bash -c 'porechop -t 1 -i "$1" -o "${2}/porechop/$(basename "$1" | sed "s/\(\.fastq\(\.gz\)\?\)$/-porechop.fastq/")" || { echo "ERROR: Adapter removal failed for $1" >&2; exit 1; }' _ {} "${output_dir}" || exit 1
    echo "  ✓ Adapter removal completed"
fi
echo ""

#################################
# Length Filtering               #
#################################

if [[ -n "$(find "${output_dir}/nanofilt" -maxdepth 1 -type f -name '*.fastq' 2>/dev/null)" ]]; then
    echo "✓ Step 2: Skipping length filtering (output already exists)"
else
    echo "► Step 2: Length filtering"
    echo "  Filtering: quality >= 12, length: 350-1500 bp"
    
    time find "${output_dir}/porechop" -maxdepth 1 -type f -name '*.fastq' -print0 | \
      xargs -0 -I {} -P "${threads}" bash -c 'NanoFilt -q 12 -l 350 --maxlength 1500 < "$1" > "${2}/nanofilt/$(basename "$1" | sed "s/porechop/nanofilt/")" || { echo "ERROR: Length filtering failed for $1" >&2; exit 1; }' _ {} "${output_dir}" || exit 1
    echo "  ✓ Length filtering completed"
fi
echo ""

###################
# Clustering (Optional) with VSEARCH #
###################

if [[ "${enable_clustering}" == "true" ]] || [[ "${enable_clustering}" == "True" ]] || [[ "${enable_clustering}" == "TRUE" ]]; then
    if [[ -n "$(find "${output_dir}/clusters" -maxdepth 1 -type f \( -name '*.fasta' -o -name '*.fastq' \) 2>/dev/null)" ]]; then
        echo "✓ Step 3: Skipping clustering (output already exists)"
        clustering_input="${output_dir}/clusters"
    else
        echo "► Step 3: Clustering with VSEARCH"
        echo "  Identity threshold: ${identity}"
        echo "  Word length: ${wordlength}"
        
        time find "${output_dir}/nanofilt" -maxdepth 1 -type f -name '*.fastq' -print0 | \
          xargs -0 -I {} -P "${threads}" bash -c 'vsearch --cluster_fast "$1" --id '"${identity}"' --wordlength '"${wordlength}"' --centroids "${2}/clusters/$(basename "$1" .fastq)-centroids.fasta" 2>/dev/null || { echo "ERROR: Clustering failed for $1" >&2; exit 1; }' _ {} "${output_dir}" || exit 1
        echo "  ✓ Clustering completed"
        clustering_input="${output_dir}/clusters"
    fi
else
    echo "✓ Step 3: Clustering skipped (enable_clustering=false)"
    clustering_input="${output_dir}/nanofilt"
fi
echo ""

###################
# Chimera Removal  #
###################

if [[ -n "$(find "${output_dir}/chimera" -maxdepth 1 -type f -name '*.fasta' 2>/dev/null)" ]]; then
    echo "✓ Step 4: Skipping chimera removal (output already exists)"
else
    echo "► Step 4: Chimera removal"
    echo "  Using minimap2 ava-ont (all-vs-all alignment) + yacrd + scrubb"
    echo "  yacrd thresholds: min overlap coverage=4, min overlap len=0.4"
    
    time find "${clustering_input}" -maxdepth 1 -type f -name '*.fasta' -print0 | \
      xargs -0 -I {} bash -c '
        input_file="$1"
        output_dir="$2"
        threads="$3"
        basename_noext=$(basename "$input_file" .fasta)
        paf_file="/tmp/${basename_noext}.paf"
        yacrd_file="/tmp/${basename_noext}.yacrd"
        
        # Self-alignment for chimera detection
        minimap2 -x ava-ont -g 500 -t "${threads}" "$input_file" "$input_file" > "$paf_file" 2>/dev/null || { echo "ERROR: minimap2 failed for $input_file" >&2; exit 1; }
        
        # Mark chimeras
        yacrd -i "$paf_file" -o "$yacrd_file" -c 4 -n 0.4 scrubb -i "$input_file" -o "${output_dir}/chimera/${basename_noext}.scrubb.fasta" 2>/dev/null || { echo "ERROR: yacrd/scrubb failed for $input_file" >&2; exit 1; }
        
        # Cleanup temp files
        rm -f "$paf_file" "$yacrd_file"
      ' _ {} "${output_dir}" "${threads}" || exit 1
    
    echo "  ✓ Chimera removal completed"
fi
echo ""

####################
# Quality Filtering #
####################

#if [ "$(ls -A ${output_dir}/refilt)" ]; then
#    echo "Skipping quality filtering as ${output_dir}/refilt is already filled."
#else
#    echo ""
#    echo "5. Quality filtering"
#    echo "-------------------"
#    time ls ${output_dir}/chimera/*.scrubb.fastq | parallel "NanoFilt -q 12 < {} > ${output_dir}/refilt/{/.}-refilt.fastq" || { echo "Quality filtering failed"; exit 1; }
# fi

##########################
# Read Mapping (minimap2) #
##########################

if [[ -n "$(find "${output_dir}/mapped" -maxdepth 1 -type f -name '*.paf' 2>/dev/null)" ]]; then
    echo "✓ Step 5: Skipping mapping (output already exists)"
else
    echo "► Step 5: Read mapping with minimap2"
    echo "  Preset: map-ont (optimized for long reads ~800bp metabarcodes)"
    echo "  Reference: ${ref_db}"
    
    time find "${output_dir}/chimera" -maxdepth 1 -type f -name '*.fasta' -print0 | \
      xargs -0 -I {} bash -c 'minimap2 -x map-ont -Q -t "$3" --secondary=no -K 10M "$2" "$1" > "${4}/mapped/$(basename "$1" .fasta).paf" 2>/dev/null || { echo "ERROR: Mapping failed for $1" >&2; exit 1; }' _ {} "${ref_db}" "${threads}" "${output_dir}" || exit 1
    
    echo "  ✓ Mapping completed"
fi
echo ""

#############
# Filtering  #
#############

#############
# Filtering  #
#############

if [[ -f "${output_dir}/filteredPAFs/filtered_otu.tsv" ]]; then
    echo "✓ Step 6: Skipping PAF filtering with confidence (output already exists)"
else
    echo "► Step 6: Filtering PAF files with confidence"
    echo "  Script: ${filter_script}"
    
    if [[ ! -f "${filter_script}" ]]; then
      echo "ERROR: Filter script not found at ${filter_script}" >&2
      exit 1
    fi
    
    mkdir -p "${output_dir}/filteredPAFs"
    time "${filter_script}" -i "${output_dir}/mapped" -o "${output_dir}/filteredPAFs/filtered_otu.tsv" -c "${confidence}" || { echo "ERROR: PAF filtering with confidence failed" >&2; exit 1; }
    
    echo "  ✓ PAF filtering completed"
fi
echo ""

##################################
# Create OTU Table from Filtered PAFs #
##################################

if [[ -f "${output_dir}/otu_table.csv" ]]; then
    echo "✓ Step 7: Skipping OTU table creation (output already exists)"
else
    echo "► Step 7: Creating OTU table from filtered PAFs"
    echo "  Input: ${output_dir}/filteredPAFs/filtered_otu.tsv"
    
    if [[ ! -f "${output_dir}/filteredPAFs/filtered_otu.tsv" ]]; then
      echo "ERROR: Filtered OTU file not found at ${output_dir}/filteredPAFs/filtered_otu.tsv" >&2
      exit 1
    fi
    
    # Convert TSV to CSV format - keep all columns including confidence scores
    awk 'BEGIN {FS="\t"; OFS=","} {$1=$1; print}' "${output_dir}/filteredPAFs/filtered_otu.tsv" > "${output_dir}/otu_table.csv" 2>/dev/null || { echo "ERROR: OTU table conversion failed" >&2; exit 1; }
    
    echo "  ✓ OTU table created: ${output_dir}/otu_table.csv"
fi
echo ""

#####################################
# Create Phyloseq Taxonomy Table     #
#####################################

if [[ -f "${output_dir}/phyloseq_taxonomy.csv" ]]; then
    echo "✓ Step 8: Skipping taxonomy table creation (output already exists)"
else
    echo "► Step 8: Creating taxonomy table"
    echo "  Script: ${script_dir}/taxonomyTable.py"
    
    if [[ ! -f "${script_dir}/taxonomyTable.py" ]]; then
      echo "ERROR: taxonomyTable.py not found at ${script_dir}/taxonomyTable.py" >&2
      exit 1
    fi
    
    "${script_dir}/taxonomyTable.py" -i "${output_dir}/otu_table.csv" -t "${taxonomy_db}" > "${output_dir}/phyloseq_taxonomy.csv" 2>/dev/null || { echo "ERROR: Taxonomy table creation failed" >&2; exit 1; }
    
    echo "  ✓ Taxonomy table created: ${output_dir}/phyloseq_taxonomy.csv"
fi
echo ""
echo "==============================================="
echo "Pipeline completed successfully!"
echo "==============================================="
echo "Output files:"
echo "  OTU table:        ${output_dir}/otu_table.csv"
echo "  Taxonomy table:   ${output_dir}/phyloseq_taxonomy.csv"
echo "==============================================="
