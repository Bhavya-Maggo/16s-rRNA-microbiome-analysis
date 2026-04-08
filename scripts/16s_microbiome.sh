set -e
# CHANGE THIS PATH
BASE_PATH="your_path"
# Threads (adjust based on your CPU)
THREADS=8
# View qiime (.qzv) files here:
# https://view.qiime2.org/
# Metadata file
METADATA="$BASE_PATH/metadata.tsv"
################################################################
echo "Starting pipeline at $(date)"

# CREATE DIRECTORIES
mkdir -p $BASE_PATH/results/fastqc
mkdir -p $BASE_PATH/results/fastqc_trimmed
mkdir -p $BASE_PATH/results/trimmed
mkdir -p $BASE_PATH/results/multiqc
mkdir -p $BASE_PATH/results/multiqc_trimmed
mkdir -p $BASE_PATH/qiime
mkdir -p $BASE_PATH/exports
mkdir -p $BASE_PATH/fastq_raw
mkdir -p $BASE_PATH/R_results

# STEP 1: FASTQC ON RAW READS
echo "Running FastQC on raw reads..."
fastqc -t $THREADS $BASE_PATH/fastq_raw/*.fastq -o $BASE_PATH/results/fastqc
echo "Done..."

# STEP 2: MULTIQC SUMMARY
echo "Running MultiQC..."
multiqc $BASE_PATH/results/fastqc -o $BASE_PATH/results/multiqc -n multiqc_report
echo "Done..."

# STEP 3: TRIMMING WITH FASTP
echo "Running fastp trimming..."
for R1 in $BASE_PATH/fastq_raw/*_R1_001.fastq; do
  base=$(basename "$R1" _R1_001.fastq)
  R2="$BASE_PATH/fastq_raw/${base}_R2_001.fastq"
  out_R1="$BASE_PATH/results/trimmed/${base}_R1_trimmed.fastq"
  out_R2="$BASE_PATH/results/trimmed/${base}_R2_trimmed.fastq"
  report="$BASE_PATH/results/trimmed/${base}_fastp.html"
  echo "Processing $base"
  # -q 20 → quality threshold
  # -u 10 → % of low quality bases allowed
  # -l 30 → minimum length
  fastp -i "$R1" -I "$R2" \
    -o "$out_R1" -O "$out_R2" \
    -q 20 -u 10 -l 30 \
    -h "$report"
done
echo "Done..."

# STEP 4: FASTQC ON TRIMMED READS
echo "Running FastQC on trimmed reads..."
for fq in $BASE_PATH/results/trimmed/*_trimmed.fastq; do
  fastqc "$fq" -o $BASE_PATH/results/fastqc_trimmed
done
echo "Done..."

# STEP 5: MULTIQC SUMMARY
echo "Running MultiQC on Trimmed reads..."
multiqc $BASE_PATH/results/fastqc_trimmed -o $BASE_PATH/results/multiqc_trimmed -n multiqc_report
echo "Done..."

# STEP 6: ACTIVATE QIIME2 ENVIRONMENT
echo "Activating QIIME2..."
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate qiime2-2024.10

# STEP 7: CREATE MANIFEST FILE
echo "Creating manifest file..."
MANIFEST="$BASE_PATH/qiime/manifest.csv"
echo "sample-id,absolute-filepath,direction" > $MANIFEST
for f in $BASE_PATH/results/trimmed/*_R1_trimmed.fastq; do
  sample=$(basename $f | cut -d'_' -f1)
  r2=${f/_R1_/_R2_}
  echo "${sample},${f},forward" >> $MANIFEST
  echo "${sample},${r2},reverse" >> $MANIFEST
done
echo "Done..."

# STEP 8: IMPORT INTO QIIME2
echo "Importing data into QIIME2..."
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $MANIFEST \
  --output-path $BASE_PATH/qiime/demux.qza \
  --input-format PairedEndFastqManifestPhred33

qiime demux summarize \
  --i-data $BASE_PATH/qiime/demux.qza \
  --o-visualization $BASE_PATH/qiime/demux.qzv
echo "Done..."

################################################################
# IMPORTANT CHECKPOINT
# Open demux.qzv and VERIFY quality before proceeding.
# Where median quality falls below:
#   Q30 (good)
#   Q25 (acceptable)
#   <Q20 (bad → cut here)
# Ensure overlap of forward + reverse reads:
#   Forward trunc + Reverse trunc - Amplicon length >= 20-30 bp
################################################################

# STEP 9: DADA2 DENOISING
echo "Running DADA2..."
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $BASE_PATH/qiime/demux.qza \
  --p-trunc-len-f 245 \
  --p-trunc-len-r 238 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-n-threads $THREADS \
  --o-table $BASE_PATH/qiime/table.qza \
  --o-representative-sequences $BASE_PATH/qiime/rep-seqs.qza \
  --o-denoising-stats $BASE_PATH/qiime/stats.qza

qiime metadata tabulate \
  --m-input-file $BASE_PATH/qiime/stats.qza \
  --o-visualization $BASE_PATH/qiime/stats.qzv
echo "Done..."

################################################################
# IMPORTANT CHECKPOINT
# Open stats.qzv to check:
#   Reads retained after filtering: good rate is usually >70%
#   Reads merged successfully: >80% is ideal
################################################################

# Check representative sequences and feature table
qiime feature-table summarize \
  --i-table $BASE_PATH/qiime/table.qza \
  --o-visualization $BASE_PATH/qiime/table.qzv \
  --m-sample-metadata-file $BASE_PATH/metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data $BASE_PATH/qiime/rep-seqs.qza \
  --o-visualization $BASE_PATH/qiime/rep-seqs.qzv

# STEP 10: TAXONOMY CLASSIFICATION
echo "Classifying taxonomy..."
CLASSIFIER="$BASE_PATH/silva-138-99-nb-classifier.qza"
qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads $BASE_PATH/qiime/rep-seqs.qza \
  --p-n-jobs $THREADS \
  --o-classification $BASE_PATH/qiime/taxonomy.qza

qiime metadata tabulate \
  --m-input-file $BASE_PATH/qiime/taxonomy.qza \
  --o-visualization $BASE_PATH/qiime/taxonomy.qzv

qiime taxa barplot \
  --i-table $BASE_PATH/qiime/table.qza \
  --i-taxonomy $BASE_PATH/qiime/taxonomy.qza \
  --m-metadata-file $BASE_PATH/metadata.tsv \
  --o-visualization $BASE_PATH/qiime/taxa-bar-plots.qzv
echo "Done..."

# STEP 11: FILTER MITO/CHLORO FROM TABLE AND REP-SEQS
echo "Filtering mitochondria and chloroplast..."

qiime taxa filter-table \
  --i-table $BASE_PATH/qiime/table.qza \
  --i-taxonomy $BASE_PATH/qiime/taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table $BASE_PATH/qiime/table-clean.qza

qiime taxa filter-seqs \
  --i-sequences $BASE_PATH/qiime/rep-seqs.qza \
  --i-taxonomy $BASE_PATH/qiime/taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-sequences $BASE_PATH/qiime/rep-seqs-clean.qza
echo "Done..."

# STEP 12: PHYLOGENY
echo "Building phylogenetic tree from filtered sequences..."
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $BASE_PATH/qiime/rep-seqs-clean.qza \
  --o-alignment $BASE_PATH/qiime/aligned-clean.qza \
  --o-masked-alignment $BASE_PATH/qiime/masked-clean.qza \
  --o-tree $BASE_PATH/qiime/unrooted-clean.qza \
  --o-rooted-tree $BASE_PATH/qiime/rooted-clean.qza \
  --p-n-threads $THREADS
echo "Done..."


# Rarefaction curve BEFORE choosing sampling depth
echo "Rarefaction curve (choose sampling depth here)..."
qiime diversity alpha-rarefaction \
  --i-table $BASE_PATH/qiime/table-clean.qza \
  --i-phylogeny $BASE_PATH/qiime/rooted-clean.qza \
  --p-max-depth 8666 \
  --m-metadata-file $BASE_PATH/metadata.tsv \
  --o-visualization $BASE_PATH/qiime/alpha-rarefaction.qzv
echo "Done..."

# Sampling depth (set AFTER checking alpha-rarefaction.qzv)
SAMPLING_DEPTH=2500

# Filter out Mock samples
qiime feature-table filter-samples \
  --i-table $BASE_PATH/qiime/table-clean.qza \
  --m-metadata-file $BASE_PATH/metadata.tsv \
  --p-where "\"Description\" != 'Mock'" \
  --o-filtered-table $BASE_PATH/qiime/table-no-mock.qza

# STEP 13: DIVERSITY ANALYSIS
echo "Running diversity analysis..."
qiime diversity core-metrics-phylogenetic \
  --i-table $BASE_PATH/qiime/table-no-mock.qza \
  --i-phylogeny $BASE_PATH/qiime/rooted-clean.qza \
  --p-sampling-depth $SAMPLING_DEPTH \
  --m-metadata-file $METADATA \
  --output-dir $BASE_PATH/qiime/core-metrics \
  --p-n-jobs-or-threads $THREADS

# Alpha diversity significance test
qiime diversity alpha-group-significance \
  --i-alpha-diversity $BASE_PATH/qiime/core-metrics/shannon_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization $BASE_PATH/qiime/core-metrics/shannon-significance.qzv
  
qiime diversity alpha-group-significance \
  --i-alpha-diversity $BASE_PATH/qiime/core-metrics/faith_pd_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization $BASE_PATH/qiime/core-metrics/faith_pd_vector-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity $BASE_PATH/qiime/core-metrics/observed_features_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization $BASE_PATH/qiime/core-metrics/observed_features_vector-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity $BASE_PATH/qiime/core-metrics/evenness_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization $BASE_PATH/qiime/core-metrics/evenness_vector-significance.qzv

# Beta diversity significance test
# NOTE: --m-metadata-column must match a column in metadata.tsv
# Change "time" below if grouping column has a different name
qiime diversity beta-group-significance \
  --i-distance-matrix $BASE_PATH/qiime/core-metrics/bray_curtis_distance_matrix.qza \
  --m-metadata-file $METADATA \
  --m-metadata-column time \
  --o-visualization $BASE_PATH/qiime/core-metrics/bray-curtis-significance.qzv \
  --p-pairwise
  
# Weighted UniFrac significance
qiime diversity beta-group-significance \
  --i-distance-matrix $BASE_PATH/qiime/core-metrics/weighted_unifrac_distance_matrix.qza \
  --m-metadata-file $BASE_PATH/metadata.tsv \
  --m-metadata-column time \
  --o-visualization $BASE_PATH/qiime/core-metrics/weighted-unifrac-significance.qzv

# Unweighted UniFrac significance
qiime diversity beta-group-significance \
  --i-distance-matrix $BASE_PATH/qiime/core-metrics/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $BASE_PATH/metadata.tsv \
  --m-metadata-column time \
  --o-visualization $BASE_PATH/qiime/core-metrics/unweighted-unifrac-significance.qzv

# STEP 14: EXPORT FOR R
echo "Exporting results..."

# Feature table 
qiime tools export \
  --input $BASE_PATH/qiime/table-no-mock.qza \
  --output $BASE_PATH/exports

# Convert biom to TSV
biom convert \
  -i $BASE_PATH/exports/feature-table.biom \
  -o $BASE_PATH/exports/feature-table.tsv \
  --to-tsv

# Taxonomy
qiime tools export \
  --input-path $BASE_PATH/qiime/taxonomy.qza \
  --output-path $BASE_PATH/exports

# Filtered rep seqs (FASTA) 
qiime tools export \
  --input-path $BASE_PATH/qiime/rep-seqs-clean.qza \
  --output-path $BASE_PATH/exports

# Filtered rooted tree (Newick format for phyloseq)
qiime tools export \
  --input-path $BASE_PATH/qiime/rooted-clean.qza \
  --output-path $BASE_PATH/exports

################################################################
echo "Pipeline completed successfully at $(date)"