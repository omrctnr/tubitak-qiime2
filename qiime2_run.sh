#!/bin/bash

## komutlar sonucu elde edilen ana dosyalar ve excel çıktıları aşağıdaki linke eklenmiştir ##
    ## https://drive.google.com/drive/folders/19GzNsgi4OkEfDIX9mbbCOGixFwgE_fPN?usp=sharing ##

conda activate qiime2

# qiime --version >> q2cli version 2022.2.0

# raw dataların olduğu ilk artifactı oluşturma
qiime tools import \
    --input-path manifest-1.tsv \ 
    --output-path reads.qza \
    --type 'SampleData[SequencesWithQuality] \
    --input-format SingleEndFastqManifestPhred33V2
    
# veri özeti
qiime demux summarize \
    --i-data reads.qza \
    --o-visualization demux.qzv
    
# denoising/clustering adımı
qiime dada2 denoise-single \ 
    --i-demultiplexed-seqs reads.qza \
    --p-trunc-len 0 \ 
    --o-representative-sequences rep-seqs.qza \
    --o-table table.qza \
    --o-denoising-stats stats.qza
    
# dada2 sonuçlarını export et  
mkdir dada2_results
qiime tools export \
   --input-path table.qza \
   --output-path dada2_results
   
# biom dosyasını tsv'ye çevir
biom convert \
   -i dada2_results/feature-table.biom \
   -o dada/otu_table.tsv \
   --to-tsv

sed -i '1d' dada2_results/otu_table.tsv
sed -i 's/#OTU ID//' dada2_results/otu_table.tsv

# representative sequenceleri export et
qiime tools export \
  --input-path rep-seqs.qza \
  --output-path dada2_results

# taksonomik sınıflandırma
qiime feature-classifier classify-sklearn \ 
  --i-classifier gg-13-8-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
  
# taksonomi bilgisini excele aktar
qiime tools export \
  --input-file taxonomy.qza \
  --output-path phyloseq/
