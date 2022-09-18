#!/bin/bash

conda activate qiime2

# raw dataların olduğu ilk artifactı oluşturma
qiime tools import \
    --input-path manifest-1.tsv \ 
    --output-path reads.qza \
    --type 'SampleData[SequencesWithQuality]' \
    --input-format SingleEndFastqManifestPhred33V2
    
# veri özeti
qiime demux summarize \
    --i-data reads.qza \
    --o-visualization demux.qzv
    
###--p-trunc-len argümanundan emin olmadığım için 0 yaptım.
### Doğruluğundan şüphe etmekle birlikte internette okuduğum bir bilgiye göre amplicon verilerindeki gen bölgeleri korunmuş olduğu için onları tekrardan kesmek istemiyormuşuz.

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
    # gg-13-8-99-nb-classifier.qza dosyası şu linkten alınmıştır >> https://docs.qiime2.org/2022.8/data-resources/#taxonomy-classifiers-for-use-with-q2-feature-classifier

qiime feature-classifier classify-sklearn \ 
  --i-classifier gg-13-8-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
  
# taksonomi bilgisini excele aktar
qiime tools export \
  --input-file taxonomy.qza \
  --output-path phyloseq/
