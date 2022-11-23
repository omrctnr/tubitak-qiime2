# tubitak2209-project

### qiime --version 
    2cli version 2022.2.0

# Importing and Denoising
Hocam dada2’nin trunc-len argümanları için gönderdiğiniz kaynaklardan ve internetten faydalanarak bir çok kombinasyon denedim.
En son şu ayarları kullanarak analizi tamamladım `--p-trunc-len-f 258 --p-trunc-len-r 220 --p-trim-left-f 0 --p-trim-left-r 8`. Şu ana kadar elde edebildiğim en yüksek **Reads** ve **Feature** sayılarını  elde ettim. **(~410.000 okuma 1199 Features)** 

    qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path . --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path reads.qza
    
    qiime dada2 denoise-paired --i-demultiplexed-seqs reads_deleted.qza --p-trunc-len-f 258 --p-trunc-len-r 220 --p-trim-left-f 0 --p-trim-left-r 8  --p-n-reads-learn 50000 --o-representative-sequences rep-seqs.qza --o-table table.qza --o-denoising-stats stats.qza --p-n-threads 6

# Taxonomic Annotation
**GreenGenes Database**’yi kullanarak taksonomik atamayı gerçekleştirdim. 

    qiime feature-classifier classify-sklearn --i-classifier ../taxonomy_db/gg-13-8-99-nb-classifier.qza --i-reads rep-seqs.qza --o-classification tax.qza
    
# Visualization and Calculating Diversity Metrics 

    qiime taxa barplot --i-table table.qza --i-taxonomy tax.qza --m-metadata-file qiime2_metadata.tsv --o-visualization taxa-bar-plots.qzv

    qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza
    qiime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza
    qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
    qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza
    qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table.qza --p-sampling-depth 2000 --m-metadata-file qiime2_metadata.tsv --output-dir core-metrics-results
    qiime diversity alpha-rarefaction --i-table table.qza --i-phylogeny rooted-tree.qza --p-max-depth 2000 --m-metadata-file qiime2_metadata.tsv --o-visualization alpha-metrics-results/alpha-rarefaction.qzv
    
# Creating a table for use in R    
Bu kısımda bahsettiğiniz grafikleri ve çeşitlilik analizlerini R'da gerçekleştirebilmek için analiz verilerini dışarıya aktardım. Böylece R analizlerinde kullanılacak bir Excel dosyası oluşturdum. Size gönderdiğim **Figürlerden** kontrol edebilirsiniz.

    qiime tools export --input-path table.qza --output-path .
    biom convert -i feature-table.biom -o otu_table.tsv --to-tsv
    
    qiime metadata tabulate --m-input-file tax.qza --o-visualization tax.qzv
    qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file qiime2_metadata.tsv
    qiime feature-table transpose --i-table table.qza --o-transposed-feature-table transpozed-table.qza 
    qiime metadata tabulate --m-input-file rep-seqs.qza --m-input-file tax.qza --m-input-file transpozed-table.qza --o-visualization merged-data.qzv
    qiime tools export --input-path merged-data.qzv --output-path merged-data
    
