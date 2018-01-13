seqtab.nochim<-readRDS('~/Dropbox/16S dada2 output/new/seqtab_final.rds')
write.csv(seqtab.nochim, file = "~/Dropbox/16S dada2 output/dada2_to_qiime/seqtab.nochim.csv")


export_taxa_table_and_seqs = function(seqtab.nochim, file_seqtab, file_seqs) {
  seqtab.t = as.data.frame(t(seqtab.nochim))
  seqs = row.names(seqtab.t)
  row.names(seqtab.t) = paste0("SV", 1:nrow(seqtab.t))
  outlist = list(data_loaded = seqtab.t)
  mctoolsr::export_taxa_table(outlist, file_seqtab)
  seqs = as.list(seqs)
  seqinr::write.fasta(seqs, row.names(seqtab.t), file_seqs)
}
export_taxa_table_and_seqs(seqtab.nochim,"~/Dropbox/16S dada2 output/dada2_to_qiime/seqtab.nochim.csv.otutable",
                           "~/Dropbox/16S dada2 output/dada2_to_qiime/seqtab.nochim.csv.repset")


#after running dada_to_qiime.sh return to R and create your phyloseq object

biom <- import_biom("~/Dropbox/16S OREI/16S dada2 output/dada2_to_qiime/otus.biom", parseFunction = parse_taxonomy_greengenes)
tree = read.tree("~/Dropbox/16S OREI/16S dada2 output/dada2_to_qiime/tree.tre")
map <- import_qiime_sample_data("~/Dropbox/16S OREI/16S dada2 output/dada2_to_qiime/16S_Metadata_final.txt")
from_qiime_16S <- merge_phyloseq(biom, map, tree)
sample_data(from_qiime_16S)
saveRDS(from_qiime_16S , file="~/Dropbox/16S OREI/16S dada2 output/dada2_to_qiime/from_qiime_16S.rds")
