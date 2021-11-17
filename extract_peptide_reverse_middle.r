

library(ensembldb)
library(EnsDb.Hsapiens.v86)
#library(Gviz)

data = read.table("/Users/4472414/References/genome/Exon_Annotation_All_Surfaceome/exon_coordinates_for_R_20211114.txt", sep="\t", header=T);
for (i in (length(data[,1]) - 30000):1) {
	chr = data[i,2];
	coord = data[i,1];
	#chr = "3";
	#coord = "3:194357927-194361046"
	gnm <- GRanges(coord)
	edbx <- filter(EnsDb.Hsapiens.v86, filter = ~ seq_name == chr)
	gnm_tx <- genomeToTranscript(gnm, edbx)
	gnm_prt <- genomeToProtein(gnm, edbx)	
	
	if (start(gnm_prt[[1]]) != -1) {
		prt <- proteins(edbx, filter = ProteinIdFilter(names(gnm_prt[[1]])))
		seq = substr(prt$protein_sequence, start(gnm_prt[[1]]), end(gnm_prt[[1]]))
		output = toString(c(coord, seq, prt[[1]], prt[[2]]));
		write(output, file = "/Users/4472414/References/genome/Exon_Annotation_All_Surfaceome/exon_coordinates_all_protein_reverse_middle.txt", append = TRUE)
		
	} else {
		seq = "NA"
		output = toString(c(coord, seq));
		write(output, file = "/Users/4472414/References/genome/Exon_Annotation_All_Surfaceome/exon_coordinates_all_protein_reverse_middle.txt", append = TRUE)
	}
}

library(Gviz)

gnm <- GRanges(coord)
## Since we're using Ensembl chromosome names we have to set:
options(ucscChromosomeNames = FALSE)

## Define a genome axis track
gat <- GenomeAxisTrack(range = gnm)

## Get all genes in that region
gnm_gns <- getGeneRegionTrackForGviz(edbx, filter = GRangesFilter(gnm))
gtx <- GeneRegionTrack(gnm_gns, name = "tx", geneSymbol = TRUE,
                       showId = TRUE)

## Generate a higlight track
ht <- HighlightTrack(trackList = list(gat, gtx), range = gnm)
## plot the region
plotTracks(list(ht))
