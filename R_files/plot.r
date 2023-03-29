library(Gviz)

# This is the bam file - the bam index file (.bai) must also be in the same dir
bam_file <- "R_files/SRR2273740_Aligned.sortedByCoord.out.bam"

# needed otherwise the chromosome names are not recognized:
options(ucscChromosomeNames = FALSE)

# This is the coverage - the blue histograms (type = "h")
coverage <- DataTrack(
    range = bam_file, type = "h",
    name = "Coverage", window = -1
    # optionally... transformation = e.g. function(x) { log2(x + 1) }
    # (can be any numeric function to transform the data)
)

# This is the axis for the chromosome (NW_017739545.1)
genome_track <- GenomeAxisTrack()

# This is the "transcripts" track with the arrows
transcripts_track <- AnnotationTrack(
    # I pulled these coordinates by manually reading off the GTF from stringtie
    # there will of course be a way to do this automatically
    start = c(608763, 611644, 610722, 611644),
    width = c(360, 257, 357, 257),
    chromosome = "NW_017739545.1",
    strand = rep("+", 4),
    group = rep(c(
        "XM_019714805.1",
        "XM_019714804.1"
    ), c(2, 2)),
    name = "Transcripts"
)

# do the plot
plotTracks(
    list( # order in this list is top to bottom on the graph
        transcripts_track,
        genome_track,
        coverage
    ),
    groupAnnotation = "group",
    # set the range of the plot:
    from = 605000,
    to = 614000
)
