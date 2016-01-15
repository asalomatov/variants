require("Rsamtools")
require("rbamtools")

### extracting information from bam files using rbamtools

bam <- bamReader("/mnt/scratch/asalomatov/data/SSC/wes/bam/14699.p1.realigned.recal.bam", idx=TRUE)
isOpen(bam)
indexInitialized(bam)
filename(bam)
align<-getNextAlign(bam)
name(align)
position(align)
str(getRefData(bam))
getRefCount(bam)
getRefName(bam)
getRefCoords(bam,"Y")
getRefCoords(bam,"Y")[1]
rg<-bamRange(bam,coords)
# for any numeric coordinates
coords <- as.integer(c(0, 802319, 802320))
rg <- bamRange(bam, coords)
size(rg)


rewind(rg)
align<-getNextAlign(rg)

#alignment field accessors
name(align)
flag(align)
refID(align)
position(align)
mapQuality(align)
cigarData(align)
nCigar(align)
mateRefID(align)
matePosition(align)
alignSeq(align)
alignQual(align)
alignQualVal(align)

#flag accessors
paired(align)
properPair(align)
unmapped(align)
mateUnmapped(align)
reverseStrand(align)
mateReverseStrand(align)
firstInPair(align)
secondInPair(align)
secondaryAlign(align)
failedQC(align)
pcrORopt_duplicate(align)

nucStats(rg)
#nucStats(bam)
coerce
pop_front
pop_back
push_front
push_back

# accsess given position in alignment
library(stringr)
n <- position(align)
myseq <- alignSeq(align)
mypos <- coords[2] - n + 1
substr(myseq, mypos, mypos)

ntAlign <- function(pos, align){
    #
    # align is extracted with getNextAlign
    # pos is a zero based position in reference genome
    #
    n <- position(align)
    mypos <- pos - n + 1
    substr(alignSeq(align), mypos, mypos)
}

ntAlign1 <- function(pos, align){
    #
    # align is extracted with getNextAlign
    # pos is a zero based position in reference genome
    #
    n <- position(align)
    mypos <- pos - n + 1
    strsplit(alignSeq(align), "")[[1]][mypos]
}

ntAlign2 <- function(pos, align){
    #
    # align is extracted with getNextAlign
    # pos is a zero based position in reference genome
    #
    n <- position(align)
    mypos <- pos - n + 1
    str_split(alignSeq(align), "")[[1]][mypos]
}

ntAlign(coords[2], align)
ntAlign1(coords[2], align)
ntAlign2(coords[2], align)
library(compiler)
ntAlign2 <- cmpfun(ntAlign1)
ntAlign4 <- cmpfun(ntAlign)
#align1 is the best
#install.packages("~/Downloads/microbenchmark_1.4-2.1.tar.gz")
require("microbenchmark")
microbenchmark(ntAlign(coords[2], align),ntAlign1(coords[2], align),
               ntAlign2(coords[2], align),ntAlign4(coords[2],align), times = 1000)

while(!is.null(align))
{
 print(reverseStrand(align))
 align<-getNextAlign(rg)
}

bamClose(bam)
#
# Open reader and initialize BAM-index
reader<-bamReader(bam,idx=TRUE)


### access pos in reference genome
#getSeq
#FaFile
#require("seqinr")
#read.fasta
hg19.fa <- "/mnt/xfs1/bioinfo/data/bcbio_nextgen/150607/genomes/Hsapiens/GRCh37/seq/GRCh37.fa" 
fa <- open(FaFile(hg19.fa))                   # open
# define range to extract
#GRanges
chrom <- '1'
pos <- 802320
k <- 0
gr <- GRanges(seqnames=Rle(c(chrom), c(1)),
            IRanges((pos-k):(pos+k), width=1), strand='*')
#(idx <- scanFaIndex(fa))
dna <- scanFa(fa, gr)
as.character(dna)
ntRef <- function(chr, pos, fa.file){
    gr <- GRanges(seqnames=Rle(c(chr), c(1)),
                IRanges(pos:pos, width=1), strand='*')
    as.character(scanFa(fa.file, gr))
}
#ntRef1 <- function(chr, pos, fa.file){
#    gr <- GRanges(seqnames=Rle(c(chr), c(1)),
#                IRanges(pos:pos, width=1), strand='*')
#    fa <- open(FaFile(fa.file))                   # open
#    as.character(scanFa(fa, gr))
#}
#ntRef("1", 802320, fa)
#ntRef1("1", 802320, hg19.fa)
#ntRef2 <- cmpfun(ntRef)
#microbenchmark(ntRef("1", 802320, fa), ntRef1("1", 802320, hg19.fa), ntRef2("1", 802320, fa), times = 1000)
# ntRef is best

# DNAshape
# library(devtools)
# install_github(repo = "TsuPeiChiu/DNAshapeR")
require(DNAshapeR)
browseVignettes("DNAshapeR")
getShape
getFasta
chrom <- 'chr1'
pos <- 802320
k <- 0
gr <- GRanges(seqnames=Rle(c(chrom), c(1)),
            IRanges((pos-k):(pos+k), width=1), strand='*')
library(BSgenome.Hsapiens.UCSC.hg19)
getFasta(gr, BSgenome = BSgenome.Hsapiens.UCSC.hg19, width = 100, filename = "tmp.fa")
fn <- "tmp.fa"
pred <- getShape(fn)
pred
str(pred)



dna[[1]]
DNAStringSet
ranges(idx) <- narrow(ranges(idx), -10)  # last 10 nucleotides
(dna <- scanFa(fa, param=idx[1:2]))



