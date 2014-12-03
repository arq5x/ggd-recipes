mysql --user=genome --host=genome-mysql.cse.ucsc.edu \
   -A \
   -N \
   -B \
   -D hg19 \
   -e "SELECT chrom, chromStart, chromEnd, perGc, name from cpgIslandExt”
