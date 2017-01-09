# trim: A sliding window trimming tool for FASTQ format with Sanger encoded quality
## About
This tool is used to trim Next Generation Sequencing reads provided in FASTQ format from Illumina MiSeq in order to obtain a higher quality dataset that can be used for furhter downstream analyses. The reads are assumed to have a high quality base calls in the 5' end and a declining quality in the 3' end. Any deviation from this assumption might cause too many reads to be discarded.

Trim uses an algorithm that calculates an average quality value in a window by sliding over each quality record in a fastq file. Each read will be trimmed at the point in the 3' end where the smoothed average falls below a specified quality threshold value.

The script comes with a second option to run a cumulative sum algorithm for trimming. For each base in a read, the difference between a specified base error tolerance (default 0.05) and the base error probablility is added up to series of cumulative sums for each base position. A subsequence of the read is selected that is located between the first position that lies above the quality threshold in the 5' end and the position where the maximal cumulative sum value is found in the 3' end. 
