# trim: A sliding window trimming tool for FASTQ files using phred quality scores
## About
This tool is used to trim Next Generation Sequencing reads provided in FASTQ format from Illumina MiSeq in order to obtain a higher quality dataset that can be used for furhter downstream analyses. The reads are assumed to have a high quality base calls in the 5' end and a declining quality in the 3' end. Any deviation from this assumption might cause too many reads to be discarded.

Trim uses an algorithm that slides over each quality record in a fastq file and calculates a "smooth" average. Each read will be trimmed at the point in the 3' end where the smoothed average falls below a specified quality threshold value. 
