MDS/IES - Annotation Sequence Software (MI-ASS)

Main Algorithm:

For each MAC contig C do:
	- Mask all telomeric sequences in C
	- BLAST C against all MIC contigs to obtain high scoring pairs list HSP (the word size is at least 28 bp)
		* Line up all high scoring pairs in HSP list against C and check for overlaps between different high scoring pairs
		* If any two high scoring pairs overlap at least 50%, merge them
		* Continue merging high scoring pairs until no more can be merged
		* The modified HSP list is considered to be initial list of MDSs
	- If there is no gaps in the obtained annotation, then done
	- Otherwise, BLAST C against all MIC contigs again to get shorter high scoring pairs (the word size is at least 12 bp)
	- Try to fill the gaps in the annotation with the shorter high scoring pairs
	- Output final MDS annotation for MAC contig C


Additional software and libraries used:
Basic Local Alignment Search Tool (BLAST) of version 2.2.29+ or higher
The download page is:
http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

Python libraries used:
1) Pyfasta
2) Numpy (as a part of Pyfasta)
3) Regex 

Make sure that the above libraries are installed (ex. with pip) before running MI-ASS

The main executable file is AnnotateCiliate.py.
EssentialFunctions.py contains all functions that are required during the execution of MI-ASS.
settings.py consists of the parameters that can be changed before running the program. Make sure that 5' and 3' telomeric regular expressions corresponding to your organism are set properly.

To run the program you need to have MAC sequence file, MIC sequence file (both are usually fasta files), and output direction specified.
The flag for specifying MAC sequence file is -mac, or --mac
The flag for specifying MIC sequence file is -mic, or --mic
The flag for specifying the output directory is -o, or --o
Additionaly, there is a flag to reblast previously blasted sequences: -reblast, or --rb

Example usage:

python AnnotateCiliate.py -mac /path/to/file/mac.fa -mic /path/to/file/mic.fa -o /path/to/output

Disable Debugging before running the program.