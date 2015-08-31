MDS/IES - Annotation Sequence Software (MI-ASS)

Software used:
Basic Local Alignment Search Tool (BLAST) of version 2.2.29+ or higher
The download page is:
http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

Python libraries used:
1) Pyfasta
2) Numpy (as a part of Pyfasta)
3) Regex 

Make sure that the above libraries are installed (ex. with pip) before running MI-ASS

The main executable file is AnnotateCiliate.py
EssentialFunctions.py contains all functions that ar erequired during the execution of MI-ADD
settings.py consists of the parameters that can be changed before running the program. Make sure that 5' and 3' telomeric regular expressions corresponding to your organism are set.

To run the program you need to have MAC sequence file, MIC sequence file (both are usually fasta files), and output direction specified.
The flag for specifying MAC sequence file is -mac, or --mac
The flag for specifying MIC sequence file is -mic, or --mic
The flag for specifying the output directory is -o, or --o
Additionaly, there is a flag to reblast previously blasted seqeunces: -reblast, or --rb

Example usage:

python AnnotateCiliate.py -mac /path/to/file/mac.fa -mic /path/to/file/mic.fa -o /path/to/output

Disable Debugging before running the program.