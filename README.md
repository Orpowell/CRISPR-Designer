# CRISPR Designer
CRISPR Designer is a tool for the rapid design of guide sequences and repair templates to generate point mutations within 60 nucleotides of a PAM site in a target protein.  The primers and sequences generated are designed for use with the protocol outlined by Ryan et. al 2016 titled, “Crispr–cas9 genome engineering in Saccharomyces cerevisiae cells”.

Dependencies for CRISPR designer can be found in [requirements.txt](https://github.com/Orpowell/CRISPR-Designer/requirements.txt), Python >=3.7 is also required.
## Installation on MacOSX and Linux
At the command line, change directory to the directory where CRISPR-designer.py was downloaded, E.g , using the full path name.

	cd <download-directory>

Now move the file to where you normally keep your binaries. This directory should be in your path. Note: you may require administrative privileges to do this (either switching user to root or by using sudo).

As root:

	mv CRISPR-designer.py /usr/local/bin/

As regular user:

	sudo mv CRISPR-designer.py /usr/local/bin/

After installation, CRISPR Designer can be run directly from the shell or Terminal using the following command:

	CRISPR-designer.py

Alternatively, CRISPR Designer can be run from an IDE.

## Using CRISPR Designer
![Tutorial](https://github.com/Orpowell/CRISPR-Designer/blob/master/tutorial.gif)

Before using this tool please familiarise yourself with the protocol (see bibliography)
### Inputs
CRISPR Designer requires the following inputs while running:
- The path to the FASTA file containing the nucleotide sequence of the target protein
- The position  of the amino acid to be mutated (e.g 174).
- The amino acid to subsititute at the position provided (single letter code only) (e.g A)
- (Optional) A directory to store the output file in. If no directory is provided the file will be saved to the current working directory

CRISPR designer will then automatically design sgRNA sequences, repair templates and primers based on the information provided by the user.  Codons are automatically mutated to implement the request mutation and in the case of 60mer protocol a synonymous mutation is also automatically introduced.

If no appropriate PAM site can be found within 60 nucleotides region upstream of the target mutation, the programme will exit with no output. Although, It may still be possible to implement the mutation the on the reverse-strand or by using a larger search region upstream of nearby PAM sites. However, CRISPR-designer is currently unable to assist in these scenarios. (Please refer to the original protocol by Ryan et. al, 2016) 

### Output
All data is written to a single output file. The file is named based on mutation generated by the user. E.g If you mutated glutamine to alanine at position 140 of your protein, the output file would be E140A.txt .

Data is formatted similiar to a fasta file with each header displaying the sequence name and length in base pairs (n bp). Mutations to the provided sequence are highlighted by lower case lettering. 

	> sgRNA forward primer n bp
	sequence

	> sgRNA reverse primer n bp
	sequence

	> repair template n bp
	sequence

	> forward repair template primer	 n bp
	sequence

	> reverse repair template primer n bp
	sequence

	> full repair template n bp
	sequence

## Future Updates
Future updates may include:
- Support for Windows 
- Expansion of tool to include other methods outlined by Ryan et. al 2016 such as DNA bar coding (Please request if you’re interested!).

## Feedback
Any and all feedback is welcome, just raise an issue and I'll get back to you!

## Bibliography
<div class="csl-entry">Ryan, O. W., Poddar, S., &#38; Cate, J. H. D. (2016). Crispr–cas9 genome engineering in Saccharomyces cerevisiae cells. <i>Cold Spring Harbor Protocols</i>, <i>2016</i>(6), 525–533. https://doi.org/10.1101/pdb.prot086827</div>



