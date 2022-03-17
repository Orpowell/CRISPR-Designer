# CRISPR Designer
CRISPR Designer is a tool for the rapid design of guide sequences and repair templates to generate point mutations within 60 nucleotides of a PAM site in a target protein.  The primers and sequences generated are designed for use with the protocol outlined by Ryan et. al 2016 titled, “Crispr–cas9 genome engineering in Saccharomyces cerevisiae cells”.

Dependencies for CRISPR designer can be found in [requirements.txt](https://github.com/Orpowell/CRISPR-Designer/blob/master/requirements.txt), Python >=3.7 is also required.

## Future Updates
Future updates may include:
- Codon replacement to facilitate easier design of repair templates
- Expansion of tool to include other methods outlined by Ryan et. al 2016 such as DNA bar coding (Please request if you’re interested!).

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
- The position of the nucleotide within the target mutation codon that you wish to mutate (1, 2 or 3).
- What nucleotide to replace the nucleotide at the chosen position in the target codon (A —> G)
- The name of output text file (The file is saved to the working directory)

If no appropriate PAM site can be found within 60 nucleotides region upstream of the target mutation, the programme will exit with no output. Although, It may still be possible to implement the mutation the on the reverse-strand or by using a larger search region upstream of nearby PAM sites. (Please refer to the original protocol by Ryan et. al, 2016)

CRISPR-designer is currently unable to assist in these scenarios. 
### Outputs
All outputs are written to a single text file formatted as shown below. The 20-mer guide RNAs and any single nucleotide in repair template and primer mutations are shown in lower case letters. 

	FILENAME.txt
	
	>Nucleotide Sequence
	sequence

	>Protein Sequence
	sequence

	> sgRNA forward primer
	sequence

	> sgRNA reverse primer
	sequence

	> repair template	
	sequence

	> forward repair template primer	
	sequence

	> reverse repair template primer
	sequence

	> full repair template
	sequence

## How CRISPR Designer works
The tool function in the following broad steps:

1. Ask user for the target protein nucleotide sequence in FASTA format
2. Read file and Identify all PAM sites in the sequence and the 20 nucleotide region upstream.
3. Ask user what mutation they would like to make
4. Filter 20mers sites for those containing the target mutation site
5. Generate sgRNA sequences and repair template with primers and save as a text file
6. If no 20mers are found repeat steps 2-5 using a 60 nucleotide upstream search region instead. 
7. If no appropriate PAM sites are found within 60 nucleotides the programme will exit with no output.

## Feedback
Any and all feedback is welcome, just raise an issue and I'll get back to you!

## Bibliography
<div class="csl-entry">Ryan, O. W., Poddar, S., &#38; Cate, J. H. D. (2016). Crispr–cas9 genome engineering in Saccharomyces cerevisiae cells. <i>Cold Spring Harbor Protocols</i>, <i>2016</i>(6), 525–533. https://doi.org/10.1101/pdb.prot086827</div>



