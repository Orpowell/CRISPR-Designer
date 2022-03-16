# CRISPR Designer
CRISPR Designer is a tool for the rapid design of guide sequences and repair templates to generate point mutations within 60 nucleotides of a PAM site in a target protein.  The primers and sequences generated are designed for the protocol outlined in “CRISPR-Cas9 Genome Engineering in Saccharomyces cerevisiae Cells” (Ryan et. al 2016).

# Installation on MacOSX and Linux
At the command line, change directory to the directory where ava.py was downloaded, E.g , using the full path name.

	cd <download-directory>

Now move the file to where you normally keep your binaries. This directory should be in your path. Note: you may require administrative privileges to do this (either switching user to root or by using sudo).

As root:

	mv CRISPR-designer.py /usr/local/bin/

As regular user:

	sudo mv CRISPR-designer.py /usr/local/bin/

After installation, A.V.A can be run directly from the shell or Terminal using the follwing command:

	CRISPR-designer.py

Alternatively, A.V.A can be run from an IDE.

## Using CRISPR Designer

[Tutorial](https://github.com/Orpowell/CRISPR-Designer/blob/master/tutorial.gif)




