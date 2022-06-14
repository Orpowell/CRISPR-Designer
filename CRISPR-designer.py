#!/usr/local/bin/python3

import re
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

# Codon Tab;e
codon_table = {
    'TCA': 'S',  # Serina
    'TCC': 'S',  # Serina
    'TCG': 'S',  # Serina
    'TCT': 'S',  # Serina
    'TTC': 'F',  # Fenilalanina
    'TTT': 'F',  # Fenilalanina
    'TTA': 'L',  # Leucina
    'TTG': 'L',  # Leucina
    'TAC': 'Y',  # Tirosina
    'TAT': 'Y',  # Tirosina
    'TAA': '*',  # Stop
    'TAG': '*',  # Stop
    'TGC': 'C',  # Cisteina
    'TGT': 'C',  # Cisteina
    'TGA': '*',  # Stop
    'TGG': 'W',  # Triptofano
    'CTA': 'L',  # Leucina
    'CTC': 'L',  # Leucina
    'CTG': 'L',  # Leucina
    'CTT': 'L',  # Leucina
    'CCA': 'P',  # Prolina
    'CCC': 'P',  # Prolina
    'CCG': 'P',  # Prolina
    'CCT': 'P',  # Prolina
    'CAC': 'H',  # Histidina
    'CAT': 'H',  # Histidina
    'CAA': 'Q',  # Glutamina
    'CAG': 'Q',  # Glutamina
    'CGA': 'R',  # Arginina
    'CGC': 'R',  # Arginina
    'CGG': 'R',  # Arginina
    'CGT': 'R',  # Arginina
    'ATA': 'I',  # Isoleucina
    'ATC': 'I',  # Isoleucina
    'ATT': 'I',  # Isoleucina
    'ATG': 'M',  # Methionina
    'ACA': 'T',  # Treonina
    'ACC': 'T',  # Treonina
    'ACG': 'T',  # Treonina
    'ACT': 'T',  # Treonina
    'AAC': 'N',  # Asparagina
    'AAT': 'N',  # Asparagina
    'AAA': 'K',  # Lisina
    'AAG': 'K',  # Lisina
    'AGC': 'S',  # Serina
    'AGT': 'S',  # Serina
    'AGA': 'R',  # Arginina
    'AGG': 'R',  # Arginina
    'GTA': 'V',  # Valina
    'GTC': 'V',  # Valina
    'GTG': 'V',  # Valina
    'GTT': 'V',  # Valina
    'GCA': 'A',  # Alanina
    'GCC': 'A',  # Alanina
    'GCG': 'A',  # Alanina
    'GCT': 'A',  # Alanina
    'GAC': 'D',  # Acido Aspartico
    'GAT': 'D',  # Acido Aspartico
    'GAA': 'E',  # Acido Glutamico
    'GAG': 'E',  # Acido Glutamico
    'GGA': 'G',  # Glicina
    'GGC': 'G',  # Glicina
    'GGG': 'G',  # Glicina
    'GGT': 'G'  # Glicina
}


# sgRNA designer class
class sgRNA:
    def __init__(self, nucleotide_sequence, amino_acid_position):
        self.sequence = nucleotide_sequence
        self.codon_position = (amino_acid_position * 3) - 3

        # Data outputs
        self._sgRNA_forward = None
        self._sgRNA_reverse = None
        self._60mer_switch = None

    # Make sgRNAs
    def make_sgRNAs(self):

        # Search for Thymine quadruplets in a given sequence and return boolean value
        def quadruple_nucleotide_check(nucleotide_sequence) -> bool:
            quadruplet = re.compile("([T])\\1\\1\\1")  # Regular expression for TTTT
            check = re.search(quadruplet, nucleotide_sequence)  # Search given sequence for RE

            if check is not None:  # If RE is present returns false
                return False
            else:  # If RE not found return true
                return True

        # Identify all PAM sites within range of target codon and return list of positions
        def find_pam_sites(xmer) -> list:
            # define search region for PAM sites within "xmer" nucleotides upstream of target codon
            search_region = self.sequence[self.codon_position: self.codon_position + xmer]
            i = 0
            pam_positions = []
            for _ in search_region:
                codon = search_region[i:i + 2]
                i += 1

                # If a PAM site is found store the position of the second base pair in PAM
                if codon == 'GG':
                    pam_position = i + 1  # Store position of last base of PAM site
                    pam_positions.append(pam_position + self.codon_position)  # Store value in a list

            return pam_positions

        # isolates 20mer from PAM sites and returns 20mer as a string
        def find_20mers(site, nucleotide_sequence) -> str:
            identified_20mer = nucleotide_sequence[0:site][-23:]  # Store sequence upstream of PAM site
            the_20mer = identified_20mer[:-3]  # Remove PAM site

            # Check for 4 consecutive thymine bases in 20mer
            quadruple = quadruple_nucleotide_check(str(identified_20mer))

            # If no quadruple bases found return start and end position and sequence of 20mer
            if quadruple:
                return the_20mer

        # isolates 60mer from PAM sites and returns list of 60mer sequence and PAM site (for synonymous mutation)
        def find_60mers(site, nucleotide_sequence) -> list:
            identified_60mer = nucleotide_sequence[0:site][-63:]  # Store sequence upstream of PAM site
            the_60mer = identified_60mer[:-3]  # Remove PAM site

            # Check for 4 consecutive thymine bases in 20mer
            quadruple = quadruple_nucleotide_check(str(identified_60mer))

            # If no quadruple bases found return sequence of 20mer and PAM site
            if quadruple:
                return [the_60mer, site]

        # Attempt to generate sgRNA using 20mer approach
        try:
            pam_sites = find_pam_sites(22)

            valid_20mers = [find_20mers(x, self.sequence) for x in pam_sites]  # Find all valid 20mer sequences

            result = [oligo for oligo in valid_20mers if oligo is not None]  # Remove None values

            oligo = result[-1]  # Choose last 20mer

            front = 'CGGGTGGCGAATGGGACTTT'  # Front primer
            back = 'GTTTTAGAGCTAGAAATAGC'  # back primer

            indentified_sgRNA = front + oligo.lower() + back  # generate forward sgRNA

            print('\n>PAM site found within 20 nucleotides...\n\n>Forward and Reverse sgRNA sequences generated...\n')
            self._sgRNA_forward = indentified_sgRNA
            self._sgRNA_reverse = str(Seq(indentified_sgRNA).reverse_complement())  # generate reverse sgRNA

        # Attempt to generate sgRNA using 60mer approach if 20mer fails
        except IndexError:
            try:
                print('\n>No PAM site within 20 nucletoides of target codon searching within 60 nucleotides... \n')

                pam_sites = find_pam_sites(62)

                valid_20mers = [find_60mers(x, self.sequence) for x in pam_sites]  # Find all valid 60mer sequences

                result = [oligo for oligo in valid_20mers if oligo[0] is not None]  # Remove None values

                oligo = result[-1][0][-20:]  # Target synonymous mutation site not actual mutation FIX ASAP

                front = 'CGGGTGGCGAATGGGACTTT'  # Front primer
                back = 'GTTTTAGAGCTAGAAATAGC'  # back primer

                indentified_sgRNA = front + oligo.lower() + back

                print('>PAM site found within 60 nucleotides\n')
                self._sgRNA_forward = indentified_sgRNA
                self._sgRNA_reverse = str(Seq(indentified_sgRNA).reverse_complement())
                self._60mer_switch = result[-1][1] // 3  # position of PAM site as amino acid

            # If 60mer method cannot design sgRNA then excit program
            except TypeError:
                print('error: No pam site found within 60 nucletoides of the target codon\n\n>Shutting down...\n')

                sys.exit(1)

    # get sgRNAs
    def get_sgRNA(self):
        return self._sgRNA_forward, self._sgRNA_reverse

    # get 60mer switch
    def get_switch(self):
        return self._60mer_switch


# Repair Template Designer class
class RepairTemplate:
    def __init__(self, nucleotide_sequence, amino_acid_position, amino_acid_mutation):
        # Initialise inputs
        self.sequence = nucleotide_sequence
        self.amino_acid_position = amino_acid_position
        self.codon_position = (self.amino_acid_position * 3) - 3
        self.amino_acid_mutation = amino_acid_mutation
        self._switch = None
        self.mutation = None

        # Data outputs
        self.repair_template_sequence = None
        self.forward_primer_sequence = None
        self.reverse_primer_sequence = None
        self.full_template_sequence = None

    # Make Repair template for 20mer
    def make_repair_template_for_20mer(self):
        # Convert sequence to list of codons
        codon_list = [self.sequence[i:i + 3] for i in range(0, len(self.sequence), 3)]

        # Save description of mutation in format intital aa, position, mutated aa e.g E140A
        self.mutation = Seq(codon_list[self.amino_acid_position - 1]).translate() + str(
            self.amino_acid_position) + self.amino_acid_mutation + '.txt'

        # Implement mutation at target site
        codon_list[self.amino_acid_position - 1] = (list(codon_table.keys())[
            list(codon_table.values()).index(self.amino_acid_mutation)]).lower()

        mutant_gene = "".join(codon_list)

        # generate 60 nt core template, two 70 nt primers and 160 nt full repair template
        repair_template_sequence = mutant_gene[
                                   self.codon_position - 30:self.codon_position + 30]
        forward_primer_sequence = mutant_gene[(self.codon_position - 80):(self.codon_position - 10)]
        reverse_primer_sequence = Seq(
            mutant_gene[(self.codon_position + 10):(self.codon_position + 80)]).reverse_complement()
        full_template = mutant_gene[(self.codon_position - 80):(self.codon_position + 80)]

        print('>Repair template and primers generated\n')
        self.repair_template_sequence = repair_template_sequence
        self.forward_primer_sequence = forward_primer_sequence
        self.reverse_primer_sequence = reverse_primer_sequence
        self.full_template_sequence = full_template

    # Make repair template for 60mer
    def make_repair_template_for_60mer(self):

        # Find a synonymous codon for synonymous mutation
        def synonymous_mutator(target, present_codon) -> string:
            synonymous_codons = [k for k, v in codon_table.items() if v == target and k != present_codon]
            return synonymous_codons[0]

        # Convert sequence to list of codons
        codon_list = [self.sequence[i:i + 3] for i in range(0, len(self.sequence), 3)]

        # Save description of mutation in format intital aa, position, mutated aa e.g E140A
        self.mutation = Seq(codon_list[self.amino_acid_position - 1]).translate() + str(
            self.amino_acid_position) + self.amino_acid_mutation + '.txt'

        # Implement mutation at target site
        codon_list[self.amino_acid_position - 1] = (list(codon_table.keys())[
            list(codon_table.values()).index(self.amino_acid_mutation)]).lower()

        # Implement synonymous mutation within 20 nucleotides of PAM site
        synonymous_mutation_site = self._switch - 5
        aa = Seq(codon_list[synonymous_mutation_site]).translate()
        codon_list[synonymous_mutation_site] = synonymous_mutator(aa, codon_list[synonymous_mutation_site]).lower()

        mutant_gene = "".join(codon_list)

        # generate 60 nt core template, two 70 nt primers and 160 nt full repair template
        core_start = self.codon_position - 12
        core_end = self.codon_position + 48

        repair_template_sequence = mutant_gene[core_start:core_end]
        forward_primer_sequence = mutant_gene[core_start - 50:core_start + 20]
        reverse_primer_sequence = Seq(mutant_gene[core_end - 20:core_end + 50]).reverse_complement()
        full_template = mutant_gene[(core_start - 50):(core_end + 50)]

        self.repair_template_sequence = repair_template_sequence
        self.forward_primer_sequence = forward_primer_sequence
        self.reverse_primer_sequence = reverse_primer_sequence
        self.full_template_sequence = full_template

    # Design 20mer or 60mer template based on presence of the 60mer switch
    def design_template(self):
        if type(self._switch) is int:
            self.make_repair_template_for_60mer()

        else:
            self.make_repair_template_for_20mer()

    # get templates
    def get_template(self):
        return self.forward_primer_sequence, self.reverse_primer_sequence, self.repair_template_sequence, self.full_template_sequence

    # set 60mer switch
    def set_switch(self, x):
        self._switch = x

    # get mutation description
    def get_mutation(self):
        return self.mutation


# Sequencing Primer designer class
class SequencingPrimer:
    def __init__(self, nucleotide_sequence, amino_acid_position):
        self.codon_position = (amino_acid_position * 3) - 3
        self.sequence = nucleotide_sequence
        self.primer = None

    # Design primer 50nt upstream of target codon
    def make_primer(self) -> str:
        primer = self.sequence[self.codon_position - 70:self.codon_position - 50]
        print(len(primer))
        self.primer = primer

    # get primer sequence
    def get_primer(self):
        return self.primer


# Output data (all sequences) class
class Output:

    def __init__(self, output_directory):
        # Data outputs
        self.sgRNA_forward = None
        self.sgRNA_reverse = None
        self.repair_template_sequence = None
        self.full_repair_template_sequence = None
        self.forward_primer_sequence = None
        self.reverse_primer_sequence = None
        self.sequencing_primer = None
        self.output_name = None
        self.output = output_directory

    # set sgRNA sequences from sgRNA class
    def set_sgRNA_sequences(self, forward, reverse):
        self.sgRNA_forward = forward
        self.sgRNA_reverse = reverse

    # set repair template sequences from RepairTemplate class
    def set_template_sequences(self, forward_primer, reverse_primer, core_sequence, full_template):
        self.reverse_primer_sequence = reverse_primer
        self.forward_primer_sequence = forward_primer
        self.repair_template_sequence = core_sequence
        self.full_repair_template_sequence = full_template

    def set_sequencing_primer(self, primer):
        self.sequencing_primer = primer

    # Set output file name
    def set_output_name(self, x):
        self.output_name = x

    # Generate output file with all data
    def create_output_file(self):

        print('>Writing output file...\n')

        # Saves output file in current working directory unless an alternate directory is given
        if self.output is None:
            file_path = self.output_name
        else:
            file_path = self.output + '/' + self.output_name

        with open(file_path, 'w+') as file:
            file.write(f' \n> sgRNA forward primer {len(self.sgRNA_forward)} bp\n{self.sgRNA_forward} \n')
            file.write(f' \n> sgRNA reverse primer {len(self.sgRNA_reverse)} bp\n{self.sgRNA_reverse}\n')
            file.write(f' \n> repair template {len(self.repair_template_sequence)}\n{self.repair_template_sequence}\n')
            file.write(
                f" \n> forward repair template primer {len(self.forward_primer_sequence)} bp\n{self.forward_primer_sequence}\n")
            file.write(
                f' \n> reverse repair template primer {len(self.reverse_primer_sequence)} bp\n{self.reverse_primer_sequence}\n')
            file.write(
                f' \n> full repair template {len(self.full_repair_template_sequence)} bp\n{self.full_repair_template_sequence}\n')
            file.write(f' \n> sequencing primer {len(self.sequencing_primer)} bp\n{self.sequencing_primer} \n')

        print(f'>Output file written to: {file_path}\n\n>Shutting down...\n')


# Wrapper class for mutant design
class MutantDesigner:
    def __init__(self, nucleotide_sequence, amino_acid_position, amino_acid_mutation, output):
        self.sgRNAs = sgRNA(nucleotide_sequence, amino_acid_position)
        self.seq_primer = SequencingPrimer(nucleotide_sequence, amino_acid_position)
        self.template = RepairTemplate(nucleotide_sequence, amino_acid_position, amino_acid_mutation)
        self.output_file = Output(output)

    # Designs and generates output file with all sequences
    def design(self):
        self.sgRNAs.make_sgRNAs()
        self.seq_primer.make_primer()
        self.template.set_switch(self.sgRNAs.get_switch())
        self.template.design_template()
        self.output_file.set_output_name(str(self.template.get_mutation()))
        self.output_file.set_sgRNA_sequences(*self.sgRNAs.get_sgRNA())
        self.output_file.set_template_sequences(*self.template.get_template())
        self.output_file.set_sequencing_primer(self.seq_primer.get_primer())
        self.output_file.create_output_file()


# Command Line Interface with Argparse
def cmd_lineparser():
    parser = argparse.ArgumentParser(prog='CRISPR-Designer', add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)

    group_inputs = parser.add_argument_group('Inputs')

    # Get path to fasta file containing nucleotide sequence
    group_inputs.add_argument('-s', '--sequence', metavar='\b', type=str, action='store',
                              help='path to fasta file',
                              default=None, required=True)

    # Get position of amino acid target
    group_inputs.add_argument('-p', '--position', metavar='\b', type=int, action='store',
                              help='position of target amino acid',
                              default=None, required=True)

    # Get amino acid to mutate to at target position
    group_inputs.add_argument('-m', '--mutant', metavar='\b', type=str, action='store',
                              help='mutation of target amino acid',
                              default=None, required=True)

    group_output = parser.add_argument_group('Outputs')
    # Get output directory (optional)
    group_output.add_argument('-o', '--output', metavar='\b', type=str, action='store',
                              help='directory to store output file', default=None)

    group_options = parser.add_argument_group('Options')
    # Get Version
    group_options.add_argument('-v', '--version', action='version', version='%(prog)s v2.1.0')
    # Get help
    group_options.add_argument("-h", "--help", action="help", help="show this help message and exit\n ")

    # Parse arguments
    arguments = parser.parse_args()
    input_list = [arguments.sequence, arguments.position, arguments.mutant, arguments.output]

    # If all arguments are None display help text
    if input_list.count(input_list[0]) == len(input_list):
        parser.parse_args(['-h'])

    # Check file provided is a fasta file
    if not (arguments.sequence.endswith('.fsa') or arguments.sequence.endswith('.fasta')):
        parser.error('--sequence requires .fsa or .fasta file as input')

    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*']

    # Check a valid single letter amino acid is provided
    if arguments.mutant not in amino_acids:
        parser.error('--mutant requires valid single amino acid code or "*" ')

    # Check an integer has been provided
    if type(arguments.position) is not int:
        parser.error('--position requires valid integer ')

    # check output is a real directory if provided
    if arguments.output is not None:
        if os.path.isdir(arguments.output) is False:
            parser.error('--output requires valid directory')

    return arguments


# Open fasta file from path provided
def open_fasta(path):
    # Parse fasta file into Biopython Seq object
    for dna_sequence in SeqIO.parse(path, "fasta"):
        return str(dna_sequence.seq)


if __name__ == '__main__':
    # Test interface
    args = cmd_lineparser()  # command line interface
    sequence = open_fasta(args.sequence)  # process fasta file
    mutant = MutantDesigner(sequence, args.position, args.mutant, args.output)
    mutant.design()
