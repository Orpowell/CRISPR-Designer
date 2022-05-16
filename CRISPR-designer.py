#!/usr/local/bin/python3.7

import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq

codon_table = codontab = {
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


class Output:

    def __init__(self, output_path):
        # Data outputs
        self.sgRNA_forward = None
        self.sgRNA_reverse = None
        self.repair_template_sequence = None
        self.full_repair_template_sequence = None
        self.forward_primer_sequence = None
        self.reverse_primer_sequence = None
        self.output = output_path

    def set_sgRNA_sequences(self, forward, reverse):
        self.sgRNA_forward = forward
        self.sgRNA_reverse = reverse

    def set_template_sequences(self, forward_primer, reverse_primer, core_sequence, full_template):
        self.reverse_primer_sequence = reverse_primer
        self.forward_primer_sequence = forward_primer
        self.repair_template_sequence = core_sequence
        self.full_repair_template_sequence = full_template

    def create_output_file(self):
        with open(self.output, 'w+') as file:
            file.write(f' \n> sgRNA forward primer {len(self.sgRNA_forward)} bp\n{self.sgRNA_forward} \n')
            file.write(f' \n> sgRNA reverse primer {len(self.sgRNA_reverse)} bp\n{self.sgRNA_reverse}\n')
            file.write(f' \n> repair template {len(self.repair_template_sequence)}\n{self.repair_template_sequence}\n')
            file.write(
                f" \n> forward repair template primer {len(self.forward_primer_sequence)} bp\n{self.forward_primer_sequence}\n")
            file.write(
                f' \n> reverse repair template primer {len(self.reverse_primer_sequence)} bp\n{self.reverse_primer_sequence}\n')
            file.write(f' \n> full repair template {len(self.full_repair_template)} bp\n{self.full_repair_template}\n')

    def check_ouput_file(self):
        print(f' \n> sgRNA forward primer {len(self.sgRNA_forward)} bp\n{self.sgRNA_forward} \n')
        print(f' \n> sgRNA reverse primer {len(self.sgRNA_reverse)} bp\n{self.sgRNA_reverse}\n')
        print(f' \n> repair template {len(self.repair_template_sequence)}\n{self.repair_template_sequence}\n')
        print(
            f" \n> forward repair template primer {len(self.forward_primer_sequence)} bp\n{self.forward_primer_sequence}\n")
        print(
            f' \n> reverse repair template primer {len(self.reverse_primer_sequence)} bp\n{self.reverse_primer_sequence}\n')
        print(
            f' \n> full repair template {len(self.full_repair_template_sequence)} bp\n{self.full_repair_template_sequence}\n')


class sgRNA:
    def __init__(self, nucleotide_sequence, amino_acid_position):
        self.sequence = nucleotide_sequence
        self.__codon_position = (amino_acid_position * 3) - 3
        self.__sgRNA_forward = None
        self.__sgRNA_reverse = None

    # Make sgRNAs
    def make_sgRNAs(self):

        def quadruple_nucleotide_check(nucleotide_sequence):
            quadruplet = re.compile("([T])\\1\\1\\1")  # Regular expression for TTTT
            check = re.search(quadruplet, nucleotide_sequence)  # Search given sequence for RE

            if check is not None:  # If RE is present returns false
                return False
            else:  # If RE not found return true
                return True

        def find_20mers(site, nucleotide_sequence):
            identified_20mer = nucleotide_sequence[0:site][-23:]  # Store upstream sequence
            the_20mer = identified_20mer[:-3]  # Remove PAM site

            quadruple = quadruple_nucleotide_check(
                str(identified_20mer))  # Check for 4 consecutive thymine bases in 20mer

            # If no quadruple bases found return start and end position and sequence of 20mer
            if quadruple:
                return the_20mer

        # define search region for PAM sites within 22 nt upstream of codon
        search_region = self.sequence[self.__codon_position: self.__codon_position + 22]
        i = 0
        pam_sites = []
        for _ in search_region:
            codon = search_region[i:i + 2]
            i += 1

            # If a PAM site is found store the position of the second base pair in PAM
            if codon == 'GG':
                pam_position = i + 1  # Store position of last base of PAM site
                pam_sites.append(pam_position + self.__codon_position)  # Store value in a list

        test = [find_20mers(x, self.sequence) for x in pam_sites]  # Identify 20mers upstream of PAM sites

        result = [oligo for oligo in test if oligo is not None]  # Remove None values

        oligo = result[-1]

        front = 'CGGGTGGCGAATGGGACTTT'  # Front primer
        back = 'GTTTTAGAGCTAGAAATAGC'  # back primer

        indentified_sgRNA = front + oligo.lower() + back

        self.__sgRNA_forward = indentified_sgRNA
        self.__sgRNA_reverse = str(Seq(indentified_sgRNA).reverse_complement())

    # get sgRNAs
    def get_sgRNA(self):
        return self.__sgRNA_forward, self.__sgRNA_reverse


class RepairTemplate:
    def __init__(self, nucleotide_sequence, amino_acid_position, amino_acid_mutation):
        # Initialise inputs
        self.sequence = nucleotide_sequence
        self.amino_acid_position = amino_acid_position
        self.codon_position = (self.amino_acid_position * 3) - 3
        self.amino_acid_mutation = amino_acid_mutation

        # Data outputs
        self.repair_template_sequence = None
        self.forward_primer_sequence = None
        self.reverse_primer_sequence = None
        self.full_template_sequence = None

    # Make Repair template
    def make_repair_template(self):
        codon_list = [self.sequence[i:i + 3] for i in range(0, len(self.sequence), 3)]

        codon_list[self.amino_acid_position - 1] = (list(codon_table.keys())[
            list(codon_table.values()).index(self.amino_acid_mutation)]).lower()

        mutant_gene = "".join(codon_list)
        repair_template_sequence = mutant_gene[
                                   self.codon_position - 30:self.codon_position + 30]  # Create 60mer with 30 nt either side of mutation

        forward_primer_sequence = mutant_gene[(self.codon_position - 80):(self.codon_position - 30)] + mutant_gene[
                                                                                                       (
                                                                                                               self.codon_position - 30):(
                                                                                                               self.codon_position - 10)]
        reverse_primer_sequence = Seq(mutant_gene[(self.codon_position + 10):(self.codon_position + 30)] + mutant_gene[(
                                                                                                                               self.codon_position + 30):(
                                                                                                                               self.codon_position + 80)]).reverse_complement()
        full_template = mutant_gene[(self.codon_position - 80):(self.codon_position + 80)]

        self.repair_template_sequence = repair_template_sequence
        self.forward_primer_sequence = forward_primer_sequence
        self.reverse_primer_sequence = reverse_primer_sequence
        self.full_template_sequence = full_template

    def get_template(self):
        return self.forward_primer_sequence, self.reverse_primer_sequence, self.repair_template_sequence, self.full_template_sequence


if __name__ == '__main__':
    # S288C_YMR202W_ERG2_coding.fsa

    sequence = 'ATGAAGTTTTTCCCACTCCTTTTGTTGATTGGTGTTGTAGGCTACATTATGAACGTATTGTTCACTACCTGGTTGCCAACCAATTACATGTTCGATCCAAAAACTTTGAACGAAATATGTAACTCGGTGATTAGCAAACACAACGCAGCAGAAGGTTTATCCACTGAAGACCTGTTACAGGATGTCAGAGACGCACTTGCCTCTCATTACGGGGACGAATACATCAACAGGTACGTCAAAGAAGAATGGGTCTTCAACAATGCTGGTGGTGCGATGGGCCAAATGATCATCCTACACGCTTCCGTATCCGAGTACTTAATTCTATTCGGAACCGCTGTTGGTACTGAAGGGCACACAGGTGTTCACTTTGCTGACGACTATTTTACCATCTTACATGGTACGCAAATCGCAGCATTGCCATATGCCACTGAAGCCGAAGTTTACACTCCTGGTATGACTCATCACTTGAAGAAGGGATACGCCAAGCAATACAGCATGCCAGGTGGTTCCTTTGCCCTTGAATTGGCTCAAGGCTGGATTCCATGTATGTTGCCATTCGGGTTTTTGGACACTTTCTCCAGTACTCTTGATTTATACACTCTATATAGAACTGTCTACCTGACTGCCAGGGACATGGGTAAGAACTTGTTGCAAAACAAAAAGTTCTAA'

    # Make sgRNA
    sgRNAs = sgRNA(sequence, 174)
    sgRNAs.make_sgRNAs()

    # Make template
    template = RepairTemplate(sequence, 174, 'Q')
    template.make_repair_template()

    # Build Template
    output = Output('Home')
    output.set_sgRNA_sequences(*sgRNAs.get_sgRNA())
    output.set_template_sequences(*template.get_template())
    output.check_ouput_file()
