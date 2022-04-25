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


class Designer:
    def __init__(self, sequence, amino_acid_position, amino_acid_mutation, output):
        # Initialise inputs
        self.sequence = sequence
        self.amino_acid_position = amino_acid_position
        self.codon_position = (self.amino_acid_position * 3) - 3
        self.amino_acid_mutation = amino_acid_mutation
        self.output = output

        # Data outputs
        self.sgRNA_forward = None
        self.sgRNA_reverse = None
        self.repair_template_sequence = None
        self.forward_primer_sequence = None
        self.reverse_primer_sequence = None

    # Make sgRNAs
    def make_sgRNAs(self):

        def quadruple_nucleotide_check(sequence):
            quadruplet = re.compile("([T])\\1\\1\\1")  # Regular expression for TTTT
            check = re.search(quadruplet, sequence)  # Search given sequence for RE

            if check is not None:  # If RE is present returns false
                return False
            else:  # If RE not found return true
                return True

        def find_20mers(site, sequence):
            identified_20mer = sequence[0:site][-23:]  # Store upstream sequence
            the_20mer = identified_20mer[:-3]  # Remove PAM site

            quadruple = quadruple_nucleotide_check(
                str(identified_20mer))  # Check for 4 consecutive thymine bases in 20mer

            # If no quadruple bases found return start and end position and sequence of 20mer
            if quadruple:
                return [the_20mer]

        # define search region for PAM sites within 22 nt upstream of codon
        search_region = self.sequence[self.codon_position: self.codon_position + 22]

        i = 0
        pam_sites = []
        for _ in search_region:
            codon = search_region[i:i + 2]
            i += 1

            # If a PAM site is found store the position of the second base pair in PAM
            if codon == 'GG':
                pam_position = i + 1  # Store position of last base of PAM site
                pam_sites.append(pam_position + self.codon_position)  # Store value in a list

        test = [find_20mers(x, self.sequence) for x in pam_sites]  # Identify 20mers upstream of PAM sites

        result = filter(lambda x: x is not None, test)  # Remove None values
        oligos = list(result)

        oligo = oligos[-1][0]

        front = 'CGGGTGGCGAATGGGACTTT'  # Front primer
        back = 'GTTTTAGAGCTAGAAATAGC'  # back primer

        sgRNA = front + oligo.lower() + back

        self.sgRNA_forward = sgRNA
        self.sgRNA_reverse = str(Seq(sgRNA).reverse_complement())

    # get sgRNAs
    def get_sgRNA(self):
        print(self.sgRNA_forward)
        print(self.sgRNA_reverse)

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

    def get_template(self):
        print(self.repair_template_sequence)
        print(self.forward_primer_sequence)
        print(self.reverse_primer_sequence)



def find_60mers(site, sequence):
    identified_60mer = sequence[0:site][-63:]  # Store upstream sequence
    the_60mer = identified_60mer[:-3]  # Remove PAM site

    quadruple = quadruple_nucleotide_check(str(identified_60mer))  # Check for 4 consecutive thymine bases in 20mer

    # If no quadruple bases found return start and end position and sequence of 20mer
    if quadruple:
        return [site - 63, site - 3, the_60mer]


def within_60mer(codon_position):
    # Generate list of all 60mers containing target codon
    possible = list(
        filter(lambda x: x is not None, [check_in_oligomer(codon_position, i, mode=2) for i in range(len(result))]))

    try:
        # select and return best 60mer
        best_match = possible[-1][0]  # Return best sequence with codon within 60mer
        minim = possible[-1][2]  # return start position of 60mer
        maxim = possible[-1][1]  # return end position of 60mer

        print('Codon is present within 60nt of a PAM site.\n')  # Confirm to user
        return best_match, maxim, minim

    except IndexError:
        print('Error: Target codon is not within 20 or 60 nucleotides of a PAM site...')


# Create 60mer repair template containing desired mutation
def make_repair_template(amino_acid_position, nucleotide_target, selected_nucleotide, mode=1):
    # Used to make repair template if mutant codon is within 20 nt of a PAM site
    if mode == 1:
        mutant_seq = list(record.seq)  # convert entire protein sequence to list
        mutant_seq[
            nucleotide_target] = selected_nucleotide.lower()  # implement mutation as lower case to highlight in final result
        mutant_gene = "".join(mutant_seq)
        repair_template_sequence = mutant_gene[
                                   amino_acid_position - 30:amino_acid_position + 30]  # Create 60mer with 30 nt either side of mutation

        forward_primer_sequence = mutant_gene[(amino_acid_position - 80):(amino_acid_position - 30)] + mutant_gene[
                                                                                                       (
                                                                                                               amino_acid_position - 30):(
                                                                                                               amino_acid_position - 10)]
        reverse_primer_sequence = Seq(mutant_gene[(amino_acid_position + 10):(amino_acid_position + 30)] + mutant_gene[(
                                                                                                                               amino_acid_position + 30):(
                                                                                                                               amino_acid_position + 80)]).reverse_complement()
        full_template = mutant_gene[(amino_acid_position - 80):(amino_acid_position + 80)]

    # Used to make repair template if mutant codon is within 60 nt of a PAM site
    else:
        mutant_seq = list(record.seq)  # convert entire protein sequence to list
        distance_from_20mer = (
                                      maximum - 20) - amino_acid_position  # calculate distance between mutation and position 20 nt from PAM site
        codon_start = distance_from_20mer + 6 + amino_acid_position  # identify codon 2 positions within 20 mer
        codon = record.seq[codon_start:codon_start + 3]  # Find codon

        print('Make synonymous mutation in this codon')  # Request user makes synonymous mutation
        mutant, site = codon_mutator(codon, codon_start)  # Make synonymous mutation
        mutant_seq[nucleotide_target] = selected_nucleotide.lower()  # implement target mutation
        mutant_seq[site] = mutant.lower()  # implement synonymous mutation
        mutant_gene = "".join(mutant_seq)
        repair_template_sequence = mutant_gene[minimum:maximum]  # Make repair template

        forward_primer_sequence = mutant_gene[(minimum - 50):minimum] + mutant_gene[minimum:(minimum + 20)]
        reverse_primer_sequence = Seq(
            mutant_gene[(maximum - 20):maximum] + mutant_gene[maximum:(maximum + 50)]).reverse_complement()
        full_template = mutant_gene[(minimum - 50):(maximum + 50)]

    return repair_template_sequence, forward_primer_sequence, reverse_primer_sequence, full_template


if __name__ == '__main__':
    '''
    record = get_sequence()  # Convert Fasta file to Seq object

    pam_sites = find_PAM(record.seq)  # Identify PAM sites

    test = [find_20mers(x, record.seq) for x in pam_sites]  # Identify 20mers upstream of PAM sites

    result = list(filter(lambda x: x is not None, test))  # Remove None values

    target_codon, position = get_codon()  # Collect target codon information

    best_20mer, switch = within_20mer(position)  # Find best 20mer containing target codon

    # Use 20mer method to generate sgRNAs and repair template
    if switch == 0:

        nucleotide, target = codon_mutator(target_codon, position)  # Allow user to mutate codon

        sgRNA_forward, sgRNA_reverse = make_sgRNA(best_20mer)  # Make sgRNA 60mers

        repair_template, forward_primer, reverse_primer, full_repair_template = make_repair_template(position, target,
                                                                                                     nucleotide)  # Make repair template 60mer and primers

    # Use 60mer method to generate sgRNAs and repair template
    else:

        test = [find_60mers(x, record.seq) for x in pam_sites]  # Identify 60mers upstream of PAM sites

        result = list(filter(lambda x: x is not None, test))  # Remove None values

        try:
            best_60mer, maximum, minimum = within_60mer(position)  # Need start/finish pos for synonymous mutation

        except TypeError:
            sys.exit()

        guide = str(best_60mer[-20:]).lower()  # Isolate sgRNA insert from 60mer

        nucleotide, target = codon_mutator(target_codon, position)  # allow user to mutate target codon

        repair_template, forward_primer, reverse_primer, full_repair_template = \
            make_repair_template(position, target,
                                 nucleotide,
                                 mode=0)  # maker repair template 60mer and primers

        sgRNA_forward, sgRNA_reverse = make_sgRNA(guide)  # Make sgRNA 60mers

    # Write all data to a .txt file in fasta format
    file_name = input('Please input file name: ') + '.txt'

    with open(file_name, 'w+') as file:
        file.write(f'{file_name}\n')
        file.write(f' \n> Protein sequence {len(record.seq.translate())} aa\n{record.seq.translate()}\n')
        file.write(f' \n> Nucleotide sequence {len(record.seq)} bp\n{record.seq}\n')
        file.write(f' \n> sgRNA forward primer {len(sgRNA_forward)} bp\n{sgRNA_forward} \n')
        file.write(f' \n> sgRNA reverse primer {len(sgRNA_reverse)} bp\n{sgRNA_reverse}\n')
        file.write(f' \n> repair template {len(repair_template)}\n{repair_template}\n')
        file.write(f" \n> forward repair template primer {len(forward_primer)} bp\n{forward_primer}\n")
        file.write(f' \n> reverse repair template primer {len(reverse_primer)} bp\n{reverse_primer}\n')
        file.write(f' \n> full repair template {len(full_repair_template)} bp\n{full_repair_template}\n')
    '''
    # S288C_YMR202W_ERG2_coding.fsa

    sequence = 'ATGAAGTTTTTCCCACTCCTTTTGTTGATTGGTGTTGTAGGCTACATTATGAACGTATTGTTCACTACCTGGTTGCCAACCAATTACATGTTCGATCCAAAAACTTTGAACGAAATATGTAACTCGGTGATTAGCAAACACAACGCAGCAGAAGGTTTATCCACTGAAGACCTGTTACAGGATGTCAGAGACGCACTTGCCTCTCATTACGGGGACGAATACATCAACAGGTACGTCAAAGAAGAATGGGTCTTCAACAATGCTGGTGGTGCGATGGGCCAAATGATCATCCTACACGCTTCCGTATCCGAGTACTTAATTCTATTCGGAACCGCTGTTGGTACTGAAGGGCACACAGGTGTTCACTTTGCTGACGACTATTTTACCATCTTACATGGTACGCAAATCGCAGCATTGCCATATGCCACTGAAGCCGAAGTTTACACTCCTGGTATGACTCATCACTTGAAGAAGGGATACGCCAAGCAATACAGCATGCCAGGTGGTTCCTTTGCCCTTGAATTGGCTCAAGGCTGGATTCCATGTATGTTGCCATTCGGGTTTTTGGACACTTTCTCCAGTACTCTTGATTTATACACTCTATATAGAACTGTCTACCTGACTGCCAGGGACATGGGTAAGAACTTGTTGCAAAACAAAAAGTTCTAA'

    obj1 = Designer(sequence, 174, 'Q', 'Home')
    obj1.make_sgRNAs()
    obj1.get_sgRNA()
    obj1.make_repair_template()
    obj1.get_template()
