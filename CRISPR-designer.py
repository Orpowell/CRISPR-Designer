import re

from Bio import SeqIO

# Get nucleotide sequence of target protein
from Bio.Seq import Seq


def get_sequence(path):
    # Parse fasta file into Biopython Seq object
    for sequence in SeqIO.parse(path, "fasta"):
        print(f'Nucleotide Sequence:\n\n{sequence.seq}\n')  # Print Nucleotide sequence
        print(f'Protein Sequence:\n\n{sequence.seq.translate()}\n')  # Print protein sequence
        return sequence


# Search for all PAM sites in sequence, returns list of PAM sites
def find_PAM(sequence):
    i = 0
    positions = []  # List for pam sites positions

    # Search through sequence 2 base pairs at a time
    for _ in sequence:
        codon = sequence[i:i + 2]
        i += 1

        # If a PAM site is found store the position of the second base pair in PAM
        if codon == 'GG':
            pam_position = i + 1  # Store position of last base of PAM site
            positions.append(pam_position)  # Store value in a list

    return positions  # Return complete list of PAM sites in sequence


# Check for stretch of sequences containing 4 or more of the same nucleotide
def quadruple_nucleotide_check(x):
    p = re.compile("([T])\\1\\1\\1")  # Regular expression for TTTT
    a = re.search(p, x)  # Search given sequence for RE

    if a is not None:  # If RE is present returns false
        return False
    else:  # If RE not found return true
        return True


# Identify the 20 nucleotides upstream of a PAM site, check its contain less than 4 consecutive thymine bases
def find_20mers(site, sequence):
    identified_20mer = sequence[0:site][-23:]  # Store upstream sequence
    the_20mer = identified_20mer[:-3]  # Remove PAM site

    quadruple = quadruple_nucleotide_check(str(identified_20mer))  # Check for 4 consecutive thymine bases in 20mer

    # If no quadruple bases found return start and end position and sequence of 20mer
    if quadruple:
        return [site - 23, site - 3, the_20mer]


def find_60mers(site, sequence):
    identified_60mer = sequence[0:site][-63:]  # Store upstream sequence
    the_60mer = identified_60mer[:-3]  # Remove PAM site

    quadruple = quadruple_nucleotide_check(str(identified_60mer))  # Check for 4 consecutive thymine bases in 20mer

    # If no quadruple bases found return start and end position and sequence of 20mer
    if quadruple:
        return [site - 63, site - 3, the_60mer]


# Get position of target amino acid, identify corresponding base and desired mutation
def get_codon():
    residue_number = input('Input residue position: ')  # Collect amino acid position

    try:
        residue_number = int(residue_number)
        codon_start = (residue_number * 3) - 3  # Calculate start position of codon in nucleotide sequence
        codon = record.seq[codon_start:codon_start + 3]  # Identify full codon sequence (3 Base pairs)

        print(f'Target is: {codon.translate()} @ position {residue_number}')  # Print amino acid and location
        return codon, codon_start

    except ValueError:
        print('Error: position must be an integer')
        retry = get_codon()
        return retry


# Check the target codon is within 20 or 60 nucleotides of a PAM site
def check_in_oligomer(codon_start, i, mode=1):
    minimum_oligomer = result[i][0]
    maximum_oligomer = result[i][1]

    # if codon is within 20 nucleotides of a PAM site
    if mode == 1:
        if (codon_start < maximum_oligomer) & (minimum_oligomer < codon_start):
            return result[i][2]

    # if codon is within 60 nucleotides of a PAM site
    else:
        if (codon_start < maximum_oligomer) & (minimum_oligomer < codon_start):
            return [result[i][2], maximum_oligomer, minimum_oligomer]


# Find all 20mers that contain the target codon and determine best one (last one, so more in the middle), exit program if none are found.
def within_20mer(codon_position):
    # Generate list of all 20mers containing target codon
    possible = list(filter(lambda x: x is not None, [check_in_oligomer(codon_position, i) for i in range(len(result))]))

    # select and return best 20mer
    try:
        best_match = possible[-1]
        print('Codon is present within 20nt of a PAM site\n')
        return best_match, 0

    # If no 20mer is found, exit program
    except IndexError:
        print('gotta do the 60 bro')
        return None, 1


def within_60mer(codon_position):
    # Generate list of all 60mers containing target codon
    possible = list(
        filter(lambda x: x is not None, [check_in_oligomer(codon_position, i, mode=2) for i in range(len(result))]))

    # select and return best 60mer
    best_match = possible[-1][0]  # Return best sequence with codon within 60mer
    minim = possible[-1][2]  # return start position of 60mer
    maxim = possible[-1][1]  # return end position of 60mer

    print('Codon is present within 60nt of a PAM site\n')  # Confirm to user
    return best_match, maxim, minim


# Provide interface for mutating codon, and determines which base in the codon will mutate and to what base
def codon_mutator(codon, codon_position):
    # get position of base to mutate within codon (1,2,3)
    def get_position():
        base_position = input('Which nucleotide do you want to mutate? (1,2 or 3)')  # collect base position

        # Check input is an integer
        try:
            base_1 = int(base_position)

            # Check input is a number between 1 and 3 or ask for input again
            if base_1 not in [1, 2, 3]:
                print('Please select a number from 1-3')
                restart = get_position()
                return restart

            # Convert input from Base1 to Base0 and return value
            else:
                base_0 = base_1 - 1
                return base_0

        # If value is not an integer ask for input again
        except ValueError:
            print('Error: Please use integer')
            restart = get_position()
            return restart

    # Collect what mutation to make within codon
    def get_nucleotide():
        # Get the base for the desired mutation
        nucleotide_change = input('What nucleotide you do want to replace your target with? (A,T,G or C)')

        # Check input is A, T, G or C or ask for input again
        if nucleotide_change not in ['A', 'T', 'G', 'C']:
            print('Please select a base from A, T, G or C')
            restart = get_nucleotide()
            return restart

        # Return nucleotide if correct
        else:
            return nucleotide_change

    codon = str(codon)  # Convert codon to string from Seq object for manipulation
    print(codon, '\n123')  # Show user target codon174
    base = get_position()  # Collect position of mutation site within codon
    mutant_nucleotide = get_nucleotide()  # collect what mutation to make at site within codon

    codon_list = list(codon)  # Separate target codon into list of 3 characters
    codon_list[base] = mutant_nucleotide  # Mutate codon as requested
    mutant_codon = "".join(codon_list)  # Create string of mutant codon

    print(
        f"Mutating {Seq(codon).translate()} ---> {Seq(mutant_codon).translate()}\n")  # Show effect of mutation on amino acid output
    mutation_site = base + codon_position  # calculate site of mutation in global sequence

    return mutant_nucleotide, mutation_site  # return mutated nucleotide and position


# Generate sgRNA 60mer for insertion into pCAS
def make_sgRNA(oligo):
    front = 'CGGGTGGCGAATGGGACTTT'  # Front primer
    back = 'GTTTTAGAGCTAGAAATAGC'  # back primer

    sgRNA = front + oligo.lower() + back  # Generate complete 60mer
    print('sgRNA 60mer:', sgRNA)  # print complete 60mer


# Create 60mer repair template containing desired mutation
def make_repair_template(position, target, nucleotide, mode=1):
    # Used to make repair template if mutant codon is within 20 nt of a PAM site
    if mode == 1:
        mutant_seq = list(record.seq)  # convert entire protein sequence to list
        mutant_seq[target] = nucleotide.lower()  # implement mutation as lower case to highlight in final result
        repair_template = "".join(mutant_seq)[
                          position - 30:position + 30]  # Create 60mer with 30 nt either side of mutation
        print(f"Repair template: {repair_template}")  # Display repair template to user

    # Used to make repair template if mutant codon is within 60 nt of a PAM site
    else:
        mutant_seq = list(record.seq)  # convert entire protein sequence to list
        distance_from_20mer = (
                                          maximum - 20) - position  # calculate distance between mutation and position 20 nt from PAM site
        codon_start = distance_from_20mer + 6 + position  # identify codon 2 positions within 20 mer
        codon = record.seq[codon_start:codon_start + 3]  # Find codon

        print('Make synonymous mutation in this codon')  # Request user makes synonymous mutation
        mutant, site = codon_mutator(codon, codon_start)  # Make synonymous mutation
        mutant_seq[target] = nucleotide.lower()  # implement target mutation
        mutant_seq[site] = mutant.lower()  # implement synonymous mutation
        repair_template = "".join(mutant_seq)[minimum:maximum]  # Make repair template
        print(f"Repair template: {repair_template}")  # print repair template


if __name__ == '__main__':

    record = get_sequence("S288C_YMR202W_ERG2_coding.fsa")  # Convert Fasta file to Seq object

    pam_sites = find_PAM(record.seq)  # Identify PAM sites

    test = [find_20mers(x, record.seq) for x in pam_sites]  # Identify 20mers upstream of PAM sites

    result = list(filter(lambda x: x is not None, test))  # Remove None values

    target_codon, position = get_codon()  # Collect target codon information

    best_20mer, switch = within_20mer(position)  # Find best 20mer containing target codon

    if switch == 0:

        nucleotide, target = codon_mutator(target_codon, position)  # Allow user to mutate codon

        make_sgRNA(best_20mer)  # Make sgRNA 60mer

        make_repair_template(position, target, nucleotide)  # Make repair template 60mer

    else:

        test = [find_60mers(x, record.seq) for x in pam_sites]  # Identify 60mers upstream of PAM sites

        result = list(filter(lambda x: x is not None, test))  # Remove None values

        best_60mer, maximum, minimum = within_60mer(position)  # Need start/finish pos for synonymous mutation

        guide = str(best_60mer[-20:]).lower()  # Isolate sgRNA insert from 60mer

        nucleotide, target = codon_mutator(target_codon, position)  # allow user to mutate target codon

        make_repair_template(position, target, nucleotide, mode=0)  # maker repair template 60mer

        make_sgRNA(guide)  # Make sgRNA 60mer
