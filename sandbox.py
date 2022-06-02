
from Bio import SeqIO
from Bio.Seq import Seq

sequence = 'ATGAAGTTTTTCCCACTCCTTTTGTTGATTGGTGTTGTAGGCTACATTATGAACGTATTGTTCACTACCTGGTTGCCAACCAATTACATGTTCGATCCAAAAACTTTGAACGAAATATGTAACTCGGTGATTAGCAAACACAACGCAGCAGAAGGTTTATCCACTGAAGACCTGTTACAGGATGTCAGAGACGCACTTGCCTCTCATTACGGGGACGAATACATCAACAGGTACGTCAAAGAAGAATGGGTCTTCAACAATGCTGGTGGTGCGATGGGCCAAATGATCATCCTACACGCTTCCGTATCCGAGTACTTAATTCTATTCGGAACCGCTGTTGGTACTGAAGGGCACACAGGTGTTCACTTTGCTGACGACTATTTTACCATCTTACATGGTACGCAAATCGCAGCATTGCCATATGCCACTGAAGCCGAAGTTTACACTCCTGGTATGACTCATCACTTGAAGAAGGGATACGCCAAGCAATACAGCATGCCAGGTGGTTCCTTTGCCCTTGAATTGGCTCAAGGCTGGATTCCATGTATGTTGCCATTCGGGTTTTTGGACACTTTCTCCAGTACTCTTGATTTATACACTCTATATAGAACTGTCTACCTGACTGCCAGGGACATGGGTAAGAACTTGTTGCAAAACAAAAAGTTCTAA'

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

aa_target = 173

codon_list = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]

print(f'initial codon {codon_list[aa_target]}')

aa = Seq(codon_list[aa_target]).translate()

print(f'initial amino acid {aa}')

codon_list[aa_target] = (list(codon_table.keys())[list(codon_table.values()).index(aa)])

print(f'\nswitched codon {codon_list[aa_target]}')

new_aa = Seq(codon_list[aa_target]).translate()
print(f'switched amino acid {new_aa}')

print(174*3-3)
print(519//3)
print(520//3)
print(521//3)
print(522//3)

