file_name = input('Please input file name: ') + '.txt'

protein_sequence = "ADSGHKFKLRHSKDIGDKDPGSBDKFJADYBFDSKHGDFIGUOFADIHAFSUVFIHUFDHJIHAS"
nucleotide_sequence = protein_sequence*3
sgRNA_forward = 'AbC'
sgRNA_reverse = 'CbA'

repair_template = '123ABC456'
forward_primer = '123123'
reverse_primer = '456456'
full_repair_template = '123123ABC456456'

with open(file_name, 'w+') as file:
    file.write(f'> Protein sequence\n{protein_sequence}\n')
    file.write(f' \n> Nucleotide sequence\n{nucleotide_sequence}\n')
    file.write(f' \n> sgRNA forward\n{sgRNA_forward} \n')
    file.write(f' \n> sgRNA reverse\n{sgRNA_reverse}\n')
    file.write(f' \n> repair template\n{repair_template}\n')
    file.write(f' \n> forward repair template primer\n{forward_primer}\n')
    file.write(f' \n> reverse repair template primer\n{reverse_primer}\n')
    file.write(f' \n> full repair template\n{full_repair_template}\n')
