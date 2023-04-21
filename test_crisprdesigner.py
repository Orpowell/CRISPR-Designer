import pytest

# Import classes and functions for testing
from crisprdesigner import open_fasta
from crisprdesigner import sgRNA
from crisprdesigner import RepairTemplate
from crisprdesigner import SequencingPrimer

'''
#######################
  Testing Information
#######################

sgRNA and RepairTemplate classes are tested using the sequence found in genomic_test.fsa. 
These sequences were created with v1.0.0 of crispr designer but have been shown to work experimentally 
and cover most use cases. Sequences were manually adpated to account for changes to sequence 
design in more recent versions. E.g codon mutation replacing single point mutations, and small
shifts in positions of sequences.

The RepairTemplate class is also tested using an altered sequence to check methioinine 
synonymous mutations are handled correctly by CRISPR Designer (found in genomic_test_meth.fsa)

SequencingPrimer is tested using synthetic data. This is because it was added much later in development and the 
sequences have yet to be experimentally validated.

'''

# Erg2 sequence for testing
erg2_sequence = open_fasta('genomic_test.fsa')

# Synthetic data for testing SequencingPrimer class
sequencing_test_sequence = 'A' * 100 + ('T' * 5 + 'C' * 5 + 'T' * 10) + (
        'G' * 60 + ('A' * 5 + 'C' * 5 + 'G' * 10) + 'C' * 280) + ('T' * 5 + 'C' * 5 + 'G' * 10) + 'A' * 100


# Test sequences are loaded correctly into crispr designer
@pytest.mark.open_fasta
def test_open_fasta():
    seq = open_fasta('genomic_test.fsa')
    assert seq == erg2_sequence


# Test sgRNA primers are designed correctly for 20mer implementation
@pytest.mark.sgRNA
def test_20mer_sgRNA_design():
    test = sgRNA(erg2_sequence, 819)
    test.make_sgRNAs()
    test_primers = ['CGGGTGGCGAATGGGACTTTcccttgaattggctcaaggcGTTTTAGAGCTAGAAATAGC',
                    'GCTATTTCTAGCTCTAAAACgccttgagccaattcaagggAAAGTCCCATTCGCCACCCG']

    assert ([*test.get_sgRNA()] == test_primers)


# Test sgRNA primers are designed correctly for 60mer implementation
@pytest.mark.sgRNA
def test_60mer_sgRNA_design():
    test = sgRNA(erg2_sequence, 702)
    test.make_sgRNAs()
    test_primers = ['CGGGTGGCGAATGGGACTTTgaagccgaagtttacactccGTTTTAGAGCTAGAAATAGC',
                    'GCTATTTCTAGCTCTAAAACggagtgtaaacttcggcttcAAAGTCCCATTCGCCACCCG']
    assert ([*test.get_sgRNA()] == test_primers)


# Test 60mer is activated correctly in sgRNA class
@pytest.mark.sgRNA
def test_60mer_sgRNA_switch():
    test = sgRNA(erg2_sequence, 702)
    test.make_sgRNAs()
    assert (test.get_switch() == 250)


# Test core and full templates are correctly designed for 20mer implementation
@pytest.mark.RepairTemplate
def test_20mer_core_and_full_template():
    test = RepairTemplate(erg2_sequence, 174, 'M', 819)
    test.design_template()
    test_repair_core = 'TACAGCATGCCAGGTGGTTCCTTTGCCCTTatgTTGGCTCAAGGCTGGATTCCATGTATG'
    test_full_teplate = 'TTTACACTCCTGGTATGACTCATCACTTGAAGAAGGGATACGCCAAGCAATACAGCATGCCAGGTGGTTCCTTTGCCCTTatgTTGGCTCAAGGCTGGATTCCATGTATGTTGCCATTCGGGTTTTTGGACACTTTCTCCAGTACTCTTGATTTATACAC'

    assert ([*test.get_template()][-2:] == [test_repair_core, test_full_teplate])


# Test repair template primers are correctly designed for 20mer implementation
@pytest.mark.RepairTemplate
def test_20mer_template_primers():
    test = RepairTemplate(erg2_sequence, 174, 'M', 819)
    test.design_template()
    test_forward_primer = 'TTTACACTCCTGGTATGACTCATCACTTGAAGAAGGGATACGCCAAGCAATACAGCATGCCAGGTGGTTC'
    test_reverse_primer = 'GTGTATAAATCAAGAGTACTGGAGAAAGTGTCCAAAAACCCGAATGGCAACATACATGGAATCCAGCCTT'

    assert ([*test.get_template()][:2] == [test_forward_primer, test_reverse_primer])


# Test core templates is correctly designed for 60mer implementation
@pytest.mark.RepairTemplate
def test_60mer_core_template():
    test = RepairTemplate(erg2_sequence, 135, 'M', 702)
    test.set_switch(250)
    test.design_template()

    test_repair_core = 'TTACATGGTACGatgATCGCAGCATTGCCATATGCCACTGAAGCCgagGTTTACACTCCT'

    assert ([*test.get_template()][2] == test_repair_core)


# Test full templates is correctly designed for 60mer implementation
@pytest.mark.RepairTemplate
def test_60mer_full_template():
    test = RepairTemplate(erg2_sequence, 135, 'M', 702)
    test.set_switch(250)
    test.design_template()
    test_full_teplate = 'GTACTGAAGGGCACACAGGTGTTCACTTTGCTGACGACTATTTTACCATCTTACATGGTACGatgATCGCAGCATTGCCATATGCCACTGAAGCCgagGTTTACACTCCTGGTATGACTCATCACTTGAAGAAGGGATACGCCAAGCAATACAGCATGCC'

    assert ([*test.get_template()][3] == test_full_teplate)


# Test forward template primer is correctly designed for 60mer implementation
@pytest.mark.RepairTemplate
def test_60mer_forward_template_primers():
    test = RepairTemplate(erg2_sequence, 135, 'M', 702)
    test.set_switch(250)
    test.design_template()
    test_forward_primer = 'GTACTGAAGGGCACACAGGTGTTCACTTTGCTGACGACTATTTTACCATCTTACATGGTACGatgATCGC'
    print(test.get_template())
    assert ([*test.get_template()][0] == test_forward_primer)


# Test reverse template primer is correctly designed for 60mer implementation
@pytest.mark.RepairTemplate
def test_60mer_reverse_template_primers():
    test = RepairTemplate(erg2_sequence, 135, 'M', 702)
    test.set_switch(250)
    test.design_template()

    test_reverse_primer = 'GGCATGCTGTATTGCTTGGCGTATCCCTTCTTCAAGTGATGAGTCATACCAGGAGTGTAAACctcGGCTT'
    print(test.get_template())
    assert ([*test.get_template()][1] == test_reverse_primer)


# Methionine synonymous mutation test case: full template
@pytest.mark.RepairTemplateMethionine
@pytest.mark.RepairTemplate
def test_methioine_full_template():
    test_full_teplate = 'GTACTGAAGGGCACACAGGTGTTCACTTTGCTGACGACTATTTTACCATCTTACATGGTACGatgATCGCAGCATTGCCATATGCCACTGAAGCCATGgtaTACACTCCTGGTATGACTCATCACTTGAAGAAGGGATACGCCAAGCAATACAGCATGCC'
    test = RepairTemplate(open_fasta('genomic_test_meth.fsa'), 135, 'M', 702)
    test.set_switch(250)
    test.design_template()
    assert ([*test.get_template()][3] == test_full_teplate)


# Methionine synonymous mutation test case: core template
@pytest.mark.RepairTemplateMethionine
@pytest.mark.RepairTemplate
def test_methioine_core_template():
    test_core_teplate = 'GTACTGAAGGGCACACAGGTGTTCACTTTGCTGACGACTATTTTACCATCTTACATGGTACGatgATCGC'
    test = RepairTemplate(open_fasta('genomic_test_meth.fsa'), 135, 'M', 702)
    test.set_switch(250)
    test.design_template()
    assert ([*test.get_template()][0] == test_core_teplate)


# Methionine synonymous mutation test case: forward primer
@pytest.mark.RepairTemplateMethionine
@pytest.mark.RepairTemplate
def test_methioine_forward_primer_template():
    test_forward_primer = 'GGCATGCTGTATTGCTTGGCGTATCCCTTCTTCAAGTGATGAGTCATACCAGGAGTGTAtacCATGGCTT'
    test = RepairTemplate(open_fasta('genomic_test_meth.fsa'), 135, 'M', 702)
    test.set_switch(250)
    test.design_template()
    assert ([*test.get_template()][1] == test_forward_primer)


# Methionine synonymous mutation test case: reverse primer
@pytest.mark.RepairTemplateMethionine
@pytest.mark.RepairTemplate
def test_methioine_reverse_primer_template():
    test_reverse_primer = 'TTACATGGTACGatgATCGCAGCATTGCCATATGCCACTGAAGCCATGgtaTACACTCCT'
    test = RepairTemplate(open_fasta('genomic_test_meth.fsa'), 135, 'M', 702)
    test.set_switch(250)
    test.design_template()
    assert ([*test.get_template()][2] == test_reverse_primer)


# Test amplified region is correctly defined
@pytest.mark.SequencingPrimer
def test_amplified_region():
    test = SequencingPrimer(sequencing_test_sequence, 300)
    test.make_region_primers()
    assert ([*test.get_primers()][0] == sequencing_test_sequence[100:500])


# Test amplified region primers are correctly designed
@pytest.mark.SequencingPrimer
def test_amplified_region_primers():
    test_amplified_forward_primer = ('T' * 5 + 'C' * 5 + 'T' * 10)
    test_amplified_reverse_primer = ('C' * 10 + 'G' * 5 + 'A' * 5)

    test = SequencingPrimer(sequencing_test_sequence, 300)
    test.make_region_primers()
    assert ([*test.get_primers()][1:3] == [test_amplified_forward_primer, test_amplified_reverse_primer])


# Test sequencing primer for the amplified region is correctly designed
@pytest.mark.SequencingPrimer
def test_sequencing_primer():
    seq_primer = ('A' * 5 + 'C' * 5 + 'G' * 10)
    test_seq_primer = SequencingPrimer(sequencing_test_sequence, 300)
    test_seq_primer.make_seq_primer()
    assert ([*test_seq_primer.get_primers()][3] == seq_primer)
