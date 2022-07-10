#!/usr/local/bin/python3

# Required Libraries for CRISPR designer
import re
import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq

# Import classes from CRISPR designer
from crisprdesigner import sgRNA
from crisprdesigner import RepairTemplate
from crisprdesigner import sgRNA
from crisprdesigner import Output
from crisprdesigner import SequencingPrimer


def test_ShortSequence():
    test_sequence = 'AAAAAGG'*20
    obj = sgRNA(test_sequence, 22)
    obj.make_sgRNAs()
    print(obj.get_sgRNA())

test_ShortSequence()
