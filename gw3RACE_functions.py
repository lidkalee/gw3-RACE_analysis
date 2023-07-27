# gw3RACE_functions.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re

## TEST RNA TYPE

pattern_mRNA = re.compile("spac|spap|spbp|spcc|spcp|spbc")
pattern_tRNA = re.compile("trna")
pattern_mito = re.compile("spmit")
pattern_rRNA = re.compile("rrna")
pattern_snRNA = re.compile("snrna")
pattern_snoRNA = re.compile("snorna")
pattern_ncRNA = re.compile("ncrna")


def test_RNA_type(feature_name):
    """tests for the type of tail """
    if pattern_mRNA.search(feature_name): # for example AATT
        return "mRNA"
    elif pattern_tRNA.search(feature_name): # for example AATT
        return "tRNA"
    elif pattern_mito.search(feature_name): # for example AATT
        return "mito_RNA"
    elif pattern_rRNA.search(feature_name): # for example AATT
        return "rRNA"
    elif pattern_snRNA.search(feature_name): # for example AATT
        return "snRNA"
    elif pattern_snoRNA.search(feature_name): # for example AATT
        return "snoRNA"
    elif pattern_ncRNA.search(feature_name): # for example AATT
        return "ncRNA"
    else: #only A's
        return "other_type"

assert test_RNA_type("spncrna") == "ncRNA"
assert test_RNA_type("sspcc") == "mRNA"
assert test_RNA_type("AAAAA") == "other_type"
assert test_RNA_type("TAATTT") == "other_type"

########################################### TAKE TAIL FROM CIGAR

def take_tail_fromcigar8(strand, cig, seq):
    '''
    Function returns tail seqence based on the CIGAR and sequence of read
    Input:
        cig: CIGAR from sam file after alignment qith soft-clipping  [str]
        seq: sequence of read started from 3' end, having posibble tail   [str]
    Output:
        sequence of the tail extracted from read sequence based on CIGAR [str]
    Assumptions:
        1. For tailed reads, CIGAR starts from \d+S (softcliping code). I's about the number
        of softclipped nucleotides. We assume that this is also the length of the polyA tail.
        2. Sequence seq starts from 3' tail of read.
    '''

    cigar2 = str(cig) # convert cigar to string
    if not "S" in cigar2: # return nothing if CIGAR does not contain S (soft-clipping)
        return ""
    else:  # analyze data if CIGAR starts from number and S (3' soft-clipping). If not - return nothing
        
        
        if strand =='+':      
            search_result_minus = re.findall("[0-9]+S$", cigar2)
            if search_result_minus: # find CIGARs starting from 3' soft-clipping
                number_S = int(search_result_minus[-1].split('S')[0])
                # extract number of soft-clipped nucleotides
                tail = seq[-number_S:] # extract tail (based of number of soft-clipped nucleotides)
                return tail
            else:
                return ""
        elif strand == '-':
            search_result_plus = re.findall("^[0-9]+S", cigar2)
            if search_result_plus: # find CIGARs starting from 3' soft-clipping
                number_S = int(cigar2.split('S')[0]) # extract number of soft-clipped nucleotides
                tail = seq[0:number_S] # extract tail (based of number of soft-clipped nucleotides)
                return tail
            else:
                return ""
        else:
             return ""
            
            
            
#########################

pattern_polyAU_forward3 = re.compile("^A{1,}T{2,}") # reverse transcribed
pattern_oligoU_forward3 = re.compile("^A{1,}$")
pattern_polyA_forward3 = re.compile("^T{2,}$")

pattern_polyAU_reverse3 = re.compile("A{2,}T{1,}$") # reverse transcribed
pattern_oligoU_reverse3 = re.compile("T{1,}$")
pattern_polyA_reverse3 = re.compile("A{2,}$") 

def test_tail_cigargrep8(cigar_grep, strand_R1, tail):
    """ Tests for the type of tail. Patterns are prepared for analysis of R2 read, reverse transcribed:
        base sequence start from the end of polyadenylation tail at 3' end of RNA,
        where T means adenine, and A meand uridine in the 3' tail"""
    if len(str(tail)) > 0:
        if cigar_grep=='cigar':
            if strand_R1 == '-':
                if "A" in tail:
                    if pattern_polyAU_forward3.fullmatch(tail): 
                        return "polyAU"
                    elif pattern_oligoU_forward3.fullmatch(tail): 
                        return "oligoU"
                    else:
                        return "mixed_tail"
                elif pattern_polyA_forward3.fullmatch(tail):
                        return "polyA"
                else: 
                    return "mixed_tail"
            elif strand_R1 == '+':
                if "T" in tail:
                    if pattern_polyAU_reverse3.fullmatch(tail): 
                        return "polyAU"
                    elif pattern_oligoU_reverse3.fullmatch(tail): 
                        return "oligoU"
                    else:
                        return "mixed_tail"
                elif pattern_polyA_reverse3.fullmatch(tail):
                        return "polyA"
                else: 
                    return "mixed_tail"
            else: 
                    return "strand_nn"
        elif cigar_grep=='grep':
            if strand_R1 == '+':
                if "A" in tail:
                    if pattern_polyAU_forward3.fullmatch(tail): 
                        return "polyAU"
                    elif pattern_oligoU_forward3.fullmatch(tail): 
                        return "oligoU"
                    else:
                        return "mixed_tail_grep"
                elif pattern_polyA_forward3.fullmatch(tail):
                        return "polyA"
                else: 
                    return "mixed_tail"
            elif strand_R1 == '-':
                if "A" in tail:
                    if pattern_polyAU_forward3.fullmatch(tail): 
                        return "polyAU"
                    elif pattern_oligoU_forward3.fullmatch(tail): 
                        return "oligoU"
                    else:
                        return "mixed_tail_grep"
                elif pattern_polyA_forward3.fullmatch(tail):
                        return "polyA"
                else: 
                    return "mixed_tail"
            else:
                return "nocigar_grep"
        else:
            return "without_tail"
    elif len(str(tail)) <1:
        if cigar_grep=='cigar':
            return "no_tail"
        elif cigar_grep =='grep':
            return "mixed_tail"
        else:
            return "jakis_bias"
            
    
    else: 
        return 'jakis_bias'
    
#######################

def grep_tail_edit_onlyfromSeq(seq):
    '''
    Function ....
    '''

    search_result_plus = re.findall("^[A]{0,}[T]{0,}", seq)
    if search_result_plus: # find CIGARs starting from 3' soft-clipping
        return search_result_plus[0]
    else:
        return ""


def tail_fromGREPorCIGAR_description(cigar, tailGrep, tailCigar):
    if cigar == '*':
        return 'grep'
    else:
        return 'cigar'
    

def tail_fromGREPorCIGAR(cigar, tailGrep, tailCigar):
    if cigar == '*':
        return tailGrep
    else:
        return tailCigar
    
##########################################333
import re

def flag_16_coordinate(cigar):
    cigar_sum_len = sum([int(s.strip('MND')) for s in re.findall(r'\d+[MND]', cigar)])
    return(cigar_sum_len)

def stop_based_on_cigar(cigar, strand_R1, coord_R2):
    if cigar == '*':
        return 0
    elif  cigar != '*':
        if strand_R1 == '-':
            stop_R2 = coord_R2 -1
            return stop_R2
        elif strand_R1 == '+' :
            distance = flag_16_coordinate(cigar)
            stop_R2 = coord_R2 -1 + distance
            return stop_R2
        else:
            return 'cigar_but_unknown_strand'
    else:
        return 0
        
def distance_to_TES(cigar, strand_R1, stop_R2, gene_start,gene_stop):
    if cigar == '*':
        return 0
    elif  cigar != '*':
        if strand_R1 == '-':
            distance_to_TES = gene_start - stop_R2
            return distance_to_TES
        elif strand_R1 == '+' :
            distance_to_TES = stop_R2 - gene_stop
            return distance_to_TES
        else:
            return 0
    else:
        return 0

