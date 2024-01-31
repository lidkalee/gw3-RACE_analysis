# gw3RACE_functions.py

import pandas as pd
import matplotlib.pyplot as plt
import re

##################### TEST RNA TYPE  #####################

pattern_mRNA = re.compile("spac|spap|spbp|spcc|spcp|spbc")
pattern_tRNA = re.compile("trna")
pattern_mito = re.compile("spmit")
pattern_rRNA = re.compile("rrna")
pattern_snRNA = re.compile("snrna")
pattern_snoRNA = re.compile("snorna")
pattern_ncRNA = re.compile("ncrna")


def test_RNA_type(feature_name):
    '''

    The function determines the RNA type based on the given feature name using predefined regular expression patterns.
    
    The function searches the feature name for specific patterns that match various RNA types, including mRNA, tRNA, mitochondrial RNA (mito_RNA), rRNA, snRNA, snoRNA, and ncRNA. If the feature name does not match any of the predefined patterns, it defaults to "other_type".
    
    Parameters:
    - feature_name (str): The name or identifier of the feature to be classified. This is typically a string that may include gene identifiers or RNA sequence names.
    
    Returns:
    - str: The determined type of RNA based on the feature name. Possible return values are "mRNA", "tRNA", "mito_RNA", "rRNA", "snRNA", "snoRNA", "ncRNA", or "other_type".
    
    Examples:
    - test_RNA_type("spac1234") returns "mRNA"
    - test_RNA_type("trnaMet") returns "tRNA"
    - test_RNA_type("spmit1234") returns "mito_RNA"
    - test_RNA_type("rrna12") returns "rRNA"
    - test_RNA_type("snrnaU1") returns "snRNA"
    - test_RNA_type("snorna123") returns "snoRNA"
    - test_RNA_type("ncrnaXist") returns "ncRNA"
    - test_RNA_type("unknown123") returns "other_type"
    
    Note:
    The function uses predefined regular expressions to match the feature name against known RNA type patterns. These patterns are designed based on common naming conventions and may need to be updated or expanded based on the specific dataset or nomenclature used.
    '''
    if pattern_mRNA.search(feature_name):
        return "mRNA"
    elif pattern_tRNA.search(feature_name):
        return "tRNA"
    elif pattern_mito.search(feature_name): 
        return "mito_RNA"
    elif pattern_rRNA.search(feature_name): 
        return "rRNA"
    elif pattern_snRNA.search(feature_name):
        return "snRNA"
    elif pattern_snoRNA.search(feature_name): 
        return "snoRNA"
    elif pattern_ncRNA.search(feature_name): 
        return "ncRNA"
    else: 
        return "other_type"

assert test_RNA_type("spncrna") == "ncRNA"
assert test_RNA_type("sspcc") == "mRNA"
assert test_RNA_type("AAAAA") == "other_type"
assert test_RNA_type("TAATTT") == "other_type"

##################### TAKE TAIL FROM CIGAR  #####################

def take_tail_fromcigar8(strand, cig, seq):
    """
    Extracts the tail sequence from a read based on the CIGAR string and the actual read sequence, taking into account the strand orientation.

    The function identifies soft-clipped regions in the CIGAR string, which indicate the presence of a tail in the read. The tail sequence is then extracted from the corresponding end of the read sequence.

    Parameters:
    - strand (str): The strand orientation of the read, indicated by '+' for forward and '-' for reverse.
    - cig (str): The CIGAR string obtained from the SAM file after alignment, which includes soft-clipping information.
    - seq (str): The nucleotide sequence of the read, starting from the 3' end, which may include a tail sequence.

    Returns:
    - str: The extracted tail sequence based on the soft-clipping information from the CIGAR string. Returns an empty string if no tail is identified or if the CIGAR string does not contain soft-clipping information.

    Assumptions:
    - Soft-clipped nucleotides (indicated by 'S' in the CIGAR string) at the beginning or end of the read are considered part of the tail sequence.
    - The function assumes the read sequence (`seq`) is provided in 5' to 3' orientation relative to the strand.

    Notes:
    - The function does not account for complex scenarios where soft-clipping occurs in the middle of the CIGAR string or for reads with multiple soft-clipped regions.
    - It is essential that the input CIGAR string and read sequence are correctly aligned and correspond to the same read.
    """
    cigar2 = str(cig)  # Convert CIGAR to string for processing
    if "S" not in cigar2:  # Check if CIGAR string contains soft-clipping
        return ""  # Return empty string if no soft-clipping is found

    # Process the CIGAR string based on strand orientation to extract the tail sequence
    if strand == '+':
        search_result_minus = re.findall("[0-9]+S$", cigar2)
        if search_result_minus:
            number_S = int(search_result_minus[-1].split('S')[0])  # Number of soft-clipped nucleotides at the end
            tail = seq[-number_S:]  # Extract tail from the end of the sequence
            return tail
    elif strand == '-':
        search_result_plus = re.findall("^[0-9]+S", cigar2)
        if search_result_plus:
            number_S = int(search_result_plus[0].split('S')[0])  # Number of soft-clipped nucleotides at the start
            tail = seq[:number_S]  # Extract tail from the start of the sequence
            return tail

    return ""  # Return empty string if no tail is identified or for unsupported strand orientations
            
            

##################### TEST TAIL TYPE  #####################

pattern_polyAU_forward3 = re.compile("^A{1,}T{2,}") # reverse transcribed
pattern_oligoU_forward3 = re.compile("^A{1,}$")
pattern_polyA_forward3 = re.compile("^T{2,}$")

pattern_polyAU_reverse3 = re.compile("A{2,}T{1,}$") # reverse transcribed
pattern_oligoU_reverse3 = re.compile("T{1,}$")
pattern_polyA_reverse3 = re.compile("A{2,}$") 

def test_tail_cigargrep8(cigar_grep, strand_R1, tail):
    """
    Determines the type of tail in a sequence based on regex patterns, considering the sequence's strand orientation and whether the analysis is based on CIGAR or grep.

    Parameters:
    - cigar_grep (str): Indicates the method used for tail identification ('cigar' or 'grep').
    - strand_R1 (str): The strand orientation ('+' for forward, '-' for reverse).
    - tail (str): The sequence being analyzed for tail content.

    Returns:
    - str: The determined tail type. Possible values include "polyAU", "oligoU", "polyA", "mixed_tail", "strand_nn", "mixed_tail_grep", "without_tail", "no_tail", or "bias".

    The function uses several regex patterns to identify different tail types:
    - "polyAU": A sequence starting with one or more 'A's followed by two or more 'T's (or vice versa for the reverse strand).
    - "oligoU": A sequence composed entirely of 'A's (interpreted as uridines in the context of RNA).
    - "polyA": A sequence of two or more 'T's (adenines in the context of RNA).
    - "mixed_tail": A sequence that does not match any of the above patterns.
    - "strand_nn": Returned if the strand orientation is not recognized.
    - "mixed_tail_grep": A special case of "mixed_tail" when the analysis is based on grep.
    - "without_tail": Returned if the sequence is empty or does not contain a tail.
    - "no_tail": Returned if the analysis is based on CIGAR and no tail is detected.
    - "bias": A fallback value indicating an unhandled case or error.

    Notes:
    - The function assumes the tail sequence is provided in the 5' to 3' orientation relative to the strand.
    - The interpretation of 'A' and 'T' in the tail is context-dependent, considering the reverse transcription process where 'T' represents adenine and 'A' represents uridine in the RNA sequence.
    """
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
            return "bias"
            
    
    else: 
        return 'bias'
    
##################### TAKE TAIL USING GREP  #####################

def grep_tail_edit_onlyfromSeq(seq):
    """
    Extracts a leading sequence of 'A's followed by 'T's from the given sequence. This function is useful for identifying simple poly-A tails or mixed A/T tails at the beginning of a sequence, which are common in RNA sequencing data.

    Parameters:
    - seq (str): The nucleotide sequence to be analyzed, expected to start from the 3' end.

    Returns:
    - str: The extracted sequence of 'A's followed by 'T's from the start of the given sequence. Returns an empty string if no such pattern is found.

    Note:
    - This function specifically looks for sequences starting with 'A's followed by 'T's, which could represent poly-A tails or regions with mixed A and T nucleotides. It does not account for 'T's followed by 'A's or other variations.
    - The function uses regular expressions to identify the pattern and is designed to work with sequences in the 5' to 3' orientation.
    - The function returns only the first match found, which is the entire sequence of 'A's followed by 'T's at the start of the given sequence if present.
    """
    # Regular expression to find a sequence of 'A's optionally followed by 'T's at the beginning of the sequence
    search_result_plus = re.findall("^[A]{0,}[T]{0,}", seq)

    # Return the found sequence, or an empty string if no match is found
    if search_result_plus:
        return search_result_plus[0]
    else:
        return ""



def tail_fromGREPorCIGAR_description(cigar, tailGrep, tailCigar):
    """
    Determines the method used to identify the tail sequence based on the CIGAR string.

    Parameters:
    - cigar (str): The CIGAR string associated with the sequence. A '*' in the CIGAR string indicates that the tail was identified using the grep method.
    - tailGrep (str): The tail sequence identified using grep.
    - tailCigar (str): The tail sequence identified using CIGAR string analysis.

    Returns:
    - str: 'grep' if the tail was identified using grep (indicated by '*' in the CIGAR string); otherwise, 'cigar'.

    Note:
    - This function does not use `tailGrep` and `tailCigar` parameters in its current implementation but requires them for consistency with `tail_fromGREPorCIGAR` function or potential future modifications.
    """
    return 'grep' if cigar == '*' else 'cigar'
    

def tail_fromGREPorCIGAR(cigar, tailGrep, tailCigar):
    """
    Selects the appropriate tail sequence based on the method used for its identification, indicated by the CIGAR string.

    Parameters:
    - cigar (str): The CIGAR string associated with the sequence. A '*' in the CIGAR string indicates that the tail was identified using the grep method.
    - tailGrep (str): The tail sequence identified using grep.
    - tailCigar (str): The tail sequence identified using CIGAR string analysis.

    Returns:
    - str: The tail sequence identified by grep if the CIGAR string is '*', otherwise the tail sequence identified by CIGAR string analysis.

    Note:
    - This function uses the presence of '*' in the CIGAR string as an indicator that the grep method was used to identify the tail sequence. Otherwise, it assumes the tail was identified through CIGAR string analysis.
    """
    return tailGrep if cigar == '*' else tailCigar
    
##################### OTHER FUNCTIONS  #####################


def flag_16_coordinate(cigar):
    """
    Calculates the sum of lengths for specific operations in a CIGAR string.

    This function processes a CIGAR string and sums up the lengths of match (M), deletion (D), and skipped region (N) operations to calculate a total length. This length can be used to determine genomic coordinates from sequencing data.

    Parameters:
    - cigar (str): The CIGAR string from the alignment of a read to a reference genome.

    Returns:
    - int: The sum of lengths for match, deletion, and skipped region operations in the CIGAR string.

    Note:
    - The function assumes that the input CIGAR string is well-formed and contains only the operations it is designed to process (M, D, N).
    - This function does not account for insertions (I) and other operations that do not consume reference sequence.
    """
    cigar_sum_len = sum([int(s.strip('MND')) for s in re.findall(r'\d+[MND]', cigar)])
    return cigar_sum_len

def stop_based_on_cigar(cigar, strand_R1, coord_R2):
    """
    Calculates the stop position of a read based on its CIGAR string, taking into account the strand orientation.

    Parameters:
    - cigar (str): The CIGAR string from the alignment of the read to a reference genome. A '*' indicates an undefined CIGAR string.
    - strand_R1 (str): The strand orientation of the read ('+' for forward, '-' for reverse).
    - coord_R2 (int): The starting coordinate of the read on the reference genome.

    Returns:
    - int or str: The calculated stop position of the read on the reference genome. Returns 'cigar_but_unknown_strand' if the strand is not recognized ('+' or '-'). Returns 0 if the CIGAR string is undefined ('*') or in case of other errors.

    Note:
    - The function uses `flag_16_coordinate` to calculate the distance covered by the read on the reference genome based on the CIGAR string.
    - For reads on the '-' strand, the function assumes the stop position is one less than the starting coordinate. For '+' strand reads, it adds the calculated distance to the starting coordinate (minus one) to find the stop position.
    - This function may not handle complex cases such as reads with multiple segments (chimeric alignments) or other unusual mapping scenarios.
    """
    if cigar == '*':
        return 0
    else:
        if strand_R1 == '-':
            return coord_R2 - 1
        elif strand_R1 == '+':
            distance = flag_16_coordinate(cigar)
            return coord_R2 - 1 + distance
        else:
            return 'cigar_but_unknown_strand'

##################### DISTANCE OF 3'END TO TES  #####################
        
def distance_to_TES(cigar, strand_R1, stop_R2, gene_start, gene_stop):
    """
    Calculates the distance from the 3' end of a read to the transcription end site (TES) of a gene.

    The function takes into account the strand orientation of the read and the gene's start and stop coordinates to calculate the distance. A '*' in the CIGAR string indicates that the read is not properly aligned, and thus the distance cannot be calculated accurately.

    Parameters:
    - cigar (str): The CIGAR string from the alignment of the read. A '*' indicates an unaligned read or undefined CIGAR string.
    - strand_R1 (str): The strand orientation of the read ('+' for forward, '-' for reverse).
    - stop_R2 (int): The stop position of the read on the reference genome.
    - gene_start (int): The start position of the gene on the reference genome.
    - gene_stop (int): The stop position of the gene on the reference genome.

    Returns:
    - int: The calculated distance from the 3' end of the read to the TES of the gene. Returns 0 for unaligned reads or in cases where the strand orientation is not specified.

    Note:
    - For reads on the '-' strand, the distance is calculated as the gene's start position minus the read's stop position.
    - For reads on the '+' strand, the distance is calculated as the read's stop position minus the gene's stop position.
    - This function assumes that the input positions and orientations are correctly provided and that the read and gene are on the same chromosome or contig.
    """
    if cigar == '*':  # Check for unaligned reads or undefined CIGAR strings
        return 0
    else:
        if strand_R1 == '-':
            return gene_start - stop_R2  # Calculate distance for '-' strand
        elif strand_R1 == '+':
            return stop_R2 - gene_stop  # Calculate distance for '+' strand
        else:
            return 0  # Return 0 for undefined strand orientations or other errors


