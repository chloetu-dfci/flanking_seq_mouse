#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 15:21:00 2022

@author: chloetu
"""


import sys
import pandas as pd
import re
import requests
from math import trunc

columns = ['Transcript','HGVSc','HGVSp']

# Current script: 
#   Step 1: Mutated genes
# Obtains the cDNA sequence of a gene from Ensembl (the difference between 
# cDNA and genomic sequence is that the cDNA sequences do not contain introns) using transcript ID.
# Creates a cDNA sequence to reflect the mutation in the gene.
# Finds the exonic sequences from the cDNA sequence in both reference and mutated 
# and turns the protein coding sequence into codons. 
# Tells you if you have gone into the UTRs (lower cased nucleotides, when mask_feature = 1) 
# Using the codons, finds the amino acid created by the reference and mutant sequence.
# Cross-checks the reference and mutant amino acids generated from Ensembl with the 
# reference and mutant amino acids from Neil's mouse neoantigen pipeline (pvactools).
# Extracts around 100~ nucleotides upstream and downstream of the mutation, while honoring the
# correct open reading frame.

    

#######

# Find amino acid corresponding to given codon
def map_codon_to_aa(codon, aa_codon_table):
    matching_codon = aa_codon_table[aa_codon_table['codons'] == codon]
    matching_aa = matching_codon.iloc[0,0]
    return matching_aa

# Import amino acid/codon translation file
def import_aa_codon_info():
    aa_codon_table = pd.read_csv("/Users/chloetu/Desktop/mouse_epitope_extraction/amino_acid_codons.tsv", sep="\t")
    return aa_codon_table


# Return cdna gene sequence string using gene ID
def get_cdna_sequence_from_transcript_id(ensembl_transcript_id):
    # Get gene sequence from ensembl
    server = "https://nov2020.rest.ensembl.org/"
    ext = "/sequence/id/" + ensembl_transcript_id + "?content-type=text/plain;mask_feature=1;type=cdna"
    r = requests.get(server+ext)
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    return r.text


# Not used
# Return genomic gene sequence string using gene ID
# def get_genomic_sequence_from_transcript_id(ensembl_transcript_id):
#     # Get gene sequence from ensembl
#     server = "https://nov2020.rest.ensembl.org/"
#     ext = "/sequence/id/" + ensembl_transcript_id + "?content-type=text/plain;mask_feature=1;type=genomic"
#     r = requests.get(server+ext)
    
#     if not r.ok:
#         r.raise_for_status()
#         sys.exit()

#     return r.text

# Updates gene sequence from ensembl to reflect the variant we see in MC38 sample
def get_mutant_sequence(exon_sequence, cdna_variant, mutant_nucleotide_position):
    cdna_variant = str(cdna_variant)
    if ">" in cdna_variant:
        mutant_nucleotide = cdna_variant.split(">")[1]
        exon_sequence = exon_sequence[:mutant_nucleotide_position] + mutant_nucleotide + exon_sequence[mutant_nucleotide_position+1:]
    elif "del" in cdna_variant:
        exon_sequence = exon_sequence[:mutant_nucleotide_position] + exon_sequence[mutant_nucleotide_position+1:]
    elif "dup" in cdna_variant:
        exon_sequence = exon_sequence[:mutant_nucleotide_position] + exon_sequence[mutant_nucleotide_position:mutant_nucleotide_position + 1] + exon_sequence[mutant_nucleotide_position:]
    elif "ins" in cdna_variant:
        mutant_nucleotide = cdna_variant.split("ins")[1]
        exon_sequence = exon_sequence[:mutant_nucleotide_position] + mutant_nucleotide + exon_sequence[mutant_nucleotide_position:]
    else:
        print("Variant was not found")
    return exon_sequence

# Remove version (.X) from transcript ID
def remove_version_from_transcript_id(ensembl_gene_id):
    ensembl_gene_id = ensembl_gene_id.split(".")[0]
    return ensembl_gene_id

# Get coding region (Capital letters)
def get_coding_sequence(exon_sequence):
    coding_sequence = re.search("[A-Z]+", exon_sequence).group()
    return coding_sequence

# Convert exon sequence into a list of triplets (codons). 
def convert_exon_sequence_into_codons(exon_sequence):
    CDS_begins_with_ATG_or_CTG_flag = True
    
    codon_list = []
    match_obj = re.search("[A-Z]+", exon_sequence)

    coding_sequence = match_obj.group()
    # Check if the ATG is at the start of the CDS
    atg_index = re.search("^ATG", exon_sequence)
    # Pass a flag if the coding sequence does not begin with ATG
    if atg_index is None:
        ctg_index = re.search("^CTG", exon_sequence)
        if ctg_index is None:
            CDS_begins_with_ATG_or_CTG_flag = False
    
    for first_nt in range(0, len(coding_sequence), 3):
        codon = coding_sequence[first_nt: first_nt+3]
        if len(codon) == 3:
            codon_list.append(codon)
    return codon_list, CDS_begins_with_ATG_or_CTG_flag

# Find index referencing the codon for a specific nucleotide 
def extract_codon_index(mutant_nucleotide_position):
    # Find which codon is the mutant nucleotide is in
    codon_index = trunc(mutant_nucleotide_position / 3)
    return codon_index

# Extract nucleotide sequence from codon
def extract_codon(mutant_nucleotide_position, codon_list):
    codon_index = extract_codon_index(mutant_nucleotide_position)
    codon_of_interest = codon_list[codon_index]
    return codon_of_interest, codon_index


# # Actually extract the exons, account for being too close to either side. UNUSED
# def extract_flanking_exon_list():
#     codon_diff = -1
#     # Extract X codons from above and below the mutation
#     # If the mutation is too close to the beginning of the exon sequence
#     if flanking_codon_lengths >= codon_index:
#         codon_diff = flanking_codon_lengths-codon_index+1
#         location_flag = str("Sequence was shifted down by " + str(codon_diff) + " codons")
#         flanking_exon_list = codon_list[0: codon_index + flanking_codon_lengths + codon_diff]
#     # If the mutation is too close to the end of the exon sequence
#     elif (flanking_codon_lengths + codon_index) >= len(codon_list):
#         codon_diff = len(codon_list) - (flanking_codon_lengths + codon_index) - 1
#         location_flag = str("Sequence was shifted up by " + str(abs(codon_diff)) + " codons")
#         flanking_exon_list = codon_list[codon_index - flanking_codon_lengths + codon_diff: len(codon_list)]
#     else:
#         flanking_exon_list = codon_list[codon_index - flanking_codon_lengths: codon_index + flanking_codon_lengths + 1]
#     return flanking_exon_list, location_flag, codon_diff


# Actually extract the exons
def extract_flanking_exon_list(codon_list, upstream_codon_index, downstream_codon_index):
    flanking_exon_list = codon_list[upstream_codon_index: downstream_codon_index]
    return flanking_exon_list

def update_flanking_region_indexes(flanking_codon_lengths, codon_index, codon_list_length):
    location_flag = "OK"
    codon_diff = 0
    
    # if too close to the beginning, the codon_diff is positive.
    if flanking_codon_lengths >= codon_index:
        codon_diff = flanking_codon_lengths-codon_index+1
        upstream_codon_index = 0
        downstream_codon_index = codon_index + flanking_codon_lengths + codon_diff
        location_flag = "Mutation is too close to sequence beginning"
    # if too close to the end, the codon_diff is negative.
    elif (flanking_codon_lengths + codon_index) >= codon_list_length:
        codon_diff = codon_list_length - (flanking_codon_lengths + codon_index) - 1
        upstream_codon_index = codon_index - flanking_codon_lengths + codon_diff
        downstream_codon_index = codon_list_length
        location_flag = "Mutation is too close to sequence end"
    else:
        upstream_codon_index = codon_index - flanking_codon_lengths
        downstream_codon_index = codon_index + flanking_codon_lengths + 1
        codon_diff = 0
    return upstream_codon_index, downstream_codon_index, codon_diff, location_flag


def find_bsai_codon_indexes(exon_sequence):
    bsai_nt_list = find_bsai_nt_indexes(exon_sequence)
    bsai_start_codon_list = []
    bsai_end_codon_list = []
    for i in range(len(bsai_nt_list)):
        bsai_start_codon_list.append(extract_codon_index(bsai_nt_list[i]))
        bsai_end_codon_list.append(extract_codon_index(bsai_nt_list[i]+5))
    return bsai_start_codon_list, bsai_end_codon_list


# Check where the exon has a BsaI site: The recognition sequence is GGTCTC (a non-palindromic cutter)
def find_bsai_nt_indexes(exon_string):
    bsai_nt_list = [index.start() for index in re.finditer("GGTCTC", exon_string)]
    return bsai_nt_list

# Return sequence surrounding the mutation, in codon list form.
def get_flanking_exon_list(codon_list, mutant_nucleotide_position, bsai_start_codon_list, bsai_end_codon_list):
    # Default values
    bsai_location_flag = "Not found"
    bsai_count_flag = 0
    bsai_shift_flag = 0
    bsai_start_codon_index = 0
    
    # Pick length of upstream/downstream sequence
    flanking_lengths = 100
    # Calculate how many codons the desired flanking length is approx. equivalent to
    flanking_codon_lengths = int((flanking_lengths - (flanking_lengths%3))/3)
    # Find which codon the mutation is in
    codon_index = extract_codon_index(mutant_nucleotide_position)
    # Get the full length of our sequence in codons
    codon_list_length = len(codon_list)
    
    # Update flanking region indexes if the mutation too close to the beginning/end of the sequence?
    upstream_codon_index, downstream_codon_index, codon_diff, location_flag = update_flanking_region_indexes(flanking_codon_lengths, codon_index, codon_list_length)
    
    # Check if Bsai codons are found inside flanking region
    for i in range(len(bsai_start_codon_list)):
        # Downstream_codon_index-1 is the TRUE index. Downstream_codon_index is technically the number of steps to move forward
        # If the Bsai sequence is found within the flanking region, can I shift the region so I keep the mutation but avoid the bsai sequence?
        # If the Bsai sequence is found between the start of the exon and the mutation codon
        if upstream_codon_index <= bsai_start_codon_list[i] < codon_index:
            bsai_count_flag = bsai_count_flag + 1
            bsai_location_flag = "BsaI sequence is upstream of the mutation"
            bsai_start_codon_index = bsai_start_codon_list[i]
            codon_index_after_bsai = bsai_end_codon_list[i] + 1
            # How many codons did I shift?
            bsai_shift_flag = codon_index_after_bsai-upstream_codon_index
            # Update the upstream codon index to begin after bsai sequence
            upstream_codon_index = codon_index_after_bsai
            # Update the downstream codon index to continue after bsai sequence IF theres enough room.
            if (downstream_codon_index + bsai_shift_flag) < codon_list_length:
                downstream_codon_index = downstream_codon_index + bsai_shift_flag
                
        elif codon_index < bsai_end_codon_list[i] <= (downstream_codon_index-1):
            bsai_location_flag = "BsaI sequence is downstream of the mutation"
            bsai_count_flag = bsai_count_flag + 1
            bsai_start_codon_index = bsai_start_codon_list[i]
            codon_index_before_bsai = bsai_start_codon_list[i]
            # How many codons did I shift?
            bsai_shift_flag = codon_index_before_bsai - downstream_codon_index
            # Update the downstream codon index to begin before bsai sequence
            downstream_codon_index = codon_index_before_bsai
            # Update the upstream codon index to continue before bsai sequence IF theres enough room.
            if (upstream_codon_index + bsai_shift_flag) >= 0:
                upstream_codon_index = upstream_codon_index + bsai_shift_flag
 
    flanking_exon_list = extract_flanking_exon_list(codon_list, upstream_codon_index, downstream_codon_index)

    return flanking_exon_list, location_flag, codon_diff, bsai_count_flag, bsai_location_flag, bsai_start_codon_index, upstream_codon_index, downstream_codon_index

# Convert codon lists to a single string
def convert_list_to_str(codon_list):
    codon_str = ''.join(codon_list)
    return codon_str

# Check ensembl's reference amino acid- does it match up with the WES data?
def compare_WES_and_ensembl_WT_aa(aa_codon_table, WT_aa_ensembl, WT_aa_WES):    
    WT_aa_match = False
    # Check the WES reference amino acid, and ensembl's reference amino acid
    if WT_aa_ensembl == WT_aa_WES:
        WT_aa_match = True
    else:
        WT_aa_match = False
    return WT_aa_match

# Check ensembl's variant amino acid- does it match up with the WES data?
def compare_WES_and_ensembl_MT_aa(mutant_nucleotide_position, MT_aa_ensembl, MT_aa_WES):
    MT_aa_match = False
    # Check the WES variant amino acid, and ensembl's variant amino acid
    if MT_aa_ensembl == MT_aa_WES:
        MT_aa_match = True
    else:
        MT_aa_match = False
    return MT_aa_match


# Check mutant nucleotide position is correct, I've only figured it out for substitutions because the data doesn't have other variants
def check_nucleotide_at_position(mutant_nucleotide_position, exon_sequence, cdna_variant):
    if ">" in cdna_variant:
        nucleotide_at_location = exon_sequence[mutant_nucleotide_position]
        reference_nucleotide = "".join(re.findall("\\d(.)>", cdna_variant))
        variant_flag = "Substitution"
        if nucleotide_at_location is not reference_nucleotide:
            variant_flag = "Reference nucleotide does not match nucleotide in ensemble sequence."
    elif "dup" in cdna_variant:
        variant_flag = "Duplication"
    elif "del" in cdna_variant:
        variant_flag = "Deletion"
    elif "ins" in cdna_variant:
        variant_flag = "Insertion"      
    else:
        variant_flag = "Error"    
    return variant_flag        



# Runs the code for getting flanking variant gene sequence 
def find_mutant_exon_seq(ensembl_gene_id, cdna_variant, mutant_nucleotide_position, aa_codon_table, WT_aa_WES, MT_aa_WES):
    # Default values:
    flanking_variant_seq = "Ensembl and WES sequences do not match"
    location_flag = ""
    flanking_variant_seq = "" 
    length_flag = -1
    codon_diff = 0
    bsai_count_flag = 0
    bsai_location_flag = 0
    bsai_codon_index = 0
    upstream_codon_index = 0
    downstream_codon_index = 0
    
    # Update nucleotide position (python being 0-based and ensembl co-ordinates being 1-based)
    mutant_nucleotide_position = int(mutant_nucleotide_position)
    mutant_nucleotide_position = mutant_nucleotide_position - 1
    
    # Get gene sequence for the reference gene
    ensembl_gene_id = remove_version_from_transcript_id(ensembl_gene_id)
    exon_sequence = get_cdna_sequence_from_transcript_id(ensembl_gene_id)
    
    # Get exon sequence
    exon_sequence = get_coding_sequence(exon_sequence)

    # Check mutant nucleotide position is correct
    variant_flag = check_nucleotide_at_position(mutant_nucleotide_position, exon_sequence, cdna_variant)
    # Get the mutant sequence in a codon list format
    mutant_exon_sequence = get_mutant_sequence(exon_sequence, cdna_variant, mutant_nucleotide_position)
    
    # protein_coding_sequence = get_exon_sequence(exon_sequence)
    reference_codon_list, CDS_begins_with_ATG_or_CTG_flag = convert_exon_sequence_into_codons(exon_sequence)
    variant_codon_list, CDS_begins_with_ATG_or_CTG_flag = convert_exon_sequence_into_codons(mutant_exon_sequence)
    
    # Extract normal and mutant codons from ensembl
    reference_codon_ensembl, mutation_codon_index = extract_codon(mutant_nucleotide_position, reference_codon_list)
    variant_codon_ensembl, mutation_codon_index = extract_codon(mutant_nucleotide_position, variant_codon_list)
    
    # Extract normal and mutant amino acids from ensembl
    WT_aa_ensembl = map_codon_to_aa(reference_codon_ensembl, aa_codon_table)
    MT_aa_ensembl = map_codon_to_aa(variant_codon_ensembl, aa_codon_table)
    
    # Compare normal and mutant ensembl amino acids to WES normal and mutant amino acids
    WT_aa_match = compare_WES_and_ensembl_WT_aa(aa_codon_table, WT_aa_ensembl, WT_aa_WES)
    MT_aa_match = compare_WES_and_ensembl_MT_aa(aa_codon_table, MT_aa_ensembl, MT_aa_WES)

    # Get flanking sequences of the mutant nucleotides only if the ensembl and WES amino acids match.       
    if ((WT_aa_match == True) & (MT_aa_match == True)):
        # Find which codons contain the BsaI sequence in our variant exon sequence
        bsai_start_codon_list, bsai_end_codon_list = find_bsai_codon_indexes(mutant_exon_sequence)
        
        # Get flanking sequence in list form
        flanking_variant_list, location_flag, codon_diff, bsai_count_flag, bsai_location_flag, bsai_codon_index, upstream_codon_index, downstream_codon_index = get_flanking_exon_list(variant_codon_list, mutant_nucleotide_position, bsai_start_codon_list, bsai_end_codon_list)

        # Convert flanking sequences to string form
        flanking_variant_seq = convert_list_to_str(flanking_variant_list)
        length_flag = len(flanking_variant_seq)
    else:
        location_flag = "Ensembl and WES sequences do not match"
        flanking_variant_seq = "Ensembl and WES sequences do not match"


    return pd.Series([variant_flag, WT_aa_ensembl, MT_aa_ensembl, WT_aa_match, MT_aa_match, 
                      flanking_variant_seq, length_flag, location_flag, codon_diff, CDS_begins_with_ATG_or_CTG_flag, bsai_count_flag, bsai_location_flag, bsai_codon_index, mutation_codon_index, upstream_codon_index, downstream_codon_index, len(mutant_exon_sequence)])

# # ASSUMMING BSAI SEQUENCE MUST BE INFRAME..
# # Check where the exon has a BsaI site: The recognition sequence is GGTCTC (a non-palindromic cutter)
# def check_bsai_sequence(exon_list):
#     bsai_codons_list = list()
#     for i in range(len(exon_list)-1):
#         # search for the recognition sequence
#         # if found, add the codon index at which it starts to the list
#         if exon_list[i] == "GGT" and exon_list[i+1] == "CTC":
#             bsai_codons_list.append(i)
#     return bsai_codons_list


# Return sequence surrounding the mutation, in codon list form.
# def get_flanking_exon_list(codon_list, mutant_nucleotide_position):
#     # Default values
#     bsai_location_flag = "Not found"
#     bsai_count_flag = 0
#     bsai_shift_flag = 0
#     bsai_codon_index = 0

    
#     # Pick length of upstream/downstream sequence
#     flanking_lengths = 100
#     # Calculate how many codons the desired flanking length is approx. equivalent to
#     flanking_codon_lengths = int((flanking_lengths - (flanking_lengths%3))/3)
#     # Find which codon the mutation is in
#     codon_index = extract_codon_index(mutant_nucleotide_position)
#     # Find which codons contain the BsaI sequence in our exon
#     bsai_codons_list = check_bsai_sequence(codon_list)
#     # Get the full length of our sequence in codons
#     codon_list_length = len(codon_list)
    
#     # Update flanking region indexes if the mutation too close to the beginning/end of the sequence?
#     upstream_codon_index, downstream_codon_index, codon_diff, location_flag = update_flanking_region_indexes(flanking_codon_lengths, codon_index, codon_list_length)
    
#     # Check if Bsai sequence is found inside flanking region
#     for i in range(len(bsai_codons_list)):
#         # Downstream_codon_index-1 is the TRUE index. Downstream_codon_index is technically the number of steps to move forward
#         # If the Bsai sequence is found within the flanking region, can I shift the region so I keep the mutation but avoid the bsai sequence?
#         # If the Bsai sequence is found between the start of the exon and the mutation codon
#         if upstream_codon_index <= bsai_codons_list[i] < codon_index:
#             bsai_count_flag = bsai_count_flag + 1
#             bsai_location_flag = "BsaI sequence is upstream of the mutation"
#             bsai_codon_index = bsai_codons_list[i]
            
#             codon_index_after_bsai = bsai_codons_list[i]+2
#             # How many codons did I shift?
#             bsai_shift_flag = codon_index_after_bsai-upstream_codon_index
#             # Update the upstream codon index to begin after bsai sequence
#             upstream_codon_index = codon_index_after_bsai
#             # Update the downstream codon index to continue after bsai sequence IF theres enough room.
#             if (downstream_codon_index + bsai_shift_flag) < codon_list_length:
#                 downstream_codon_index = downstream_codon_index + bsai_shift_flag
                
#         elif codon_index < bsai_codons_list[i] < (downstream_codon_index-1):
#             bsai_location_flag = "BsaI sequence is downstream of the mutation"
#             bsai_count_flag = bsai_count_flag + 1
#             bsai_codon_index = bsai_codons_list[i]
            
#             # How many codons did I shift?
#             bsai_shift_flag = bsai_codons_list[i] - downstream_codon_index
#             codon_index_before_bsai = bsai_codons_list[i] - 1
#             # Update the downstream codon index to begin before bsai sequence
#             downstream_codon_index = codon_index_before_bsai
#             # Update the upstream codon index to continue before bsai sequence IF theres enough room.
#             if (upstream_codon_index + bsai_shift_flag) >= 0:
#                 upstream_codon_index = upstream_codon_index + bsai_shift_flag
 
#     flanking_exon_list = extract_flanking_exon_list(codon_list, upstream_codon_index, downstream_codon_index)

#     return flanking_exon_list, location_flag, codon_diff, bsai_count_flag, bsai_location_flag, bsai_codon_index


# # Runs the code for getting flanking variant gene sequence 
# def find_mutant_exon_seq(ensembl_gene_id, cdna_variant, mutant_nucleotide_position, aa_codon_table, WT_aa_WES, MT_aa_WES):
#     # Default values:
#     flanking_variant_seq = "Ensembl and WES sequences do not match"
#     location_flag = ""
#     flanking_variant_seq = "" 
#     length_flag = -1
#     codon_diff = 0
#     bsai_count_flag = 0
#     bsai_location_flag = 0
#     bsai_codon_index = 0
    
#     # Update nucleotide position (python being 0-based and ensembl co-ordinates being 1-based)
#     mutant_nucleotide_position = int(mutant_nucleotide_position)
#     mutant_nucleotide_position = mutant_nucleotide_position - 1
    
#     # Get gene sequence for the reference gene
#     ensembl_gene_id = remove_version_from_transcript_id(ensembl_gene_id)
#     exon_sequence = get_cdna_sequence_from_transcript_id(ensembl_gene_id)
    
#     # Get exon sequence
#     exon_sequence = get_coding_sequence(exon_sequence)

#     # Check mutant nucleotide position is correct
#     variant_flag = check_nucleotide_at_position(mutant_nucleotide_position, exon_sequence, cdna_variant)
#     # Get the mutant sequence in a codon list format
#     mutant_exon_sequence = get_mutant_sequence(exon_sequence, cdna_variant, mutant_nucleotide_position)

#     # protein_coding_sequence = get_exon_sequence(exon_sequence)
#     reference_codon_list, CDS_begins_with_ATG_or_CTG_flag = convert_exon_sequence_into_codons(exon_sequence)
#     variant_codon_list, CDS_begins_with_ATG_or_CTG_flag = convert_exon_sequence_into_codons(mutant_exon_sequence)
    
#     # Extract normal and mutant codons from ensembl
#     reference_codon_ensembl, mutation_codon_index = extract_codon(mutant_nucleotide_position, reference_codon_list)
#     variant_codon_ensembl, mutation_codon_index = extract_codon(mutant_nucleotide_position, variant_codon_list)
    
#     # Extract normal and mutant amino acids from ensembl
#     WT_aa_ensembl = map_codon_to_aa(reference_codon_ensembl, aa_codon_table)
#     MT_aa_ensembl = map_codon_to_aa(variant_codon_ensembl, aa_codon_table)
    
#     # Compare normal and mutant ensembl amino acids to WES normal and mutant amino acids
#     WT_aa_match = compare_WES_and_ensembl_WT_aa(aa_codon_table, WT_aa_ensembl, WT_aa_WES)
#     MT_aa_match = compare_WES_and_ensembl_MT_aa(aa_codon_table, MT_aa_ensembl, MT_aa_WES)

#     # Get flanking sequences of the mutant nucleotides only if the ensembl and WES amino acids match.       
#     if ((WT_aa_match == True) & (MT_aa_match == True)):
#         # Get flanking sequence in list form
#         flanking_variant_list, location_flag, codon_diff, bsai_count_flag, bsai_location_flag, bsai_codon_index = get_flanking_exon_list(variant_codon_list, mutant_nucleotide_position)

#         # Convert flanking sequences to string form
#         flanking_variant_seq = convert_list_to_str(flanking_variant_list)
#         length_flag = len(flanking_variant_seq)
#     else:
#         location_flag = "Ensembl and WES sequences do not match"
#         flanking_variant_seq = "Ensembl and WES sequences do not match"


#     return pd.Series([variant_flag, WT_aa_ensembl, MT_aa_ensembl, WT_aa_match, MT_aa_match, 
#                       flanking_variant_seq, length_flag, location_flag, codon_diff, CDS_begins_with_ATG_or_CTG_flag, bsai_count_flag, bsai_location_flag, bsai_codon_index, mutation_codon_index, len(mutant_exon_sequence)])


# From pvactools' output, creates dataframe holding the transcript ID and the sequence variant
def import_pvac_data(pvac_file):
    # Call data
    neoantigen_df = pd.read_csv(pvac_file, sep="\t", usecols = columns)
    HGVSc = neoantigen_df.HGVSc.str.split(':', expand=True)[1]
    HGVSp = neoantigen_df.HGVSp.str.split(':', expand=True)[1]
    HGVSp.rename('aa_variant', inplace = True)
    HGVSc.rename('cdna_variant', inplace = True)
    neoantigen_df = pd.concat([neoantigen_df.Transcript,HGVSc,HGVSp], axis=1)
    neoantigen_df['MT_nucleotide_position'] = neoantigen_df.cdna_variant.str.extract('([0-9]+)', expand=True)
    neoantigen_df['WT_aa'] = neoantigen_df.aa_variant.str.extract('p.(.*?)[0-9]+', expand=True)
    neoantigen_df['MT_aa'] = neoantigen_df.aa_variant.str.extract('[0-9]+(.{3})', expand=True)
    return pd.DataFrame(neoantigen_df)

# Append normal AA checks, variant AA checks, flanking sequences and location flag to dataframe
def create_neoantigen_flanking_sequences(pvac_file):
    output_df = pd.read_csv(pvac_file, sep="\t")
    neoantigen_df = import_pvac_data(pvac_file)
    aa_codon_table = import_aa_codon_info()
    df_to_append = neoantigen_df.apply(lambda x: find_mutant_exon_seq(x.Transcript, x.cdna_variant, x.MT_nucleotide_position, aa_codon_table, x.WT_aa, x.MT_aa), axis=1, result_type ='expand')
    output_df = format_df(output_df, df_to_append)
    return output_df

# Reformats dataframe
def format_df(neoantigen_df, df_to_append):
    df_to_append.rename(columns={df_to_append.columns[0]: 'variant_flag',
                                  df_to_append.columns[1]: 'WT_aa_ensembl', 
                                  df_to_append.columns[2]: 'MT_aa_ensembl',
                                  df_to_append.columns[3]: 'WT_aa_match',
                                  df_to_append.columns[4]: 'MT_aa_match',
                                  df_to_append.columns[5]: 'MT_flanking_sequence',
                                  df_to_append.columns[6]: 'flanking_sequence_length',
                                  df_to_append.columns[7]: 'mutation_shift_direction',
                                  df_to_append.columns[8]: 'mutation_shift_in_codons',
                                  df_to_append.columns[9]: 'CDS_begin_with_ATG_or_CTG',
                                  df_to_append.columns[10]: 'BsaI_seq_count',
                                  df_to_append.columns[11]: 'BsaI_seq_direction',
                                  df_to_append.columns[12]: 'BsaI_codon_index',
                                  df_to_append.columns[13]: 'mutation_codon_index',
                                  df_to_append.columns[14]: 'sequence_codon_start_index',
                                  df_to_append.columns[15]: 'sequence_codon_end_index',
                                  df_to_append.columns[16]: 'CDS_length_in_codons'}, inplace=True)
    neoantigen_df = pd.concat([neoantigen_df, df_to_append], axis = 1)
    return neoantigen_df
    
#####################################################################

# Run for neoantigen file
#pvac_file = "/Users/chloetu/Desktop/mouse_epitope_extraction/formatted_data_nov.29/Lucas_MC38_Epitope_Binding_Predictions_Full_Formatted_VAF_Depth_Filtered_Multiple_Mutations_Allowed.tsv"

pvac_file = "/Users/chloetu/Desktop/mouse_epitope_extraction/formatted_data_nov.29/Correct_Lmo7.tsv"

output_df = create_neoantigen_flanking_sequences(pvac_file)


#output_df.to_excel("/Users/chloetu/Desktop/mouse_epitope_extraction/formatted_data_nov.29/Flanking_Sequences_of_Mutated_Genes_Full_VAF_Depth_Filtered_Multiple_Mutations_Allowed.xlsx", index=False)
output_df.to_excel("/Users/chloetu/Desktop/mouse_epitope_extraction/formatted_data_nov.29/Flanking_Sequences_of_correct_Lmo7.xlsx", index=False)


### To run in terminal:
    # python3 generate_mut_flanking_sequences_by_codon.py {pvac_file.tsv} {output_file_path.xlsx}
### Example:
    # python3 generate_mut_flanking_sequences_by_codon.py top_50_expressed_mutated_genes.tsv output.xlsx

# if __name__ == "__main__":
#     pvac_file = sys.argv[1]
#     output_file_path = sys.argv[2]

#     neoantigen_df = create_mut_flanking_sequences(pvac_file)
#     neoantigen_df.to_excel(output_file_path, index=FALSE)











