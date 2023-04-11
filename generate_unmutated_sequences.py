#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 16:55:42 2022

@author: chloetu
"""
import sys
import pandas as pd
import re
import requests
from math import trunc


# Current script: 
#   Step 2: Non-mutated genes
# Obtains the cDNA sequence of a gene from Ensembl (the difference between 
# cDNA and genomic sequence is that the cDNA sequences do not contain introns) using transcript ID.
# Finds the exonic sequences from the cDNA sequence.
# Checks that the beginning of the coding region starts with "ATG".
# Convert exonic sequence into codon list, starting with ATG as the first codon.
# Output the total length of the gene
# Find the center of the gene, then find a 200bp region surrounding the gene.




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

# Remove version (.X) from transcript ID
def remove_version_from_transcript_id(ensembl_gene_id):
    ensembl_gene_id = ensembl_gene_id.split(".")[0]
    return ensembl_gene_id

# Convert exon sequence into a list of triplets (codons). 
def convert_cdna_sequence_into_codons(exon_sequence):
    ATG_flag = True
    CDS_begins_with_ATG_or_CTG_flag = True
    contains_exons = True
    
    codon_list = []
    match_obj = re.search("[A-Z]+", exon_sequence)
    if match_obj is None:
        codon_list = []
        contains_exons = False
    else:
        coding_sequence = match_obj.group()
        # Check if the ATG is at the start of the CDS
        atg_index = re.search("^ATG", coding_sequence)
        # Pass a flag if the coding sequence does not begin with ATG
        if atg_index is None:
            CDS_begins_with_ATG_or_CTG_flag = False
        
        for first_nt in range(0, len(coding_sequence), 3):
            codon = coding_sequence[first_nt: first_nt+3]
            if len(codon) == 3:
                codon_list.append(codon)
            
    return codon_list, CDS_begins_with_ATG_or_CTG_flag, contains_exons


# Convert codon lists to a single string
def convert_list_to_str(codon_list):
    codon_str = ''.join(codon_list)
    return codon_str


# Return sequence surrounding the mutation, in codon list form.
def get_flanking_exon_list(codon_list, middle_codon_index):
    # Pick length of upstream/downstream sequence
    flanking_lengths = 100
    # Calculate how many codons the desired flanking length is approx. equivalent to
    codon_lengths = int((flanking_lengths - (flanking_lengths%3))/3)
    # Extract X codons from above and below the center of the gene
    if (len(codon_list) < (codon_lengths*2 + 1)):
        flanking_exon_list = codon_list
    else:
        flanking_exon_list = codon_list[middle_codon_index - codon_lengths: codon_lengths + middle_codon_index + 1]
    return flanking_exon_list

# Check where the exon has a BsaI site: The recognition sequence is GGTCTC (a non-palindromic cutter)
def check_bsai_sequence(exon_sequence):
    bsai_nt_list = [index.start() for index in re.finditer("GGTCTC", exon_sequence)]
    return bsai_nt_list
    
def extract_center_sequence(ensembl_transcript_id):
    center_sequence = ""
    center_sequence_length = -1
    middle_codon_index = 0
    bsai_nt_list = []
    
    cdna_sequence = get_cdna_sequence_from_transcript_id(ensembl_transcript_id)
    codon_list, CDS_begins_with_ATG_or_CTG_flag, contains_exons = convert_cdna_sequence_into_codons(cdna_sequence)
    
    number_of_codons = len(codon_list)
    
    if CDS_begins_with_ATG_or_CTG_flag is True:
        middle_codon_index = trunc(len(codon_list)/2)
        flanking_exon_list = get_flanking_exon_list(codon_list, middle_codon_index)
        center_sequence = convert_list_to_str(flanking_exon_list)
        bsai_nt_list = check_bsai_sequence(center_sequence)
        center_sequence_length = len(center_sequence)
        
        
    return pd.Series([center_sequence, CDS_begins_with_ATG_or_CTG_flag, center_sequence_length, contains_exons, len(bsai_nt_list), middle_codon_index, number_of_codons])


# Append the exon index, max exon index, and the trimmed exon sequence
def create_normal_flanking_sequences(selections_file):
    selections_df = import_selections_data(selections_file)
    df_to_append = selections_df.apply(lambda x: extract_center_sequence(x.Name), axis = 1, result_type ='expand')
    selections_df = format_df(selections_df, df_to_append)
    return selections_df

# Reformats dataframe
def format_df(selections_df, df_to_append):
    df_to_append.rename(columns ={df_to_append.columns[0]: 'center_sequence',
                                   df_to_append.columns[1]: 'CDS_begins_with_ATG_or_CTG',
                                   df_to_append.columns[2]: 'center_sequence_length',
                                   df_to_append.columns[3]: 'contains_exons',
                                   df_to_append.columns[4]: 'BsaI_seq_count',
                                   df_to_append.columns[5]: 'middle_codon_index',
                                   df_to_append.columns[6]: 'CDS_length_in_codons'},inplace = True)
    selections_df = pd.concat([selections_df, df_to_append], axis = 1)
    return selections_df


# For GEX transcript level summarized non-mutated abundance , creates dataframe holding the transcript ID
def import_selections_data(selections_file):
    selections_df = pd.read_csv(selections_file, sep="\t", usecols = ["Name","Gene.Name","TPM"])
    return pd.DataFrame(selections_df)

    
# Run for salmon file
selections_file = "/Users/chloetu/Desktop/mouse_epitope_extraction/formatted_data_nov.29/Top_1300_Expressed_Nonmutated_Genes_VAF_Depth_Filtered_by_Transcript_TPM.tsv"


output_df = create_normal_flanking_sequences(selections_file)


output_df.to_excel("/Users/chloetu/Desktop/mouse_epitope_extraction/formatted_data_nov.29/Top_1300_Nonmutated_Center_Sequences.xlsx", index=False)


