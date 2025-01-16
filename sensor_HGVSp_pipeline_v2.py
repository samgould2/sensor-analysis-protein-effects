#importing necessary packages
from collections import defaultdict
import numpy as np
import Bio.Seq
import pandas as pd
import os
import sys
from pathlib import Path
#package for hamming distance 
import Levenshtein as lv


#-------Loading in system args
# Parse command-line arguments to read in files of interest
if len(sys.argv) != 5:
    print("Usage: python3 sensor_HGVSp_pipeline.py <master_df> <cdks_df> <EDITOR> <sample_name>")
    sys.exit(1)

master_df = pd.read_csv(Path(sys.argv[1]))
cdks_df = pd.read_csv(Path(sys.argv[2]))
EDITOR = str(sys.argv[3])
sample_name = str(sys.argv[4])


#------FUNCTIONS--------
def hgvsp_simple(wt_seq, mut_seq):
    if wt_seq==mut_seq:
        return 'WT'
    else:
        pos_mutated = []

        for i in range(len(wt_seq)):
            if wt_seq[i] != mut_seq[i]:
                pos_mutated.append(i + 1)
        
        hg = ''
        for idx, pos in enumerate(pos_mutated):
            if idx==0:
                hg+= f'{wt_seq[pos-1]}{pos}{mut_seq[pos-1]}'
            else:
                hg+= f'_{wt_seq[pos-1]}{pos}{mut_seq[pos-1]}'
        

        return hg


def edit_classifier(a, gRNA_id, cdks, master_df, EDITOR):
    """ 
    a = allele_frequency_table.zip
    gRNA_id = gRNA_id
    cdks = cdk information (with WT tx and protein)
    master_df = information about location of protospacer relative to cDNA (for HGVSp determination)
    EDITOR = type of editor ('ABE' or 'CBE') for determination of if edit is canonical or not
    """

    #reverse complement because sensor is in rev-comp orientation
    edits = [str(Bio.Seq.Seq(i.replace('-', '')[12:32]).reverse_complement()) for i in a['Aligned_Sequence']]
    reads = list(a['#Reads'])
    edit_table = pd.DataFrame(dict(zip(['Edit', '#Reads'], [edits, reads]))).groupby('Edit').sum().sort_values(by='#Reads', ascending=False).reset_index()
    subset = master_df[master_df['gRNA_id']==gRNA_id]

    wt_protospacer = subset['protospacer'].values[0]
    proto_rc = str(Bio.Seq.Seq(wt_protospacer).reverse_complement())

    orientation = subset['orientation'].values[0]
    p_start = subset['prot_start'].values[0]
    p_end = subset['prot_end'].values[0]
    cdna_start = subset['cdna_start'].values[0]
    cdna_end = subset['cdna_end'].values[0]

    gene = subset['Gene'].values[0]
    wt_tx = cdks.loc[cdks['Gene']==gene, 'WT_tx_full'].values[0]
    wt_prot = cdks.loc[cdks['Gene']==gene, 'Protein'].values[0]

    hgs = []
    hamming = []
    edit_classification = []
    canonical = []
    window = []
    for i, val in edit_table.iterrows():
        read = val['Edit']

        #orientation refers to orientation of the read relative to the cDNA sequence
        #for = same orientation
        #rev = reverse complement
        #all reads modified to be in protospacer orientation
        if orientation == 'for':
            m = read[p_start:p_end]
            mut_seq = wt_tx[:cdna_start] + m + wt_tx[cdna_end:]
            if len(wt_protospacer)==len(m):
                hd = lv.hamming(wt_protospacer, m)
            else:
                hd = 20

        elif orientation == 'rev':
            m = read[p_start:p_end]
            m2 = str(Bio.Seq.Seq(m).reverse_complement())
            mut_seq = wt_tx[:cdna_start] + m2 + wt_tx[cdna_end:]
            if len(proto_rc) ==len(m2): #edge case of short reads
                hd = lv.hamming(proto_rc, m2)
            else:
                hd = 20

        mut_prot = str(Bio.Seq.Seq(mut_seq).transcribe().translate())

        j = hgvsp_simple(wt_prot, mut_prot)
        hgs.append(j)
        hamming.append(hd)

        #and then classification of the edits
        count = 0
        edit_classif = ''
        canonical_edit = True
        canonical_window = True

        for k, val2 in enumerate(read):
            wt_base = wt_protospacer[k]
            proto_location = int(k+1)
            if val2 != wt_base:
                count+=1
                edit_classif+=f'+{proto_location}{wt_base}>{val2},'

                if EDITOR == 'CBE':
                    if f'{wt_base}>{val2}' != 'C>T':
                        canonical_edit = False
                elif EDITOR =='ABE':
                    if f'{wt_base}>{val2}' != 'A>G':
                        canonical_edit = False

                if (proto_location<4) or (proto_location>8):
                    canonical_window=False

        if count ==0:
            edit_classif = 'No edit'


        edit_classification.append(edit_classif)
        canonical.append(canonical_edit)
        window.append(canonical_window)


    edit_table['HGVSp'] = hgs
    edit_table['Num_edits'] = hamming
    edit_table['DNA Change'] = edit_classification
    edit_table['Canonical_edit'] = canonical
    edit_table['Canonical_window'] = window
    edit_table['gRNA_id'] = gRNA_id
    
    return edit_table


#--------and then iterate over the mastertable----------

fp = './crispresso'

rows = []
for i, val in master_df.iterrows():

    gRNA_id = val['gRNA_id']

    #unique_id = val['unique_id']
    unique_id = gRNA_id
    
    output_folder_x = os.listdir(fp + '/' + sample_name + f"/CRISPResso_on_{unique_id}")


    if 'Alleles_frequency_table.zip' in output_folder_x:

        a = pd.read_csv(fp + '/' + sample_name + f"/CRISPResso_on_{unique_id}/Alleles_frequency_table.zip", sep='\t')

        edit_table = edit_classifier(a, gRNA_id, cdks_df, master_df, EDITOR)

        rows.append(edit_table)

    else:
        #in the event there's no sensor reads, just ignore this (can deal with it later on...)
        continue
    
#and then concatenate and output it
combined = pd.concat(rows).reset_index(drop=True)
combined.to_csv(f"./crispresso/{sample_name}_HGVSp_sensor_quant_v2.csv", index=False)
