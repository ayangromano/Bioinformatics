#!/usr/bin/python

from Bio import Entrez
from Bio import SeqIO

Entrez.email=''

tax_dict = {}

acc_file='top20_tig00000572_blastx.txt'

def get_accnumber(file):

    "This function grabs the accession number from a blast file. All unique accession numbers are placed in a list"


    with open (file, 'r') as acc:
        
        line = [line.split("\t")[1].split("|")[3] for line in acc]
        return list(set(line))



def esearch(accession_number):

    """esearch"""   

    search = Entrez.esearch(db="protein", term=accession_number)
    record = Entrez.read(search)
    return record['IdList'][0]



def efetch(uid):

    fetch = Entrez.efetch(db="protein", id=uid, rettype="gb", retmode="text")
    gb = SeqIO.read(fetch, 'genbank')
    if len(gb.annotations['taxonomy']) != 0:
        return gb.annotations['taxonomy'][0]+'\t'+gb.annotations['taxonomy'][1]
    else:
        return 'None'

for acc in get_accnumber(acc_file):
    search = esearch(acc)
    tax_dict[acc]=efetch(search)

with open (acc_file, 'r') as afile, open('output.txt', 'w') as tfile:
    for line in afile:
        line=line.split("\t")
        query=line[0]
        target=line[1]
        acc_num=target.split("|")[3] 
        percent_identity=line[2]
        alignment_length=line[3]
        num_mismatches=line[4]
        num_gapopens=line[5]
        start_pos_in_query=line[6]
        end_pos_in_query=line[7]
        start_pos_target=line[8]
        end_pos_target=line[9]
        eval=line[10]
        bit_score=line[11].strip()
        acc_num=target.split("|")[3]
        newline=query+'\t'+target+'\t'+percent_identity+'\t'+alignment_length+'\t'+num_mismatches+'\t'+num_gapopens+'\t'+start_pos_in_query+'\t'+end_pos_in_query+'\t'+start_pos_target+'\t'+end_pos_target+'\t'+eval+'\t'+bit_score+'\t'+tax_dict[acc_num]+'\n'
        tfile.write(newline)


"""   
for acc in get_accnumber(acc_file):
    search = esearch(acc)
    tax_dict[acc]=efetch(search)
    print acc+'\t'+tax_dict[acc]

with open(acc_file, 'r') as afile, open('manta_blastx_040316_tax_added.txt', 'w') as tfile:
    for line in afile:
        line=line.split("\t")
        query=line[0]
        target=line[1]
        percent_identity=line[2]
        alignment_length=line[3]
        num_mismatches=line[4]
        num_gapopens=line[5]
        start_pos_in_query=line[6]
        end_pos_in_query=line[7]
        start_pos_target=line[8]
        end_pos_target=line[9]
        eval=line[10]
        bit_score=line[11]
        acc_num=target.split("|")[3]
        newline=tax_dict[acc_num]+'\t'+query+'\t'+target+'\t'+percent_identity+'\t'+alignment_length+'\t'+num_mismatches+'\t'+num_gapopens+'\t'+start_pos_in_query+'\t'+end_pos_in_query+'\t'+start_pos_target+'\t'+end_pos_target+'\t'+eval+'\t'+bit_score
        print query, tax_dict[acc_num]
        
#tfile.write(newline)
"""
