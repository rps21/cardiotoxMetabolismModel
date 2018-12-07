with open('proteomics_gene_lists.csv','r') as f:
    proteomics = f.readlines()

protein_list = []
for protein in proteomics:
    gene_symbol = protein.split('|')[2].strip().split('_')[0]
    protein_list.append(gene_symbol)


with open('protein_list_manualedit.txt','r') as f:
    indra_list_pre = f.readlines()

indra_list = []
for prot in indra_list_pre:
    indra_list.append(prot.strip().strip('()'))

filtered_by_proteomics = []
for prot in indra_list:
    if prot in protein_list:
        filtered_by_proteomics.append(prot)
