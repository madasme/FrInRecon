from FPrecon import *

# You will need the input files to execute the test case. You can dowload them here:
# Fragmentation file -> https://sharing.crt-dresden.de/index.php/s/weEr9nAnScvJJMM/download
# Conservation file -> https://sharing.crt-dresden.de/index.php/s/lG6v03nVYqR9lxL/download

####### INPUTS #######
report = read_pickle("frin-report_v9.pickle") #Fragmentation file (specify the right path to file)
confrag_file = "frag_target_conser_stats.csv" #Conservation file (specify the right path to file)

####### SETTINGS #######
smile = "O=c1[nH]c(=O)n(C2CC(O)C(CO)O2)cc1C=CBr" #Smiles of the ligand to reconstruct
unirep = "Q9XZT6" #Uniprot representative of the target to reconstruct
recontype = "union" #Reconstruction type. It can be union, freq or freqbin
confrag_thres = 0.2 # Threshold for the fragments conservation
proportion = 0.2 #threshold for the ligand proportion

####### EXECUTION ########
confrag = conserved_fragments(confrag_file, confrag_thres)
recon = Reconstructor(smile, unirep, report, recontype, confrag, proportion)

####### SHOWING RESULTS ########
print('Ligand smiles:',recon.smile)
print('Ligand inchikey:',recon.inchi)
print('Fragments smiles:',recon.fragsmiles)
print('Fragments inchikeys:', recon.fraginchis)
print('Ligand reconstructed fingerprint with {0}:\n{1}'.format(recontype, recon.rec_fp))
