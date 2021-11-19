from FPrecon import *

####### INPUTS #######
report = read_pickle("frin-report_v9.pickle")
confrag_file = "frag_target_conser_stats.csv"

####### SETTINGS #######
smile = "O=c1[nH]c(=O)n(C2CC(O)C(CO)O2)cc1C=CBr" #Smiles of the ligand
unirep = "Q9XZT6" #Uniprot representative of the target
recontype = "union" #it can be union, freq or freqbin
confrag_thres = 0.2
proportion = 0.2

####### EXECUTION ########
confrag = conserved_fragments(confrag_file, confrag_thres)
recon = Reconstructor(smile, unirep, report, recontype, confrag, proportion)

####### SHOWING RESULTS ########
print('Ligand smiles:',recon.smile)
print('Ligand inchikey:',recon.inchi)
print('Fragments smiles:',recon.fragsmiles)
print('Fragments inchikeys:', recon.fraginchis)
print('Ligand reconstructed fingerprint with {0}:\n{1}'.format(recontype, recon.rec_fp))
