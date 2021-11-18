from FPrepre import *

####### INPUTS #######
smile = "O=c1[nH]c(=O)n(C2CC(O)C(CO)O2)cc1C=CBr" #Smiles of the ligand
unirep = "Q9XZT6" #Uniprot representative of the target
report = read_pickle("frin-report_v9.pickle") #Reads the full report to extract the fragments fingerprints
confrag = conserved_fragments("frag_target_conser_stats.csv", 0)
proportion = 0

####### EXECUTION ########
recontype = "union" #could be union, freq or freqbin
recon = Reconstructor(smile, unirep, report, recontype, confrag, proportion)

####### PRINTING RESULTS ########
print('Ligand smiles:',recon.smile)
print('Ligand inchikey:',recon.inchi)
print('Fragments smiles:',recon.fragsmiles)
print('Fragments inchikeys:', recon.fraginchis)
print('Ligand reconstructed fingerprint with {0}:\n{1}'.format(recontype, recon.rec_fp))
