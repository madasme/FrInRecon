#__Author__="Melissa F. Adasme"
#Generates the consensus fingerprint for a given fragment by merging all fingerprints in datasets

########################
### IMPORTED MODULES ###
########################

import pandas as pd
import pickle
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import Recap, AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') #Disable warnings
import numpy as np


#######################################
####### SUPPLEMENTARY FUNCTIONS #######
#######################################

def conserved_fragments(path, conser_threshold):
    """Read Tanimoto similarity data and return list of conserved fragments."""
    tis_data = pd.read_csv(path, sep="\t")
    tis_data = tis_data[tis_data["FPSAVG"] >= conser_threshold]
    return tis_data["FRAGMENT"].values.tolist()

def obconvert(pybelmolecule, outputformat='can', ignoreStereo=True, stripH=False):
    """Returns the canonical SMILES of a molecule without stereo info. Hydrogens are stripped upon request."""
    obc = pybel.ob.OBConversion()
    obc.SetOutFormat(outputformat)
    if ignoreStereo:
        obc.AddOption('i')  # ignore stereo info
    if stripH:
        obc.AddOption('DeletePolarH')
        obc.AddOption('DeleteNonPolarH')
    obc.AddOption('n')  # don't write the name to the output SMILES
    return obc.WriteString(pybelmolecule.OBMol).strip()

def inchikey(pybelmolecule, ignoreStereo=True, stripH=False, use_rdkit=True):
    """Returns the InChI key of a molecule without stereo info. Hydrogens are stripped upon request."""
    if not use_rdkit:
        return obconvert(pybelmolecule, outputformat='inchikey', ignoreStereo=ignoreStereo, stripH=stripH)
    else:
        smiles = obconvert(pybelmolecule, outputformat='can', ignoreStereo=ignoreStereo, stripH=stripH)
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        try:
            inchikey = Chem.inchi.InchiToInchiKey(Chem.MolToInchi(m))
        except ValueError:  # E.g. Sanitization error
            return None
        return inchikey

def remove_dummy_atom(fragsmile):
    """Cleans the fragment mol removing the * atoms"""
    if '*' in fragsmile:
        mol = Chem.MolFromSmiles(fragsmile, sanitize=False) #rdkit mol from fragment smiles
        #print (Chem.MolToSmiles(mol))
        duatom = Chem.MolFromSmiles('*') #dummy atom saved as molecule
        mol = AllChem.ReplaceSubstructs(mol,duatom,Chem.MolFromSmiles('[H]'),True)[0] #replace * by H)
        #print (Chem.MolToSmiles(mol))
        mol = Chem.RemoveHs(mol) #delete Hydrogens
        #print (Chem.MolToSmiles(mol))
        fragsmile = Chem.MolToSmiles(mol) #rdkit smiles coming from clean mol
    return fragsmile

def read_pickle(path):
    """read a pickle file and returns it"""
    pickle_file = open(path,"rb")
    data = pickle.load(pickle_file)
    pickle_file.close()
    return data

def conserved_fragments(path, threshold):
    """Read Tanimoto similarity data and return list of conserved fragments."""
    tis_data = pd.read_csv(path, sep="\t")
    tis_data = tis_data[tis_data["FPSAVG"] >= threshold]
    return tis_data["FRAGMENT"].values.tolist()




########################
##### MAIN CLASSES #####
########################

class FPrepre():
    """The class generates the representative fingerprint for a fragment with the given set of fingerprints coming from different complexes"""
    def __init__(self, liginchi, inchi, unirep, report):
        self.liginchi = liginchi
        self.inchi = inchi
        self.report = report
        self.FPs = self.get_fps(self.inchi, unirep=unirep)
        if self.FPs:
            self.FP_union = self.get_fp_union(self.FPs)
            self.FP_freq = self.get_fp_freq(self.FPs)
            self.FP_freq_bin = self.get_fp_freq_binary(self.FPs)
        else:
            self.FP_union = None
            self.FP_freq = None
            self.FP_freq_bin = None

    def get_fps(self, inchi, unirep=None):
        """Given a fragment inchikey and a possible uniprot representative, obtains all fingerprints of fragment complexes from the report"""
        if unirep: #If unirep provided, then the subsetting is based in fraginchi and uniprot rep. If not, then just in fraginchi
            complexes = self.report[(self.report["FRAGINCHI"] == inchi) & (self.report['UNIREP'] == unirep)]
        else:
            complexes = self.report[(self.report["FRAGINCHI"] == inchi)]
        complexes = complexes[complexes.LIGINCHI != self.liginchi] #Remove the complex with the drug to reconstruct, to avoid overfitting
        fps = complexes.FPSIMPLE.tolist()
        if len(fps) == 0: #If no complexes for given fragment then set to None
            fps = None
        return fps

    def get_fp_union(self, fps, string=False):
        """Given a list of fingerprints it constructs a representative one of the union: for a given position, if at least in one fp the bin is on, then in the representative it is on"""
        if len(fps) == 1: #If only one, then that one becomes the representative
            rep_FP = fps[0]
        else:
            all_FPs = zip(*fps)#Generates a tuple of bins in the same index but coming from different FPs lists:
            rep_FP = [1 if sum(map(int,item))>0 else 0 for item in all_FPs]#Iterates over all tuples. Transform each tuple into integers and sums each tuple. If the sum is bigger than 0 then adds a 1 into the consensus FP, otherwise a 0.
        return rep_FP

    def get_fp_freq(self,fps, string=False):
        """Givel a list of fingerprints it constructs a representative one of the frequency"""
        lenth = len(fps)
        if lenth == 1: #If only one, then that one becomes the representative
            rep_FP = fps[0]
        else:
            all_FPs = zip(*fps)  # Generates a tuple of bins in the same index but coming from different FPs lists:
            rep_FP = [sum(map(int,item))/lenth for item in all_FPs] #Iterates over all tuples. Transform each tuple into integers and sums the numbers in tuple. Round the result to one decimal
        return rep_FP

    def get_fp_freq_binary(self, fps, string=False):
        """Givel a list of fingerprints it constructs a representative one of the frequency. According to cutoff (>=30%=1) converts the fp into a binary"""
        lenth = len(fps)
        if lenth == 1: #If only one, then that one becomes the representative
            rep_FP = fps[0]
            return rep_FP
        else:
            all_FPs = zip(*fps)#Generates a tuple of bins in the same index but coming from different FPs lists:
            freq_FP = [round(sum(map(int,item))/lenth,1) for item in all_FPs] #Iterates over all tuples. Transform each tuple into integers and sums the numbers in tuple. Round the result to one decimal
            freq_FP_bin = [1 if per >= 0.3 else 0 for per in freq_FP] #If the frequency is 0.3 (30%) or more, then adds a 1 into the consensus FP, otherwise a 0
            return freq_FP_bin


class Reconstructor():
    """The class generates the reconstructed fingerprint for a given ligand smile."""
    def __init__(self, smile, unirep, report, recontype, confrag, threshold):
        self.smile = smile
        self.unirep = unirep
        self.report = report
        self.recontype = recontype
        self.inchi = self.get_inchikey(self.smile)
        self.confrag = confrag
        self.threshold = threshold # Proportion of ligand fragments that have an interaction fingerprint.
        self.fragsmiles = self.get_fragments(self.smile, app='recap')
        self.used_fragments = []
        if self.fragsmiles:
            self.fraginchis = [self.get_inchikey(smile) for smile in self.fragsmiles]
            self. rec_fp = self.get_rec_fp()
        else:
            self.fraginchis = None
            self.fr_fp_union = None
            self.fr_fp_freq = None
            self.fr_fp_freqbin = None
            self.rec_fp = None

    def get_inchikey(self, smile):
        """Given a ligand or fragment smiles, generates the pybel molecule and then obtains the inchikey """
        try:
            smile = remove_dummy_atom(smile)  #from supple functions on top
            pybel_frag_mol = pybel.readstring('can', smile) #gets pybel query mol of fragment coming from fragsmiles
            pybel_frag_mol.removeh()
            inchi = inchikey(pybel_frag_mol) #from supple functions on top
        except:
            inchi = None
        return inchi

    def get_fragments(self,smile, app):  # Add option for BRICKS
        """Given a ligand smiles, returns a list of canonical SMILES fragments using RECAP algorithm. If not possible it returns None"""
        fragsmiles = None
        mol = Chem.MolFromSmiles(smile)
        if mol:
            hierarch = Recap.RecapDecompose(mol)
            recap_fragments = hierarch.GetLeaves()  # Gets the last childs in hierarch (the ligand splitted simply and without redundant data)
            if len(list(recap_fragments.keys())) > 0:  # If RECAP generated at least one fragment
                fragsmiles = list(recap_fragments.keys())
        return fragsmiles

    def get_rec_fp(self):
        """For each fragment in the ligand generates an instance of the FPmerge class and obtains the representative fingerprint according to the 3 types: union, frequency and frequency as binary"""
        self.fr_fp_union = []  # Initialized as list
        self.fr_fp_freq = []
        self.fr_fp_freqbin = []
        rec_fp = None  # Initialized as none
        # For each fragment smile generates an instance of FPrepre class to get the representative fingerprint
        for fragsmile in self.fragsmiles:
            fraginchi = self.get_inchikey(fragsmile)
            if not "{}~{}".format(fraginchi, self.unirep) in self.confrag:
                continue
            self.used_fragments.append(fragsmile) # Here in case the threshold is modified and allows to reconstruct drugs with not all fragments
            fprep = FPrepre(self.inchi, fraginchi, self.unirep, self.report)
            if fprep.FP_union:  # If fragment had FPs and generated the representative, then it is added to the list. If it is None, then not.
                self.fr_fp_union.append(fprep.FP_union)
            if fprep.FP_freq:
                self.fr_fp_freq.append(fprep.FP_freq)
            if fprep.FP_freq_bin:
                self.fr_fp_freqbin.append(fprep.FP_freq_bin)
        # For each type of representative fragments a different reconstruction of the full FP for ligand is needed
        if self.recontype == 'union':
            if self.fr_fp_union and (len(self.fr_fp_union)/len(self.fragsmiles)) >= self.threshold:
                rec_fp = self.fr_fp_union[0]  # By default set to the first fragment fp in the list
                if len(self.fr_fp_union) > 1:  # If there is more than just one FP, then takes the union of them
                    zipped_fps = zip(*self.fr_fp_union)
                    rec_fp = [1 if sum(list(map(int, item))) > 0 else 0 for item in zipped_fps]
        if self.recontype == 'freq':
            if self.fr_fp_freq and (len(self.fr_fp_freq)/len(self.fragsmiles)) >= self.threshold:
                rec_fp = [round(float(bin), 1) for bin in self.fr_fp_union[0]]  # By default set to the first fragment fp in the list
                if len(self.fr_fp_freq) > 1:  # If there is more than just one FP, then takes the union of them
                    zipped_fps = zip(*self.fr_fp_freq)
                    rec_fp = [round(np.mean(list(map(float, item))), 1) for item in zipped_fps]  # If two fragments have the same bin acticated, the chances increase
        if self.recontype == 'freqbin':
            if self.fr_fp_freqbin and (len(self.fr_fp_freqbin)/len(self.fragsmiles)) >= self.threshold:
                rec_fp = self.fr_fp_freqbin[0]  # By default set to the first fragment fp in the list
                if len(self.fr_fp_freqbin) > 1:  # If there is more than just one FP, then takes the union of them
                    zipped_fps = zip(*self.fr_fp_freqbin)
                    rec_fp = [1 if sum(list(map(int, item))) > 0 else 0 for item in zipped_fps]  # If at least one activated then set to 1
        #print(rec_fp)
        return rec_fp
