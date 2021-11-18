# FrInRecon
Drugs binding mode reconstruction based on fragments interactions fingerprints.

Whereas the chemical libraries contain millions of drug compounds, the vast majority of them do not have structures to crystallized targets and is therefore impossible to characterize their binding to targets from a structural view. However, the concept of building blocks offers a novel perspective on the structural problem. A drug compound is considered a complex of small chemical blocks or fragments, which confer the relevant properties to the drug and have a high proportion of functional groups involved in protein binding. Based on this, we propose a novel approach to expand the scope of structure-based repositioning approaches by transferring the structural knowledge from a fragment to a compound level.

The FrInRecon code works with the idea to take a compound for which no structural data is available and based on the binding mode of its molecular fragments, reconstruct the full compound's binding mode in terms of non-vovalent interactions fingerprints. The reconstruction is target-based, meaning that the binding mode of fragments is defined in reference to a given target, therefore, the reconstructed binding mode of the full compound will be for the same target.

## Usage 
#### Dependencies
* Rdkit (version)
* OpenBabel (version >= 3.0.0) 
* NumPy (version )
* Pickle (version )
* Pandas (version)

OpenBabel: Many users encounter trouble setting up OpenBabel with Python bindings correctly. We therefore provide some [installation help for OpenBabel](#ob) below.

#### Clone the repository  
Open a terminal and clone this repository using
```
git clone https://github.com/madasme/FrInRecon.git
```

#### Setting the PYTHONPATH  
In your terminal, add the FrInRecon repository to your ```PYTHONPATH``` variable. You can also do it manually in your ./bashrc file.  
```
$ export PYTHONPATH=~/FrInRecon:${PYTHONPATH}
```

#### Input files
1. Fragments fingerprints
2. Binding mode conservation scores for fragments

#### Required settings
1. The canonical smiles of the compound to reconstruct
2. The Uniprot/Unirep ID of the target 
3. The reconstruction type to be used (e.g. union, frequency, or binary frequency)
4. Threshold of fragments binding mode conservation (i.e. how corserved has to be the binding mode in terms of mTIS)
5. Threshold of ligand proportion (i.e. proprtion of ligand fragments to be found in dataset)

### Test case
As a proof of concept we have tested the FrinRecon code by conducting:
1. A fragmentation process supported by PDBParser service from PharmAI GmbH and the RdKit implementations of RECAP and BRICS fragmentation algorithms. (Check pseudocode)
2. The non-covalent interactions detection by PLIP tool in plipxml format. ([Check PLIP](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index))
3. The PLIP fingerprinting and isomorphism method provided by PharmAI GmbH services. ([Check their Website](https://www.pharm.ai/))
4. Conservation Analysis of fragments binding mode using PLIP fingerprints. (Check pseudocode)

**Test case code example**
```
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
```


#### Running the test case
```
$ python3 testcase.py 
```

***The FrInRecon code can be used with any other fragmentation approach and any other fingerprint descriptor that can be transferred from a fragment to a compound level. Nonetheless, keep in mind that the input files format should be maintained.***


### Installing Open Babel <a name="ob"></a>
As many users encounter problems with installing the required OpenBabel tools, we want to provide some help here. However, we cannot offer technical support. Comprehensive information about the installation of OpenBabel for Windows, Linux, and macOS can be found in the [OpenBabel wiki](http://openbabel.org/wiki/Category:Installation) and the [OpenBabel Docs](https://open-babel.readthedocs.io/en/latest/Installation/install.html). Information about the installation of [OpenBabel Python bindings](https://open-babel.readthedocs.io/en/latest/UseTheLibrary/PythonInstall.html) can also be found there.


#### Using Conda, HomeBrew or the binary Windows Installer
Install OpenBabel using the binary from GitHub or with
```
# For Conda users
$ conda install openbabel -c conda-forge
# On macOS
$ brew install open-babel
```
Install the Python bindings with:  
```
$ pip install openbabel
```

* Note: If you have trouble, make sure the OpenBabel version matches the one for the python bindings!

#### Using your Package Manager (Example for Ubuntu 20.04)
```
$ apt-get update && apt-get install -y \
    libopenbabel-dev \
    libopenbabel6 \
    python3-openbabel \
    openbabel
```

#### From Source (Example for Ubuntu 18.04)
Clone the OpenBabel repository into the /src directory
```
$ git clone -b openbabel-3-0-0 \
https://github.com/openbabel/openbabel.git
```
Within /src/openbabel create end enter a directory /build and configure the build using
```
$ cmake .. \
-DPYTHON_EXECUTABLE=/usr/bin/python3.6 \
-DPYTHON_BINDINGS=ON \
-DCMAKE_INSTALL_PREFIX=/usr/local \
-DRUN_SWIG=ON
```
From within the same directory (/src/openbabel/build) compile and install using
```
$ make -j$(nproc --all) install
```
