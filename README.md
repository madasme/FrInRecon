# FrInRecon
Drugs binding mode reconstruction based on fragments interactions fingerprints.

Whereas the chemical libraries contain millions of drug compounds, the vast majority of them do not have structures to crystallized targets and is therefore impossible to characterize their binding to targets from a structural view. However, the concept of building blocks offers a novel perspective on the structural problem. A drug compound is considered a complex of small chemical blocks or fragments, which confer the relevant properties to the drug and have a high proportion of functional groups involved in protein binding. Based on this, we propose a novel approach to expand the scope of structure-based repositioning approaches by transferring the structural knowledge from a fragment to a compound level.

The FrInRecon code works with the idea to take a compound for which no structural data is available and based on the binding mode of its molecular fragments, reconstruct the full compound's binding mode in terms of non-vovalent interactions fingerprints. 

**Dependencies**
* Rdkit (version)
* OpenBabel (version)

**How to get the code?**


**Input files**
1. Fragments fingerprints
2. Binding mode conservation scores for fragments

***Test case***
As a proof of concept we have tested the FrinRecon code by conducting:
1. A fragmentation process supported by PDBParser service from PharmAI GmbH and the RdKit implementations of RECAP and BRICS fragmentation algorithms. (Check pseudocode)
2. The non-covalent interactions detection by PLIP tool in plipxml format. (Check PLIP)
3. The PLIP fingerprinting supported by PharmAI GmbH services. (Check their Website)
4. Conservation Analysis of fragments binding mode using PLIP fingerprints. (Check pseudocode)

**Running the test case**


***The FrInRecon code can be used with any other fragmentation approach and any other fingerprint descriptor that can be transferred from a fragment to a compound level. Nonetheless, keep in mind that the input files format should be maintained.***







