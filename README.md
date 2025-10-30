# Medicinal Compound Dataset Generation
## SMILES Extraction
In order to search for compound similarity, the SMILES (Simplified Molecular Input Line Entry System) had to be extracted from the .SDF files containing the three different reference sample. 
Tool used for extraction: RDKIT Package in Python
## Similarity Generation
Using the Tanimoto similarity score of 70%, similar structures were then extracted from various databases for each of the reference structure. For further purification, duplicate results from two or more databases were filtered out, and the final number for all three reference compounds was 3000.
Databases used: PubChem, Zinc, ChemBl, DrugBank, and Cactus NIH. 
Tools: Python and Manual Search.
Similarities found in total for each reference compound: 1000 each. 
Total = 3000
