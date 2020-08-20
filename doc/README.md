# Workflow
This file gives an overview of the workflow and functionality of aBASE and cBASE. For more details, please refer to the accompanying paper https://doi.org/10.1101/836999.

From amplified single cell cDNA, Ig sequences are obtained through Sanger sequencing, from which primary analysis and recommendations are provided for specific gene cloning using the aBASE (analyze-BASE) interface. Ig variable parts are then cloned into Ig expression vectors using a PCR-based approach. This cloning step can introduce nucleotide differences into the expression plasmid in comparison to the amplified single cell cDNA. Before recombinant expression of the mcAB, cBASE (compare-BASE) therefore confirms the identity of expression vector plasmid sequence and amplified cDNA, followed by the evaluation of suitability for expression. Sequence analysis algorithms for cBASE were designed to have a high positive predictive value for expression recommendation, i.e. to only recommend expression of a plasmid if it is 100% identical with the corresponding amplified cDNA. 

<img src="/doc/BASE_modules.jpg" alt="BASE modules"/>

# aBASE functionality
All relevant data acquired during mcAB generation is stored in an .xlsx table, in which each line represents one single cell with a unique mcAB identifier. This format allows the connection of sequence data derived information with any single cell metadata from other sources (e.g. cell population assignment from FACS analysis), which can simply be added in any customized extent and order. aBASE takes this .xlsx file as an input and initiates a standardized analysis algorithm for each sequence file. An output file is then automatically generated which combines the input information with the acquired data.

<img src="/doc/aBASE.jpg" alt="aBASE module"/>

# cBASE functionality

After Ig genes have been cloned into the respective expression vectors containing the constant Ig domain, plasmids are sequenced using a primer annealing upstream of the start codon. cBASE aligns and compares the plasmid Ig sequence with the amplified cDNA-derived Ig sequence by displaying nucleotide differences, assigns their relevance for amino acid changes and evaluates possible influences of primer binding interactions. Only plasmids with 100% amino acid sequence identity with the respective amplified cDNA are supposed to 
be released for expression. Results of this analysis are presented in a color encoded format, allowing rapid differentiation of plasmid genes suitable for mcAB expression (green), clones of which further bacterial colonies need to be analyzed for suitability (red) and genes where manual sequence inspection is recommended (brown). 

<img src="/doc/cBASE.jpg" alt="cBASE module"/>
