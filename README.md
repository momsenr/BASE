# Setting up BASE
- install igblast (see https://ncbi.github.io/igblast/cook/How-to-set-up.html)
- - download the germline V D and J genes (see https://ncbi.github.io/igblast/cook/How-to-set-up.html), and set up the blast databases and name them accordingly (e.g. V_without_orphon)
- install blast
- install the python libraries openpyxl and BioPython using pip or any other package manager
- adjust libBASE/pathconfig.py to your local configuration

# Setting up the constant_segments database for Ig subclass determination
- download https://github.com/b-cell-immunology/sciReptor_library
- go to directory constant_segments
- remove the lines with mouse genes from constant_loci.csv
- run proc_constant_loci.sh constant_loci.csv
- move human_gl_C.fasta to the libIgParser folder
- run the command makeblastdb -dbtype nucl -parse_seqids -in human_gl_C.fasta -out human_gl_C/ 

# Run with example files
- copy example files from example directory to working directory
- run aBASE-testrun.sh or cBASE-testrun.sh (linux)

# Explanation of cBASE output abbreviations

| Abbreviation      | Explanation                                                                                                                                               | Consequence for expression                                                                    |
|-------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| BQ                | Bad sequencing quality                                                                                                                                    | Repeat sequencing                                                                             |
| x SHM- FR1(P)     | x somatic hypermutations less in the FR1 region of the expression plasmid (this is due to primer artifacts in the PCR product).                           | This plasmid can be used for expression.                                                      |
| x SHM- J(P)       | x somatic hypermutations less in the J region of the expression plasmid (this is due to primer artifacts in the PCR product).                             | This plasmid can be used for expression.                                                      |
| x SHM J(P)+       | x additional somatic hypermutations in the J region of the expression plasmid (this is due to primer artifacts in the PCR product).                       | This plasmid can be used for expression.                                                      |
| x sSHM REGION     | x additional silent somatic hypermutations in the expression plasmid in region REGION. This has been artificially introduced in the cloning procedure.    | As the mutations have no effect on the protein level, the plasmid can be used for expression. |
| x nsSHM REGION    | x additional non-silent somatichypermutations in the expression plasmid in region REGION. This has been artificially introduced in the cloning procedure. | The expression plasmid should be discarded.                                                   |
| x sSHM- REGION    | x non-silent somatic hypermutations less in the expression plasmid in region REGION. This has been artificially introduced in the cloning procedure.      | As the mutations have no effect on the protein level, the plasmid can be used for expression. |
| x nsSHMchg REGION | x changed non-silent hypermutations in the expression plasmid in region REGION. This has been artificially introduced in the cloning procedure.           | The expression plasmid should be discarded.                                                   |
| x nsSHM- REGION   | x non-silent hypermutations less in the expression plasmid in region REGION. This has been artificially introduced in the cloning procedure.              | The expression plasmid should be discarded.                                                   |
| 0                 | No differences between expression plasmid and PCR product on the nucleotide level                                                                         | Congratulations!                                                                              |
