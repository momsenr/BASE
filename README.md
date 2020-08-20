# Background information about BASE
Repertoire analysis of patient-derived recombinant monoclonal antibodies is an important tool to study the role of B cells in autoimmune diseases of the human brain and beyond. Current protocols for generation of patient-derived recombinant monoclonal antibody libraries are time-consuming and contain repetitive steps, some of which can be assisted with the help of software automation. We developed BASE, an easy-to-use software for complete data analysis in single cell immunoglobulin cloning. BASE consists of two modules: aBASE for immunological annotations and cloning primer lookup, and cBASE for plasmid sequence identity confirmation before expression.

# How to use BASE
There are three ways to access BASE. Either you can install it locally on your machine (instructions below), which takes about 30 minutes. Alternatively, you find a fully functional and executable version of BASE which reproduces the results presented in the accompanying paper (https://www.biorxiv.org/content/10.1101/836999v1) under https://codeocean.com/capsule/3514767/tree. Duplicating this code capsule allows you to use BASE online yourself. As a third option, you can export a Docker file from the Code Ocean capsule for local installation.

# Installing BASE on your local machine
- First you need to install several dependencies:
- - install igblast (see https://ncbi.github.io/igblast/cook/How-to-set-up.html)
- - - download the germline V D and J genes (see https://ncbi.github.io/igblast/cook/How-to-set-up.html), and set up the blast databases and name them accordingly (e.g. V_without_orphon)
- - install blast
- - install the python libraries openpyxl and BioPython using pip or any other package manager
- Download a copy of this github repository
- Set up the constant_segments database for Ig subclass determination
- - download https://github.com/b-cell-immunology/sciReptor_library
- - go to directory constant_segments
- - open constant_loci.csv in a text editor and remove the lines with mouse genes from constant_loci.csv, i.e. all lines starting with "mouse" 
- - run proc_constant_loci.sh constant_loci.csv
- - move human_gl_C.fasta to the libBASE folder and change the working directory to the libBASE folder
- - run the command makeblastdb -dbtype nucl -parse_seqids -in human_gl_C.fasta -out human_gl_C/ 
- adjust libBASE/pathconfig.py to your local configuration

# BASE example files
To validate a successful installation, you can run BASE with the files in the example directory. In order to perform the example run, you need to
- copy the content from the example directory to a working directory
- run aBASE-testrun.sh or cBASE-testrun.sh (linux)
- compare your output with the expected output provided in the files aBASE-example-out.xlsx and cBASE-example-out.xlsx in the example folder

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
