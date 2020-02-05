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
