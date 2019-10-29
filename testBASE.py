#!/usr/bin/python
from Bio import SeqIO
import argparse
from libBASE.libBASE import SequenceFile
from libBASE.libBASE import exportDict
import sys


#####Parse command line arguments
parser = argparse.ArgumentParser(description='Load an *.ab1 sequence file and process the data.')
parser.add_argument('input', action='store',  metavar='filename(s)', nargs='+', help='name of the sequence file(s) to parse')
parser.add_argument('--quality', '-q', action='store_true', help='give parsing quality measures (debugging only)')
parser.add_argument('--export', '-e', action='store', help='store fasta sequence in the given output file')
parser.add_argument('--debug', '-j', action='store_true', help='give output)')
parser.add_argument('--igblast', action='store', help='save igblast output to file')

args = parser.parse_args()
parsed_sequences=[]
ed=[]
for filename in args.input:
    try:
        parsed_sequences.append(SequenceFile(filename,"abi",args.igblast))
    except FileNotFoundError:
        sys.exit("OOPS! File " + filename + " not found! Aborting...")
    except ValueError as my_err:
        sys.exit("OOPS! An error occured while parsing " + filename +". " + str(my_err) +". Aborting ...")
    except OSError as my_err:
        sys.exit("OOPS! An error occured while parsing " + filename +". " + str(my_err) +". Maybe the wrong filetype? Aborting ...")


if (args.quality):
    for seq in parsed_sequences:
        if(hasattr(seq, 'chain_type') and seq.chain_type is "H"):
            print(seq.comment)
            print(seq.filename + ". IgSC identification yields: " + seq.IgSubClass)

            
if(args.debug):
    for ps in parsed_sequences:
        if(hasattr(ps,'BlastedOutputDict')):
            print(ps.BlastedOutputDict)
        ed=exportDict(ps)
        for k in ed:
            print( str(k) + str(ed[k]))
    
if(args.export is not None):
    if(len(parsed_sequences)==1):
        parsed_sequences[0].writeToFasta(args.export)
        #print(str(parsed_sequences[0].seq)[::-1])
    else:
        for seq in parsed_sequences:
            seq.writeToFasta(seq.filename + "_" + args.export)


