#!/usr/bin/python
import argparse
import openpyxl
from libBASE.libBASE import SequenceFile
from libBASE.libBASE import updateExcelRow
import sys
from openpyxl.utils.cell import get_column_letter


###helper
def createExportDict(ws, begin, end):
    """returns a dict with keys='what should go in this column' and values='excel column name', created from worksheet ws
    """
    output_dict=dict()
    for row in ws[begin:end]:
        for cell in row:
            #190328: this used to be 
            #output_dict[cell.value]=cell.column
            #apparently now we have to use get_column_letter to convert the number into a letter
            output_dict[cell.value]=get_column_letter(cell.column)
    return output_dict 

#####Parse command line arguments
parser = argparse.ArgumentParser(description='aBase automatically igblasts sequencing reads and gives cloning suggestions and immunological annotations.')
parser.add_argument('input', action='store', help='input file')
parser.add_argument('output', action='store', help='output file')
parser.add_argument('--dataprefix', action='store',help='prefix for the sequence file filenames. this can be a directory name')
parser.add_argument('--hchain', action='store', help='column where the heavy chain is found')
parser.add_argument('--heavykeys', action='store', help='columns where the heavy chain output should be written to.')
parser.add_argument('--kappakeys', action='store', help='columns where the kappa chain output should be written to.')
parser.add_argument('--lambdakeys', action='store', help='columns where the lambda chain output should be written to.')
parser.add_argument('--kchain', action='store', help='column where the kappa chain is found')
parser.add_argument('--lchain', action='store', help='column where the lambda chain is found')
parser.add_argument('--overwrite', action='store_true', help='overwrite')

args = parser.parse_args()

if(args.input==args.output):
    sys.exit("input file is output file. Please do not do that. Aborting ...")

try:
    #2019-03-28: keep_vba=True was added since the file could not be saved else (see: https://bitbucket.org/openpyxl/openpyxl/issues/766/workbook-cannot-saved-twice)
    workbook = openpyxl.load_workbook(args.input, keep_vba=True)
except FileNotFoundError:
        sys.exit("File " + args.input + " not found! Aborting...")
except ValueError as my_err:
        sys.exit("OOPS! An error occured while parsing " + args.input +". " + str(my_err) +". Aborting ...")
except OSError as my_err:
        sys.exit("OOPS! An error occured while parsing " + args.input +". " + str(my_err) +". Maybe the wrong filetype? Aborting ...")

ws=workbook.active


chains={}
begin={}
end={}
columndict={}

if(args.hchain is not None):

    if(args.hchain.find(":")==-1):#we suppose this a single cell then
        args.hchain=args.hchain+":"+args.hchain

    chains['H']=workbook.active[args.hchain]
    if(args.heavykeys is not None):
        begin=args.heavykeys.split(":")[0]
        end=args.heavykeys.split(":")[1]
        columndict['H']=createExportDict(ws,begin,end)
        if('Comment' not in columndict['H'].keys()):
            sys.exit("No column titled 'Comment' in the fields specified by --heavykeys. Exiting...")

    else:
        print("--heavykeys not set. No heavy chain data will be written")

if(args.lchain is not None):
    if(args.lchain.find(":")==-1):#we suppose this a single cell then. The following line converts Z4 to Z4:Z4
        args.lchain=args.lchain+":"+args.lchain
    chains['L']=workbook.active[args.lchain]
    if(args.lambdakeys is not None):
        begin=args.lambdakeys.split(":")[0]
        end=args.lambdakeys.split(":")[1]
        columndict['L']=createExportDict(ws,begin,end)
        if('Comment' not in columndict['L'].keys()):
            sys.exit("No column titled 'Comment' in the fields specified by --lambdakeys. Exiting...")
    else:
        print("--lambdakeys not set. No lambda chain data will be written")

if(args.kchain is not None):
    if(args.kchain.find(":")==-1):#we suppose this a single cell then. The following line converts Z4 to Z4:Z4
        args.kchain=args.kchain+":"+args.kchain
    chains['K']=workbook.active[args.kchain]
    if(args.kappakeys is not None):
        begin=args.kappakeys.split(":")[0]
        end=args.kappakeys.split(":")[1]
        columndict['K']=createExportDict(ws,begin,end)
        if('Comment' not in columndict['K'].keys()):
            sys.exit("No column titled 'Comment' in the fields specified by --kappakeys. Exiting...")
    else:
        print("--kappakeys not set. No lambda chain data will be written")

parsed_sequences=[]

for ct in chains.keys():

    ###chains['H']] is loaded by chains['H']=workbook.active[args.heavy]. if args.heavy=Z (a whole row),
    ###a list is returned, but if args.heavy=Z4:Z240 (for example..), then a tuple is returen (?!?).
    ### thats why we have to iterate over seq, instead of seq###
    for seq, in chains[ct]:
        if seq.value is None:
            continue
        
        #check if this line has already been analyzed - and skip, args.overwrite is not set
        if(seq.value is not None and args.overwrite is False):
            if(ws[columndict[ct]["Confirmation"]+str(seq.row)].value!=None): 
                print(ws[columndict[ct]["Confirmation"]+str(seq.row)].value + " has already been analyzed.")
                continue

        filename=seq.value
        
        #21.03.21 the following line is a workaround for the inconsistent naming scheme of Eurofins
        filename=filename.replace("-","_")
        
        if(args.dataprefix is not None):
            filename=args.dataprefix+str(filename)+".ab1"
        try:
            my_ps=SequenceFile(filename) 
            parsed_sequences.append(my_ps)
    
            if(my_ps.successfullyParsed==False):
                ws[columndict[ct]['Comment']+str(seq.row)]=my_ps.comment 
                ws[columndict[ct]['QV']+str(seq.row)]=my_ps.mean_phred_quality
                ws[columndict[ct]['Confirmation']+str(seq.row)]="to be confirmed"
                ws[columndict[ct]['Function']+str(seq.row)]="BQ"
                ws[columndict[ct]['RL']+str(seq.row)]=my_ps.len 
            elif(my_ps.chain_type is not ct):
                ws[columndict[ct]['Comment']+str(seq.row)]=my_ps.comment+" "+filename + " has chain type " + my_ps.chain_type
                ws[columndict[ct]['QV']+str(seq.row)]=my_ps.mean_phred_quality
                ws[columndict[ct]['RL']+str(seq.row)]=my_ps.len 
                ws[columndict[ct]['Confirmation']+str(seq.row)]="to be confirmed"
                ws[columndict[ct]['Function']+str(seq.row)]="BQ"
            else:
                updateExcelRow(workbook,seq.row,columndict[ct],my_ps)
        except FileNotFoundError:
            try:
                ws[columndict[ct]['Function']+str(seq.row)]="BQ - file not found"
                ws[columndict[ct]['Comment']+str(seq.row)]="File "+str(filename)+" not found."
            except:
                sys.exit("OOPS! File " + filename + " not found! Also, either 'Comment' or 'Function' column was not found. Please make sure there are columns with that name in the range specified.")
        except ValueError as my_err:
            sys.exit("OOPS! An error occured while parsing " + filename +". " + str(my_err) +". Aborting ...")
        except OSError as my_err:
            sys.exit("OOPS! An error occured while parsing " + filename +". " + str(my_err) +". Maybe the wrong filetype? Aborting ...")

if(args.heavykeys is None or args.lambdakeys is None or args.kappakeys is None):
    sys.exit("Please set keys")

workbook.save(args.output)

