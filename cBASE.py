#!/usr/bin/python
import argparse
import openpyxl
from libBASE.libBASE import SequenceFile
from libBASE.libBASE import updateExcelRow
from libBASE.libBASE import exportDict
from Bio.Seq import Seq
import sys

from libBASE.libBASE import AlignPCRObject

from openpyxl.styles import Color, PatternFill, Font, Border
from openpyxl.styles import colors

redFill = PatternFill(start_color='FFFF0000',
                   end_color='FFFF0000',
                   fill_type='solid')
orangeFill = PatternFill(start_color='d1731e',
                   end_color='d1731e',
                   fill_type='solid')
lightgreenFill = PatternFill(start_color='92D050',
                   end_color='92D050',
                   fill_type='solid')
deepgreenFill = PatternFill(start_color='00B050',
                   end_color='00B050',
                   fill_type='solid')
purpleFill = PatternFill(start_color='826aaf',
                   end_color='826aaf',
                   fill_type='solid')
#

#####Parse command line arguments
parser = argparse.ArgumentParser(description='cBase compares the sequencing data of plasmids and PCR reads on a nucleotide per nucleotide basis.')
parser.add_argument('input', action='store', help='input file')
parser.add_argument('output', action='store', help='output file')
parser.add_argument('--dataprefix', action='store',help='prefix for the sequence file filenames. this can be a directory name')
parser.add_argument('--pcr2read', action='store', help='column where the name of the sequencing file of the 2nd pcr is found')
parser.add_argument('--plasmidread', action='store', help='column where the name of the sequencing file of the plasmid is found')
parser.add_argument('--shmanalysis', action='store', help='columns where the somatic hypermutation analysis should be written to. If four columns (separated by a comma) instead of one are given, the software will also give a more detailed analysis of the somatic hypermutations of the pcr2 read, of the plasmid, and the germline antibody.')

args = parser.parse_args()

if(args.input==args.output):
    sys.exit("input file is output file. Please do not do that. Aborting ...")

try:
    workbook = openpyxl.load_workbook(args.input)
except FileNotFoundError:
        sys.exit("File " + args.input + " not found! Aborting...")
except ValueError as my_err:
        sys.exit("OOPS! An error occured while parsing " + args.input +". " + str(my_err) +". Aborting ...")
except OSError as my_err:
        sys.exit("OOPS! An error occured while parsing " + args.input +". " + str(my_err) +". Maybe the wrong filetype? Aborting ...")

ws=workbook.active

if(args.pcr2read is not None):
    if(args.pcr2read.find(":")==-1):#we suppose this is a single cell then
        args.pcr2read=args.pcr2read+":"+args.pcr2read
    pcr2read=workbook.active[args.pcr2read]
if(args.shmanalysis is not None):
    differential_analysis_column=args.shmanalysis[0]
    if(len(args.shmanalysis)>1):
            pcr2_shm_column=args.shmanalysis.split(",")[1]
            plasmid_shm_column=args.shmanalysis.split(",")[2]
            #total_shm_column is not yet implemented as of 2019-09-04     
            total_shm_column=args.shmanalysis.split(",")[3]

to_compare=[]

###chains['H']] is loaded by chains['H']=workbook.active[args.heavy]. if args.heavy=Z (a whole row),
###a list is returned, but if args.heavy=Z4:Z240 (for example..), then a tuple is returen (?!?).
### thats why we have to iterate over seq, instead of seq###
output=""
green=False
for seq, in pcr2read:
    ###TODO: error handling!
    if seq.value is None:
        continue
    filename_pcr2=seq.value
    filename_pcr2=args.dataprefix+str(filename_pcr2)+".ab1"
    
    #the following line is a workaround for the inconsistent naming scheme of Eurofins
    filename_pcr2=filename_pcr2.replace("-","_")
    
    try:
        pcr2=SequenceFile(filename_pcr2) 
        
        filename_plasmid=ws[args.plasmidread+str(seq.row)].value
        filename_plasmid=args.dataprefix+str(filename_plasmid)+".ab1"
    
        #the following line is a workaround for the inconsistent naming scheme of Eurofins
        filename_plasmid=filename_plasmid.replace("-","_")
        
        print("Comparing: " + filename_pcr2 + " " + filename_plasmid)

        try:
            plasmid=SequenceFile(filename_plasmid)
            temp=AlignPCRObject(pcr2,plasmid)
            output=temp.output
        except FileNotFoundError as err:
            output="FileNotFound"
            print(str(err))
        except:
            output="Uncaught exception. Please inform the author of this software and provide the ab1-file(s)." 
       
        output_cell=differential_analysis_column+str(seq.row)

        ws[output_cell]=output 
        if(output.find("WARNING")!=-1):
            ws[output_cell].fill = orangeFill
        elif(output.find("nsSHM+ fr1")!=-1):
            if(output[output.find("nsSHM+ fr1")+1:].find("nsSHM+")!=-1):#if we have additional nsSHM in other regions than fr1, we mark the cell as red
                ws[output_cell].fill = redFill
            else:#else it's only orange (maybe we overlooked something about the primer?)
                ws[output_cell].fill = orangeFill
        elif(output.find("nsSHM+")!=-1):
            ws[output_cell].fill = redFill
        elif(output.find("do not match")!=-1):
            ws[output_cell].fill = redFill
        elif(output.find("index")!=-1):
            ws[output_cell].fill = redFill
        elif(output.find("Uncaught")!=-1):
            ws[output_cell].fill = redFill
        elif(output.find("FileNotFound")!=-1):
            pass
        elif(output.find("BQ")!=-1):
            ws[output_cell].fill = redFill
        elif(output.find("nsSHM-")!=-1):
            ws[output_cell].fill = orangeFill
        elif(output.find("differ")!=-1):
            ws[output_cell].fill = orangeFill
        elif(output.find("empty")!=-1):
            ws[output_cell].fill = redFill
        elif(output.find("not productive")!=-1):
            ws[output_cell].fill = purpleFill
        elif(output=="0"):
            ws[output_cell].fill = deepgreenFill
            green=True
        else:
            ws[output_cell].fill = lightgreenFill
            green=True

        if(len(args.shmanalysis)>1):
            try:
                if(temp.pcr1.successfullyParsed==True):
                    ed1=exportDict(temp.pcr1)
                    ws[pcr2_shm_column+str(seq.row)]=str(ed1['SHM'])
                else:
                    ws[pcr2_shm_column+str(seq.row)]="n/a"
            except:
                ws[pcr2_shm_column+str(seq.row)]="n/a"
            try:
                if(temp.pcr2.successfullyParsed==True):
                    ed2=exportDict(temp.pcr2)
                    ws[plasmid_shm_column+str(seq.row)]=str(ed2['SHM'])
                else:
                    ws[plasmid_shm_column+str(seq.row)]="n/a"
            except:
                ws[plasmid_shm_column+str(seq.row)]="n/a"
     
        if(green is not True):
            print(output)

        output=""
        green=False
    except FileNotFoundError as err:
        print(str(err))
    except ValueError as my_err:
            sys.exit("OOPS! An error occured while parsing " + filename_pcr2 +". " + str(my_err) +". Aborting ...")
    except OSError as my_err:
            sys.exit("OOPS! An error occured while parsing " + filename_pcr2 +". " + str(my_err) +". Maybe the wrong filetype? Aborting ...")

workbook.save(args.output)
