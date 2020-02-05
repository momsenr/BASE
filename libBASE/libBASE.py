from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import shutil
import os
import platform
import tempfile
import collections
import sys
import subprocess

import libBASE.pathconfig as cfg
from libBASE.primer import primer_IGHV, primer_IGHJ, primer_IGKV, primer_IGKJ, primer_IGLV, primer_IGLJ
from libBASE.IgBlastParser import LoadBlastedOutput

#This is to silence the warnings aboutSequences with a number of nucleotides not a multiple of 3
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)


#When comparing pcr2 and plasmid reads in cBASE/AlignPCRObject, we have a certain tolerance in the beginning and end of 
#the pcr2 sequencing read, depending on the length of the primer.
nt_tolerance_fwd={'H': 33, 'K': 32, 'L': 27} #pragmatically adjusted the correct values by +9 for H and L, +11 for K 
nt_tolerance_rev={'H': 21, 'K': 21, 'L': 0} 

#Being able to check for certain constant empty vector sequences can help troubleshooting when we habe low expression 
#cloning efficiency (needed for cBASE)
empty_vec_seq={'H': 'ACCGGTGTACACTCGAGCGTACGGTCGAC', 'K': 'ACCGGTTGACTAACTAGCCGTACG', 'L': 'ACCGGTTGACTAACTAGCCTCGAG'}

#When we did not use Gibson assembly for expression cloning, checking for certain restriction motives was necessary.
#We leave these checks for backwards-compatability (needed for aBASE)
restriction_motif ={'AgeI':'ACCGGT',
                'BsiWI':'CGTACG',
                'XhoI':'CTCGAG',
                'SalI':"GTCGAC"
                }

class SequenceFile():
        '''
        Pass the whole sequence file to this class and it will deal with it
        Arguments:
            filename - string for the sequence file to parse

        Optional arguments:
            filetype - string for the filetyp of filename, default is abi.
            igblast_output - where to write the output of igblastn to. Defaults to none. 

        Methods:
            writeToFast: Write file to fasta, takes output filename as an argument
            IgBlastMe: blasts the sequence and determes the Ig Subclass. takes a filename as optional parameter to write the output of igblastn to 
            identify_gene_region: takes parameter pos (position relative to query sequence), returns the gene region (possible return values:
                                    v,d,j, v_d_region, d_j_junction, v_j_junction)
            identify_V_gene_subregion takes parameter pos (position relative to query sequence), returns the V gene subregion (possible return values:
                                    fr1,cdr1,fr2,cdr2,fr3,cdr3,n/d)
            original_nt: takes parameter pos (position relative to V gene start), returns the nt at pos of the aligned germline gene sequence
            IsThisNTSilent: takes parameter pos relative to the beginning of the V gene, returns True, if the aligned nt is silent (ie either is present 
                            in the aligned gene, or the nucleotide change does not lead to an AA change. 
            translatedAA: takes parameter pos (position relative to the beginning of the V gene), returns the translation of the aligned sequence and the gene sequence
        '''
        def __init__(self, filename,filetype="abi",igblast_output=None):
            """
            parses the file into self.record, then the IgSubClass is analyzed and the length and quality of the sequence is calculated (len and mean_phread_quality)"
            """
            self.record=SeqIO.read(filename,filetype)
            self.filename=filename
            self.successfullyParsed=True
                  
            self.mean_phred_quality=int(sum(self.record.letter_annotations['phred_quality'])/max(len(self.record.letter_annotations['phred_quality']),1))
            self.len=len(self.record.seq)
            self.seq=self.record.seq
                
            if(self.len<50 or self.mean_phred_quality<12):
                self.comment="Quality: " + str(self.mean_phred_quality) + ", read length: "+str(self.len)+"."
                self.successfullyParsed=False
                self.chain_type="n/d"
            else:
                self.comment="ok"
                self.IgBlastMe(igblast_output)
               
                try:
                    if(self.BlastedOutputDict['top_v']=="N/A" and self.BlastedOutputDict['top_d']=="N/A" and self.BlastedOutputDict['top_j'][0]=="N/A"):
                        self.comment="No hits. "
                        self.successfullyParsed=False
                    else:
                        try:
                            self.aligned_start=int(self.BlastedOutputDict['v_hits'][0]['rank_1']['q_start'])
                        except:
                            print("Internal error while parsing "+ self.filename + ": Could not create complete aligned sequence.")
                            self.comment=" Internal error while parsing "+ self.filename + ": Could not create complete aligned sequence. "
                            self.successfullyParsed=False
                            self.chain_type="n/d"

                        try:
                            self.aligned_end=int(self.BlastedOutputDict['j_hits'][0]['rank_1']['q_end'])
                        except:
                            try:
                                try:
                                    self.aligned_end=int(self.BlastedOutputDict['d_hits'][0]['rank_1']['q_end'])
                                except:
                                    self.aligned_end=int(self.BlastedOutputDict['v_hits'][0]['rank_1']['q_end'])
                            except:
                                print("internal error")
    

                        if(self.BlastedOutputDict['strand']=="-"):
                            self.oriented_seq=self.seq.reverse_complement()
                        else:
                            self.oriented_seq=self.seq

                        self.aligned_seq=self.oriented_seq[self.aligned_start-1:self.aligned_end]
                    
                        if(len(self.aligned_seq)<80):
                            self.successfullyParsed=False
                            self.comment="IgBlast aligned less than 80 nt. Probably sequencing quality is bad."
                    
                        try:
                            self.gene_seq=""
                            temp=sorted(self.BlastedOutputDict['gene_alignments'], key=lambda x: self.BlastedOutputDict['gene_alignments'][x]['start'])
                            for reg in temp:
                                self.gene_seq+=self.BlastedOutputDict['gene_alignments'][reg]['seq']
                        except:
                            print("Internal error while parsing "+ filename + ": Could not create complete aligned sequence.")
                except AttributeError as my_err:
                    self.successfullyParsed=False
                    print("Internal error while parsing "+ filename + ". ")
                    self.comment=" Internal error while parsing "+ filename + ". "
                    self.chain_type="n/d"
                    print(my_err)
    
        def IgBlastMe(self,igblast_output=None):

            if(not os.path.isdir("./internal_data")):
                try:
                    shutil.copytree(cfg.igblast_internal_data, './internal_data')
                except OSError as my_err:
                    print("An error occured while IgBlasting: " + str(my_err))

            tmpFastaToBeBlasted=tempfile.NamedTemporaryFile(mode='w+t', delete=False)

            '''18-01-09: bugfix for running IgParser on Windows: tmpfiles must be created with delete=False to be able to be accesse twice'''
            tmpIgBlastOutput=tempfile.NamedTemporaryFile(mode='w+t', delete=False)
            tmpBlastOutput=tempfile.NamedTemporaryFile(mode='w+t', delete=False)

            self.writeToFasta(tmpFastaToBeBlasted.name)

            blast_cline=[cfg.blast_path]
            igblast_cline=[cfg.igblast_path]

            igblast_args_dict={
                "-germline_db_V ": '"'+cfg.germlinedb_V+'"',
                "-germline_db_J ": cfg.germlinedb_J,
                "-germline_db_D ": cfg.germlinedb_D,
                "-auxiliary_data ": cfg.igblast_auxiliary_data,
                "-domain_system": "imgt",
                #only show one alignement for one germline sequence
                "-num_alignments_V": "1",
                "-num_alignments_J": "1",
                "-num_alignments_D": "1",
                "-outfmt": '"7 std qseq sseq btop"',
                "-query ": tmpFastaToBeBlasted.name,
                "-out ": tmpIgBlastOutput.name
            }

            blast_args_dict={
                "-db ": cfg.constantdb,
                "-task ": "blastn",
                "-dust ": "no",
                "-outfmt": '"7 qseqid sseqid evalue bitscore"',
                "-max_target_seqs": "1",
                "-query ": tmpFastaToBeBlasted.name,
                "-out ": tmpBlastOutput.name
            }
            
            for command in igblast_args_dict:
                igblast_cline.append(str(command))
                igblast_cline.append(str(igblast_args_dict[command]))
            
            for command in blast_args_dict:
                blast_cline.append(str(command))
                blast_cline.append(str(blast_args_dict[command]))
            
            try: 
                #for some reason, cline has to be "joined", else it won't call igblastn correctly. 
                #I do not understand this, nevertheless it works now
                #also, shell=True is needed which is strongly discouraged
                subprocess.run(" ".join(igblast_cline), shell=True,check=True) 
            except subprocess.CalledProcessError as my_error:
                self.comment="executing igblastn failed"
                print("executing igblast failed.")
                print(my_error.cmd)
                print("This happened while igblasting " + self.filename +". ")
                return
            #saving igblast_output
            if(igblast_output is not None): 
                shutil.copy(tmpIgBlastOutput.name, igblast_output)


            try: 
                #for some reason, cline has to be "joined", else it won't call blastn correctly. 
                #I do not understand this, nevertheless it works now
                #also, shell=True is needed which is strongly discouraged
                subprocess.run(" ".join(blast_cline),stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL, shell=True,check=True) 
                #pass
            except subprocess.CalledProcessError as my_error:
                self.comment="executing blastn failed"
                print("executing blast failed.")
                print(my_error.cmd)
                print("This happened while igblasting " + self.filename +". ")
                return

            try:
                self.BlastedOutputDict=LoadBlastedOutput(open(tmpIgBlastOutput.name)).return_dict()
            except Exception as e:
                self.comment="parsing igblastn output failed."
                print("parsing igblastn output failed: " + repr(e))
                print("This happened while igblasting " + self.filename +". ")
                return

            for line in open(tmpBlastOutput.name):
                if "Ig" in line:
                    #bit values for bad matches are usually around 21
                    self.IgSC_new_bit=float(line.strip().split()[3])
                    self.IgSC_new_evalue=float(line.strip().split()[2])
                    
                    if(self.IgSC_new_bit>50 and self.IgSC_new_evalue<0.0001):
                        self.IgSC_new = line.strip().split()[1]
                        break
 

            tmpIgBlastOutput.close()
            tmpBlastOutput.close()
            tmpFastaToBeBlasted.close()

            self.chain_type=self.BlastedOutputDict['chain_type'].strip("V")
                    
            if(self.chain_type=="H"):
                try:
                    self.IgSubClass=self.IgSC_new
                    if(self.IgSC_new_bit<80):
                        self.comment="Ig SC determination confidence lower than usual (ebit: " + str(self.IgSC_new_bit) + "). Inspect manually or sequence again."
                except:
                    #IgSubClass could not be identified by blasting against the constant parts. We can apply some heuristics to guess the SubClass
                    self.IgSubClass=""

                    #IgG1/IgG2 heuristics
                    if(0<self.seq.find("TTGGTGGAGGC")<200):
                        if(29<self.seq.find("TTGGTGGAGGC")-self.seq.find("GGAGGGT")<33):
                            self.IgSubClass="IgG1"
                            self.comment=("Ig SC identification is based on heuristics. Sequencing quality was bad in the beginning. If you depend on this IgSC to be correct, please consider resequencing. If this is not an IgG1, it most likely is an IgG2.") 
                        elif(29<self.seq.find("TTGGTGGAGGC")-self.seq.find("GCAGGGC")<33):
                            self.IgSubClass="IgG2"
                            self.comment=("Ig SC identification is based on heuristics. Sequencing quality was bad in the beginning. If you depend on this IgSC to be correct, please consider resequencing. If this is not an IgG2, it most likely is an IgG1.") 
                        else:
                            self.comment="To properly identify this IgSubClass, it is strongly recommended to resequence the respective Ig"
                            self.IgSubClass="IgG1/IgG2. To properly identify this IgSubClass, it is strongly recommended to resequence the respective Ig"

                    #IgG3/IgG4 heuristics
                    if(self.seq.find("TTGGTGGAAG")>-1 and self.seq.find("TTGGTGGAAG")<200):
                        self.IgSubClass="IgG3/IgG4. To properly identify this IgSubClass, it is strongly recommended to resequence the respective Ig"

                    #IgA1/IgA2 heuristics
                    tgctg=self.seq.find("TGCTG")
                    if(0<tgctg < 100):
                        if(self.seq.find("TGCTGCAGAG")==tgctg):
                            self.IgSubClass="IgA1"
                        elif(self.seq.find("TGCTGTCGAG")==tgctg):
                            self.IgSubClass="IgA2"
                    elif(self.IgSubClass is ""):#this is kind of a last resort step: Ig SubClass is not set and we don't find the characteristic IgA motive.
                    #maybe we still find unique sequences for IgA1 or IgA2, usually a few bp before the tgctg motive?
                        if(self.seq.find("GGCGATGACCACGTTCCCATCTGGCTG") is not -1 and self.seq.find("GGCGATGACCACGTTCCCATCTGGCTG")< 100):
                            self.IgSubClass="IgA1"
                            self.comment+="Ig SC identification is based on heuristics. Sequencing quality was bad in the beginning. If you depend on this IgSC to be correct, please consider resequencing."
                        if(self.seq.find("TGCGACGACCACGTTCCCATCTTGGGG") is not -1 and self.seq.find("GGCGACGACCACGTTCCCATCTTGGGG")< 100):
                            self.IgSubClass="IgA2"
                            self.comment+="Ig SC identification is based on heuristics. Sequencing quality was bad in the beginning. If you depend on this IgSC to be correct, please consider resequencing."
 

                    cdr3=Seq(self.BlastedOutputDict['cdr3_sequence'])
                    self.cdr3pos=self.record.seq.find(cdr3.reverse_complement())
                    
                    #The following block implements the IgM-Subclass identification (page 7 of Jakobs presentation with a minor modification: 
                    #we don't search for the end of the aligned sequence but the cdr3 region, since this is easier to implement). 
                    #This works as follows: After we excluded, that we could identify the sequence as either IgG or IgA (negative identifaction), 
                    #we look for the beginning of the cdr3 region: In IgM, this should begin within the first 80 bases of the sequence.
                    if(self.IgSubClass==""):
                        if(self.cdr3pos==-1):
                            self.IgSubClass="n/d"
                            self.comment="IgSC is probably IgM, but cdr3 couldn't be identified. Please check manually that this is correct."

                            #Wether we should look in the first 70, 80, or 100 nt for the CDR3 is not yet determined, see the email w/ subject AW: IgM-Identifikation (16.11.2017) for some arguments
                        elif(self.cdr3pos<100):#if cdr3pos is within the first 100 nucleotides, this is an IgM
                            self.IgSubClass="IgM"
                            if(self.cdr3pos>70):
                                print("Warning! While trying to determine the IgSC, of " + self.filename + " we encountered a problem:")
                                print("cdr3pos is : " + str(self.cdr3pos))
                                print("We determine the IgSC according to the following rule: 'If cdr3pos is within the first C nt, this is an IgM'. Whether the optimal value for the constant C is 70, 80, or 100, is not yet determined (see email exchange between Jakob and Momsen 16.11.2017")
                        else:
                            self.IgSubClass="n/d"
                            tgga=self.seq.find("TGGA")
                            tgctg=self.seq.find("TGCTG")
                            if(tgga!=-1 and tgga < 90 and self.cdr3pos > 110 and self.cdr3pos < 120):
                                self.comment="IgSC is n/d, probably it is a IgG."
                            if(tgctg!=-1 and tgctg < 110 and self.cdr3pos > 130 and self.cdr3pos < 140):
                                self.comment="IgSC is n/d, probably it is a IgA."

        def __str__(self):
            return self.record.seq

        def writeToFasta(self,filename):
            SeqIO.write(self.record,filename,"fasta") 

        def identify_gene_region(self, pos):
            """returns the gene region (possible return values: fr1,cdr1,fr2,cdr2,fr3,cdr3, d,j, v_d_region, d_j_junction, v_j_junction)
            parameter pos is position relative to the query sequence
            """
            region=self.identify_V_gene_subregion(pos)
            if(region=="n/d"):
                region=self.identify_gene(pos)
            return region
        
        def identify_gene(self, pos):
            """returns the gene region (possible return values: v,d,j, v_d_region, d_j_junction, v_j_junction)
            parameter pos is position relative to the query sequence
            """
            for region in self.BlastedOutputDict['gene_alignments']:
                lower=int(self.BlastedOutputDict['gene_alignments'][region]['start'])
                upper=int(self.BlastedOutputDict['gene_alignments'][region]['end'])
                if(lower<=pos<=upper):
                    return region
            return "n/d"

        def identify_V_gene_subregion(self,pos):
            """returns the V gene subregion (possible return values: fr1,cdr1,fr2,cdr2,fr3,cdr3,n/d)
            parameter pos is position relative to the query sequence
            2018_04_18: Be aware, that this function returns "cdr3" even if the position does not lie in the V gene portion of the cdr3
            """
            for region in self.BlastedOutputDict['alignment_summaries']:
                if(region=="total"):
                    continue
                lower=int(self.BlastedOutputDict['alignment_summaries'][region]['from'])
                upper=int(self.BlastedOutputDict['alignment_summaries'][region]['to'])
                if(lower<=pos<upper):
                    return region
            return "n/d"
        
        def original_nt(self,pos):
            """returns the nt at pos of the gene sequence
            parameter pos is position relative to the V gene start
            """
            offset=int(self.BlastedOutputDict['gene_alignments']['v']['start'])
            return self.gene_seq[pos-offset]
        
        def IsThisNTSilent(self,pos):
            """returns True, if the aligned nt is silent (ie either is present in the aligned gene, or the nucleotide change does not lead 
            to an AA change. 
            parameter pos is position relative to the beginning of the V gene
            """
            tmp=self.translatedAA(pos)
            if(tmp[0]==tmp[1]):
                return True
            else:
                False
            
        def translatedAA(self,pos):
            """returns the translation of the aligned sequence and the gene sequence
            parameter pos is position relative to the beginning of the V gene
            """
            if(pos>len(self.gene_seq)):
                return ["n/d","n/d"]
            aapos=int(pos/3)
            offset=self.BlastedOutputDict['v_hits'][0]['rank_1']['s_start']-1
            tmp=int(aapos*3-offset)
            try:
                aa_gene=Seq(self.gene_seq[tmp:tmp+3]).translate()
            except TranslationError:
                #this is usually because the nucleotide is not determined at pos.
                aa_gene="*"
            try:
                aa_aligned=self.aligned_seq[tmp:tmp+3].translate()
            except TranslationError:
                aa_aligned="n/d"

            return [aa_gene,aa_aligned] 




class exportDict(collections.MutableMapping):
    """
    this class is an 'export dictionary' created from a parsed sequence passed to the constructor.
    it can be used like a dictionary (i.e. given to a DictWriter) and has all the extracted 
    information from a parsed sequence in it 
    arguments:
        parsed_sequence is a IgParser object
    """

    def __init__(self, parsed_sequence):
        self.store = dict()

        self['Seq_ID']=parsed_sequence.filename.split(".")[0]
        self['QV']=parsed_sequence.mean_phred_quality
        self['RL']=str(parsed_sequence.len)
        self['Confirmation']="NO"

        if(parsed_sequence.successfullyParsed==False  or hasattr(parsed_sequence,'BlastedOutputDict') is not True ):
            print("The parsed sequence " + parsed_sequence.filename + " was not blasted. The reason is: "+ parsed_sequence.comment)
            return

        self.chain_type=parsed_sequence.chain_type

        if(self.chain_type=="H"):
            temp_array= ["v","d","j"]
        else:
            temp_array= ["v","j"]

        for temp in ["v","d","j"]:
            try:
                ###The following snippet is legacy, but left there for a reason: we used to have more than one alignement per germline gene to parse.
                ###maybe we wanna go back to that flexibility in the future. 
                #first we need to find out if there is more than one match:
                if isinstance(parsed_sequence.BlastedOutputDict['top_'+temp],str):
                    #export the field IGH... and strip the first four letters IGH+[V,D,J]
                    self['IG'+self.chain_type.upper()+temp.upper()]=  parsed_sequence.BlastedOutputDict['top_'+temp].replace("IG"+self.chain_type.upper() +temp.upper(),"")
                else:# we have more than one match, so we take the first
                    #export the field IGH... and strip the first four letters IGH+[V,D,J]
                    self['IG'+self.chain_type.upper()+temp.upper()]=  parsed_sequence.BlastedOutputDict['top_'+temp][0].replace("IG"+self.chain_type.upper()+temp.upper(),"")
            except AttributeError:
                print("Ooops! Something went terribly wrong. Blasting the parsed file didn't work, maybe it is not an Ig?")
                print("I'll print the blasted dictionary for your convenience:")
                print(parsed_sequence.BlastedOutputDict)
                print("but you know, the chain type is: " + parsed_sequence.BlastedOutputDict['chain_type'])
            
        if(self.chain_type=="H"):#it is a heavy chain
            self['IgSC']=parsed_sequence.IgSubClass
        
        #TODO: was hat es mit partial_cdr_aa auf sich? sollte das nicht was ähnliches sein? hoffentlich nicht der reverse strand?
        
        if(self.chain_type=="H"):
            self['CDR3 IGHV']=parsed_sequence.BlastedOutputDict['cdr3_translated_sequence']
            #this is a really ugly workaround
            self['CDR3']=parsed_sequence.BlastedOutputDict['cdr3_translated_sequence']
            if(self['CDR3 IGHV']=="N/A"):
                self['CDR3L']=""
            else:
                self['CDR3L']=len(self['CDR3 IGHV'])

        else:
            self['CDR3']=parsed_sequence.BlastedOutputDict['cdr3_translated_sequence']
            if(self['CDR3']=="N/A"):
                self['CDR3L']=""
            else:
                self['CDR3L']=len(self['CDR3'])

        try:
            self['SHM']=parsed_sequence.BlastedOutputDict['alignment_summaries']['total']['mismatches']
            gaps=parsed_sequence.BlastedOutputDict['alignment_summaries']['total']['gaps'] 
            if(gaps is not 0):
                self['SHM']=str(self['SHM']) + " (+" + str(gaps) + " gaps)"
        except:
            self['SHM']="N/A"

        self['SHM IG'+self.chain_type.upper()+'V']=self['SHM']

        #for legacy reasons this is also supported
        self['SHM IGV'+self.chain_type.upper()]=self['SHM']

        if(parsed_sequence.BlastedOutputDict['productive']=="Yes"):
            self['Function']="Y"
        else:
            self['Function']="N"
 
        try:
            if(self.chain_type=="H"):
                self["5' Primer"]=primer_IGHV[self['IGHV']]
            elif(self.chain_type=="K"):
                self["5' Primer"]=primer_IGKV[self['IGKV']]
            elif(self.chain_type=="L"):
                self["5' Primer"]=primer_IGLV[self['IGLV']]
        except KeyError:
            if(self['Function']=="N"):
                self["5' Primer"]=""
            else:
                self["5' Primer"]="Primer not yet included"

        try:
            if(self.chain_type=="H"):
                self["3' Primer"]=primer_IGHJ[self['IGHJ']] 
            elif(self.chain_type=="K"):
                self["3' Primer"]=primer_IGKJ[self['IGKJ']]
            elif(self.chain_type=="L"):
                
                self["3' Primer"]= primer_IGLJ[self['IGLJ']] #"L3 2/3" 
        except KeyError:
            if(self['Function']=="N"):
                self["3' Primer"]=""
            else:
                self["3' Primer"]="Primer not yet included"

        for name, motif in restriction_motif.items():
            motif_pos= parsed_sequence.seq.find(motif) 
            if(motif_pos!=-1):
                self[name]="Y("+str(motif_pos)+")"
                
                #AgeI is included in some primers (at the end of the read: if the AgeI motif is found after pos 420 (Heavy) or 330 (lambda), and it is followed by TGC or AGC, we're okay with it.
                if(name=="AgeI"):
                    if(self.chain_type=="H"):
                        if (parsed_sequence.IgSubClass=="IgM"):
                            if(385<motif_pos==parsed_sequence.seq.find(motif+"TGC")):
                                self[name]="Y(P,"+str(motif_pos)+")"
                        else:
                            if(430<motif_pos==parsed_sequence.seq.find(motif+"TGC")):
                                self[name]="Y(P,"+str(motif_pos)+")"
                    if(self.chain_type=="L"):
                        if(330<motif_pos==parsed_sequence.seq.find(motif+"AGC")):
                            self[name]="Y(P,"+str(motif_pos)+")" 
            else:
                self[name]="N"
                    
    def __getitem__(self, key):
        return self.store[self.__keytransform__(key)]

    def __setitem__(self, key, value):
        self.store[self.__keytransform__(key)] = value

    def __delitem__(self, key):
        del self.store[self.__keytransform__(key)]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def __keytransform__(self, key):
        return key


def updateExcelRow(wb, row, columndict, parsed_sequence):
    ws=wb.active
    ed=exportDict(parsed_sequence)
    
    ed['Comment']=parsed_sequence.comment
    ed['Confirmation']="to be confirmed" 

    for key in columndict.keys():
        try:
            ws[columndict[key]+str(row)]= ed[key]
        except KeyError as my_err:
            sys.exit("KeyError was issued. exportDict does not know what to do with key: " + str(my_err) + ". This happened while exporting " + parsed_sequence.filename + ". Aborting...")



class AlignPCRObject():
    """
    The purpose of this class is to compare the two blasts of pcr1 and pcr2, one of them usually corresponding to a plasmid sequencing. pcr1 is usually done with a primer mix, so it might contain
    somatic hypermutations which are not present in pcr2 (sequenced with just one primer).
    """
    def __init__(self, pcr1, pcr2):
        self.output=""
        self.pcr1=pcr1
        self.pcr2=pcr2
        self.shmanalysis=""

        if(pcr1.mean_phred_quality<20):
            self.output="BQ  - " + pcr1.filename
            self.shmananalysis="n/a"
            return
        if(pcr2.mean_phred_quality<20):
            self.output="BQ  - " + pcr2.filename
            self.shmananalysis="n/a"
            return
        #if(pcr2.chain_type!="n/d"):
        for ct in "H","K","L":
            if(pcr2.seq.find(empty_vec_seq[ct])!=-1):
                self.output+="empty vector - " + pcr2.filename
                self.shmananalysis="n/a"
                return
        if(pcr2.chain_type=="N/A"):
            self.output=pcr2.filename + " is not an immunglobulin chain."
            self.shmananalysis="n/a"
            return            
        elif(pcr1.chain_type!=pcr2.chain_type):
            self.output="chain types differ! " + pcr1.filename + ": " + pcr1.chain_type + "C, " + pcr2.filename + ": " + pcr2.chain_type + "C."
            self.shmananalysis="n/a"
            return

        D1=pcr1.BlastedOutputDict
        D2=pcr2.BlastedOutputDict


        if(D1['top_j']!=D2['top_j'] and D1['top_v']!=D2['top_v'] ):
            self.output+="Diff " + pcr1.chain_type + "C (Both J and V genes do not match. These are most likely completely different chains)"
            self.shmananalysis="n/a"
            return
        elif(D1['top_j']!=D2['top_j']):
            self.output+="J genes do not match. "
        elif(D1['top_v']!=D2['top_v']):
            self.output+="V genes do not match. "
            
        if(D1['productive'].lower()!='yes'):
            self.output+=pcr1.filename + " is not productive. "
        if(D2['productive'].lower()!='yes'):
            self.output+=pcr2.filename + " is not productive. "

        #offset_pcr1 is the offset of the alignment of the 1st pcr with primer mix: since there is overlap by the 
        #primers, the igblast alignment sequence <-> gene usually has has a small offset of ~10 nt
        #offset_pcr2 is the same, if nothing goes wrong it is 0 (since the second pcr uses a different primer and sequences in
        #the reverse direction
        self.offset_pcr1=int(D1['v_hits'][0]['rank_1']['s_start'])-1
        self.offset_pcr2=int(D2['v_hits'][0]['rank_1']['s_start'])-1
        self.pcr2_alignment_start=int(D2['gene_alignments']['v']['start'])
        
        #offset is the offset between pcr1 and pcr2
        offset=self.offset_pcr1-self.offset_pcr2
        #offset is usually >=0, since pcr2 has better quality than pcr1 and usually extends the aligned sequence

        nt_tolerance_forward=nt_tolerance_fwd[pcr2.chain_type]
        nt_tolerance_reverse=nt_tolerance_rev[pcr2.chain_type]

        shmj_primer_canceled=0
        shmj_primer_added=0
        shmfr1_primer_canceled=0
        nonsilent_mutations_added=collections.defaultdict(int)
        silent_mutations_added=collections.defaultdict(int)
        nonsilent_mutations_canceled=collections.defaultdict(int)
        nonsilent_mutations_exchanged=collections.defaultdict(int)
        silent_mutations_canceled=collections.defaultdict(int)
              

        for index, letter in enumerate(pcr1.aligned_seq):
            if(offset+index<0):
                #offset is usually >=0, since pcr2 has better quality than pcr1 and usually extends the aligned sequence
                #if this is not the case, we skip the first nucleotides
                continue
            try:
                if(pcr2.aligned_seq[offset+index]==letter):
                    continue
            except IndexError:
                #print("offset:" + str(offset) + " index: " + str(index) + "letter:" + str(letter))
                #print(pcr2.aligned_seq[offset+index-5:offset+index-1])
                try: 
                    last_compared=str(pcr2.aligned_seq[offset+index-7:offset+index-1])
                except IndexError:
                    last_compared= str(pcr2.aligned_seq[offset+index-1])
                self.output+=("There was an index error while comparing sequences. This could be either due to an extra nucleotide in PCR2 (then usually you'll find lots of differences between the pcr2 and the plasmid sequence), or if the plasmid sequence terminates prematurely, e.g. if the sequencing quality is not sufficient. Until to the termination of the comparison, there were: ")
                break
            
            #TODO:bessere Namen wählen 
            pos_pcr2=offset+index+self.pcr2_alignment_start 
            if(pcr2.aligned_seq[offset + index]==pcr2.original_nt(pos_pcr2)):
                if(offset+index<nt_tolerance_forward):#this specific discrepancy (2nd pcr has same nt like gene, 
                    #1st pcr has mutation, nt is in the beginning) is due to primer overlap
                    shmfr1_primer_canceled=shmfr1_primer_canceled+1
                elif(offset+index> len(pcr2.gene_seq) - nt_tolerance_reverse):#this specific discrepancy (2nd pcr has same nt like gene, 
                    #1st pcr has mutation, nt is at the end of J gene) is due to primer overlap
                    shmj_primer_canceled=shmj_primer_canceled+1
                else:
                    region=pcr2.identify_gene_region(pos_pcr2)
                    #if(region=="v"):
                    #    region=pcr2.identify_V_gene_subregion(pos_pcr2)
                    
                    translation1=self.pcr1.translatedAA(index+offset)
                    translation2=self.pcr2.translatedAA(offset+index)
                        
                    if(translation2==translation1):
                        silent_mutations_canceled[region]+=1
                    else:
                        nonsilent_mutations_canceled[region]+=1 
            elif(letter==pcr2.original_nt(pos_pcr2)):
                if(offset+index> len(pcr2.gene_seq) - nt_tolerance_reverse):#this specific discrepancy (2nd pcr has same nt like gene, 
                    #1st pcr has mutation, nt is at the end of J gene) is due to primer overlap
                    shmj_primer_added=shmj_primer_added+1
                #print("pcr1 is correct at pos " + str(offset+index) + " sequence in pcr1 is: " + letter) 
                else:
                    region=pcr2.identify_gene_region(pos_pcr2)
                    #if(region=="v"):
                    #    region=pcr2.identify_V_gene_subregion(pos_pcr2)
                    
                    translation1=self.pcr1.translatedAA(index+offset)
                    translation2=self.pcr2.translatedAA(offset+index)

                    if(translation2==translation1):
                        silent_mutations_added[region]+=1
                    else:
                        if(self.pcr2.IsThisNTSilent(offset+index)):
                        #pcr2 introduced a new mutation relative to pcr1, which is nonsilent 
                        #it is "silent" compared to the gene sequence though, we count it as _silent_
                            silent_mutations_added[region]+=1
                        else:
                            nonsilent_mutations_added[region]+=1 

            else:#both pcr1 and pcr2 differ from the 'predicted' gene sequence
                region=pcr2.identify_gene_region(pos_pcr2)
                #if(region=="v"):
                #    region=pcr2.identify_V_gene_subregion(pos_pcr2)
        
                if(self.pcr2.IsThisNTSilent(offset+index)):
                    #pcr2 introduced a new mutation relative to pcr1, which is nonsilent 
                    #it is "silent" compared to the gene sequence though, we count it as _silent_
                    silent_mutations_added[region]+=1
                else:
                    #if(self.pcr1.IsThisNTSilent(index)==self.pcr2.IsThisNTSilent(offset+index)):
                    #
                    #print(index)
                    #print("letter: "+letter)
                    #print("plasmid: [" + pcr2.aligned_seq[offset + index] + "]"+ pcr2.aligned_seq[offset + index+1]+ pcr2.aligned_seq[offset + index+2] )
                    #print("gene: [" + pcr2.original_nt(pos_pcr2)+"]"+pcr2.original_nt(pos_pcr2+1)+pcr2.original_nt(pos_pcr2+2)+pcr2.original_nt(pos_pcr2+3)+pcr2.original_nt(pos_pcr2+4)+pcr2.original_nt(pos_pcr2+5))
                    nonsilent_mutations_exchanged[region]+=1
                    #nonsilent_mutations_added[region]+=1


        if(shmfr1_primer_canceled>0):
            self.output+=str(shmfr1_primer_canceled) + " SHM- FR1(P)"
        for reg in silent_mutations_added:
            self.output+=" " + str(silent_mutations_added[reg]) + " sSHM+ " + reg + " "
        for reg in nonsilent_mutations_added:
            self.output+=" " + str(nonsilent_mutations_added[reg]) + " nsSHM+ " + reg + " "
        for reg in silent_mutations_canceled:
            self.output+= " "+ str(silent_mutations_canceled[reg]) + " sSHM- " + reg + " "
        for reg in nonsilent_mutations_exchanged:
            self.output+=" " + str(nonsilent_mutations_exchanged[reg]) + " nsSHMchg " + reg + " "
        for reg in nonsilent_mutations_canceled:
            self.output+=" " + str(nonsilent_mutations_canceled[reg]) + " nsSHM- " + reg + " "
        if(shmj_primer_canceled>0):
            self.output+=" " + str(shmj_primer_canceled) + " SHM J(P)-"
        if(shmj_primer_added>0):
            self.output+=" "+ str(shmj_primer_added) + " SHM J(P)+"
        if(self.output==""):
            self.output="0" 

        #Finally, we check if the "cloning primer sequence" is contained in the sequence. This is necessary, since the Gibson HiFi Assembly might produce errors
        # which can not be caught by the above method, since they lie constant part. This is beta and only enabled on the H3 primer set
        ed1=exportDict(pcr1)
        ed2=exportDict(pcr2)


        if(pcr2.chain_type=="H"):
            if(pcr2.seq.find("ATGGGATGGTCATGTATCATCCTTTTTCTAGTAGCAACTGCAACCGGTGTACATTC")==-1):
                self.output+="CAVE: Likely Mutation in 5' primer. "
            if(pcr2.seq.    find("TCAGCGTCGACCAAGGGCCCATCGGTCTTCCCCCTGGCACCCTCC")==-1):
                self.output+="CAVE: Likely Mutation in 3' primer. "
        if(pcr2.chain_type=="K"):
            if(pcr2.seq.find("ATGGGATGGTCATGTATCATCCTTTTTCTAGTAGCAACTGCAACCGGTGTACATT")==-1 and pcr2.seq.find("ATGGGATGGTCATGTATCATCCTTTTTCTAGTAGCAACTGCAACCGGTGTACATG")==-1):
                self.output+="CAVE: Likely Mutation in 5' primer. "
            if(pcr2.seq.find("ATCAAACGTACGGTGGCTGCACCATCTGTCTTCATCTTCCCGCCA")==-1 and pcr2.seq.find("ATTAAACGTACGGTGGCTGCACCATCTGTCTTCATCTTCCCGCCA")==-1):
                self.output+="CAVE: Likely Mutation in 3' primer. "
        if(pcr2.chain_type=="L"):
            if(pcr2.seq.find("ATGGGATGGTCATGTATCATCCTTTTTCTAGTAGCAACTGCAACCGGTTC")==-1):
                self.output+="CAVE: Likely Mutation in 5' primer. "
            if(pcr2.seq.find("CACTCTGTTCCCGCCCTCGAGTGAGGAGCTTCAAGCCAACAAGGCCACACTG")==-1 and pcr2.seq.find("CACTCTGTTCCCACCCTCGAGTGAGGAGCTTCAAGCCAACAAGGCCACACTG")==-1):
                self.output+="CAVE: Likely Mutation in 3' primer. "
                
        self.shmanalysis="Total SHM according to 2nd pcr/plasmid: " + str(ed1['SHM']) + "/" + str(ed2['SHM'])
