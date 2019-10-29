#!/bin/python
import sys
try:
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
except ImportError:
	print( "Need Biopython to use the IgBlast output parser class")


"""
Originally taken from pyigblast
small modifications by Momsen Reincke in September 2017:
    -change dictionary of JSON (strip .) to be able to export JSON to mongodb)
bigger modifications in Jan 2018:
    - refactor the parsing of alignment summaries (parse all at once)
"""

class LoadBlastedOutput():
    '''The helper class to parse an individual blast result'''

    def __init__(self, query):
        _rearrangment_breaker = False
        _junction_breaker = False
        _fields_breaker = False

        # initalize the fields, some can be empty on crappy reads
        # basic
        self.query = ""
        self.domain_classification = ""

        # title fields
        self.rearrangment_summary_titles = ""
        self.alignment_summary_titles = ""
        self.junction_detail_titles = ""
        self.hit_fields = ""

        # summaries
        self.rearrangment_summary = ""
        self.junction_detail = ""
        self.alignment_summaries = {}

        # hits
        self.hits_v = []  # to be parsed in another function
        self.hits_d = []  # to be parsed in another function
        self.hits_j = []  # to be parsed in another function

        # will hold everything
        self.blast_dict = {}
        self.alignments = {}
        for line in query:
            if "Query" in line:
                self.query = line.split(":")[1].strip()
            if "Domain classification requested:" in line:
                self.domain_classification = line.split(":")[1].strip()
            if "rearrangement summary" in line:
                self.rearrangment_summary_titles = line.strip().split(
                    "(")[2].split(")")[0].split(",")
                _rearrangment_breaker = True
                continue
            if _rearrangment_breaker:
                self.rearrangment_summary = line.strip().split("\t")
                _rearrangment_breaker = False
            if "junction details" in line:
                self.junction_detail_titles = line.strip().split(
                    "(")[2].split(")")[0].split(",")
                _junction_breaker = True
                continue
            if _junction_breaker:
                self.junction_detail = line.strip().split("\t")
                _junction_breaker = False
            if "Alignment summary" in line:
                self.alignment_summary_titles = line.strip().split(
                    "(")[1].split(")")[0].split(",")
                self.alignment_summary_titles = self.alignment_summary_titles
            if line.startswith("FR1-IMGT"):
                self.alignment_summaries['fr1']= line.strip().split()[1:] 
            if line.startswith("CDR1-IMGT"):
                self.alignment_summaries['cdr1']= line.strip().split()[1:] 
            if line.startswith("FR2-IMGT"):
                self.alignment_summaries['fr2']= line.strip().split()[1:] 
            if line.startswith("CDR2-IMGT"):
                self.alignment_summaries['cdr2']= line.strip().split()[1:] 
            if line.startswith("FR3-IMGT"):
                self.alignment_summaries['fr3']= line.strip().split()[1:] 
            if line.startswith("CDR3-IMGT"):
                #here we have to start at split by 2:, since the line actually reads:
                #CDR3-IMGT (germline)	311	322	12	12	0	0	100
                #mind the (germline)
                self.alignment_summaries['cdr3']= line.strip().split()[2:] 
                
                """
                added by MR 2017-09-10
                the next continue is an elegant solution for the following problem:
                both the CDR3 line and CDR3-IMGT start with CDR3. To read out the CDR3 line, we want to not confuse
                it with the CDR3-IMGT line. We check CDR3-IMGT first, when we check CDR3 then later in the for loop,
                we can be sure it really is the CDR3 line and not the CDR3-IMGT line
                """
                continue
            if line.startswith("Total"):
                """
                added by MR 2017-09-11
                Unfortunately there are three lines starting with Total: the Total alignment summary, 
                but also "Total queries", "Total identifiable CDR3" and "Total unique clonotype"
                We need to take care of this. That's why we check total number of stripped items. We are looking
                for lines that like the following (and therefore strip down to 8 items):
                
                    Total	N/A	N/A	295	274	21	0	92.9
                added by MR 2019-05-01
                we parse the CDR3 info as quality control
                """
                if("CDR3" in line):
                    self.total_identifiable_cdr3=line.strip().split()[4]
                elif(len(line.strip().split())==8):    
                    self.alignment_summaries['total'] = line.strip().split()[1:]
              
            if line.startswith("CDR3"):
                self.cdr3_sequence = line.strip().split("\t")[1]
                self.cdr3_translated_sequence = line.strip().split("\t")[2]
                self.cdr3_start = line.strip().split("\t")[3]
                self.cdr3_end = line.strip().split("\t")[4]

            if "# Fields:" in line:
                self.hit_fields = line.strip().split(":")[1].split(",")
                _fields_breaker = True
            if _fields_breaker:
                if line.startswith("V"):
                    self.hits_v.append(line)
                elif line.startswith("D"):
                    self.hits_d.append(line)
                elif line.startswith("J"):
                    self.hits_j.append(line)

        self.process()

    def process(self):
        self.blast_dict[
            self.query] = {"domain_classification": self.domain_classification,
                           "rearrangement": self.parse_rearranment(),
                           "junction": self.parse_junction(),
                           "alignment_summaries": self.parse_alignment_summaries(),
                           "v_hits": self.parse_v_hits(),
                           "d_hits": self.parse_d_hits(),
                           "j_hits": self.parse_j_hits()
                           }
        self.translate_junction()

    def translate_junction(self):
    	self.cdr3_partial = ""
    	if self.junction_together:
    		coding_region = Seq(self.junction_together,IUPAC.ambiguous_dna)
    		self.cdr3_partial = str(coding_region.translate())

    def parse_rearranment(self):
        _return_dict = {}
        for title, value in zip(self.rearrangment_summary_titles, self.rearrangment_summary):
            if len(value.split(',')) > 1:
                #cast multiple entries for tuple, makes them easier for json
                _return_dict[title.strip()] = tuple(value.split(','))
            else:
                _return_dict[title.strip()] = value
        return _return_dict

    def parse_junction(self):
        _return_dict = {}
        self.junction_together = ""
        for title, value in zip(self.junction_detail_titles, self.junction_detail):
            if "(" in value:
                _return_dict["d-or-j_junction"] = value.split(
                    "(")[1].split(")")[0]
                self.junction_together += value.split(
                    "(")[1].split(")")[0]
            else:
                _return_dict[title.strip().lower().replace(" ", "_")] = value
               	if value != "N/A":
               		self.junction_together += value
        return _return_dict

 
    def parse_alignment_summaries(self):
        _return_dict = {}
        for region in self.alignment_summaries.keys():
            _return_dict[region]={}
            for title, value in zip(self.alignment_summary_titles, self.alignment_summaries[region]):
                title=title.strip()#added by MR 2017-09-11: no whitespaces at the beginning
                try:
                    _return_dict[region][title] = int(value)
                except ValueError:
                    try:
                        _return_dict[region][title] = float(value)
                    except ValueError:
                        _return_dict[region][title] = value

        #2018-04-18: we treat cdr3 differently: since we are actually interested in the whole cdr3 region, not the cdr3 region which comes from the top V gene hit, we modify this a bit:
        #try:
        #2019-05-01: further modified this region to evaluate if we detected a cdr3 region or not
        if(self.total_identifiable_cdr3 is not "0"):
            if(not "cdr3" in _return_dict or int(_return_dict["cdr3"]['from'])!=int(self.cdr3_start)):
            #cdr3_start comes from the "sub-region sequence details" field of igblast, whereas "cdr3[from] comes from the
            #alignment summary. Usually they coincide, if they don't this is most likely due to bad quality/problems with the read
                print("WARNING: CDR3 region cannot be uniquely identified.")
                _return_dict["cdr3"]={}
            _return_dict["cdr3"]['to']=self.cdr3_end
            _return_dict["cdr3"]['from']=self.cdr3_start
        
        return _return_dict

    def parse_v_hits(self):
        _return_dict = {}
        rank = 1
        for entry in self.hits_v:
            _entry_dict = {}
            #_entry_dict['rank'] = int(rank)
            for value, title in zip(entry.split()[1:], self.hit_fields):
                try:
                    #2017-09-10MR: stripping blanks and . to be able to export the json to a mongodb
                    _entry_dict[title.strip().replace(' ', '_').replace('.','')] = float(value)
                except:
                    _entry_dict[title.strip().replace(' ', '_').replace('.','')] = value
            _return_dict['rank_' + str(rank)] = _entry_dict
            rank += 1
        return _return_dict

    def parse_d_hits(self):
        _return_dict = {}
        rank = 1
        for entry in self.hits_d:
            _entry_dict = {}
            #_entry_dict["rank"] = int(rank)
            for value, title in zip(entry.split()[1:], self.hit_fields):
                try:
                    #2017-09-10MR: stripping blanks and . to be able to export the json to a mongodb
                    _entry_dict[title.strip().replace(' ', '_').replace('.','') ] = float(value)
                except:
                    _entry_dict[title.strip().replace(' ', '_').replace('.','')] = value
            _return_dict['rank_' + str(rank)] = _entry_dict
            rank += 1
        return _return_dict

    def parse_j_hits(self):
        _return_dict = {}
        rank = 1
        for entry in self.hits_j:
            _entry_dict = {}
            #_entry_dict['rank'] = int(rank)
            for value, title in zip(entry.split()[1:], self.hit_fields):
                try:
                    #2017-09-10MR: stripping blanks and . to be able to export the json to a mongodb
                    _entry_dict[title.strip().replace(' ', '_').replace('.','')] = float(value)
                except:
                    _entry_dict[title.strip().replace(' ', '_').replace('.','')] = value
            _return_dict['rank_' + str(rank)] = _entry_dict
            rank += 1
        return _return_dict

    def return_dict(self):
        '''Our Main Function that will return a dictionary'''
        # to be converted to a json document
        self._return_dict = {}
        
        #hits arrays if we have more than one hit we kept in the blast query
        self.v_hits_array = []
        self.j_hits_array = []
        self.d_hits_array = []
        _gene_alignments={}

        self.json_dictionary = {
            "_id": self.query,
            "format": self.blast_dict[self.query]['domain_classification']
        }

        # Most important should be considered individually
        try:
            self._return_dict["top_v"] = self.blast_dict[
                self.query]['rearrangement']['Top V gene match']
        except KeyError:
            self._return_dict["top_v"] = "N/A"
        try:
            self._return_dict["top_d"] = self.blast_dict[
                self.query]['rearrangement']['Top D gene match']
        except KeyError:
            self._return_dict["top_d"] = "N/A"
        try:
            self._return_dict["top_j"] = self.blast_dict[
                self.query]['rearrangement']['Top J gene match']
        except KeyError:
            self._return_dict["top_j"] = "N/A",
        try:
            self._return_dict["strand"] = self.blast_dict[
                self.query]['rearrangement']['Strand']
        except KeyError:
            self._return_dict["strand"] = "N/A"
        try:
            self._return_dict["chain_type"] = self.blast_dict[
                self.query]['rearrangement']['Chain type']
        except KeyError:
            self._return_dict["chain_type"] = "N/A"
        try:
            self._return_dict["stop_codon"] = self.blast_dict[
                self.query]['rearrangement']['stop codon']
        except KeyError:
            self._return_dict["stop_codon"] = "N/A"
        try:
            self._return_dict["productive"] = self.blast_dict[
                self.query]['rearrangement']['Productive']
        except KeyError:
            self._return_dict["productive"] = "N/A"
        try:
            self._return_dict["in_frame"] = self.blast_dict[
                self.query]['rearrangement']['V-J frame']
        except KeyError:
            self._return_dict["in_frame"] = "N/A"

        # add junctions. won't modify key names this time
        try:
            for junction_entry in self.blast_dict[self.query]['junction']:
                self._return_dict[junction_entry] = self.blast_dict[
                    self.query]['junction'][junction_entry]
        except KeyError:
            for junction_title in self.junction_detail_titles:
                self._return_dict[junction_title] = "N/A"

        #TODO: fill values with N/A if the alignmenets are empty (can this happen, actually?)
        #if(len(self.alignment_summaries))=0)
        #else:
        #    self.alignment_summaries = {
        #        "fr3_align": "N/A", "cdr3_align": "N/A", "total_align": "N/A"}
        self._return_dict['alignment_summaries'] = self.blast_dict[self.query]['alignment_summaries']

        # vhits
        try:
            for rank in sorted(self.blast_dict[self.query]['v_hits']):
                self.v_hits_array.append(
                    {rank: self.blast_dict[self.query]['v_hits'][rank]})
        except ValueError:
            self.v_hits_array = "N/A"

        # dhits
        try:
            for rank in sorted(self.blast_dict[self.query]['d_hits']):
                self.d_hits_array.append(
                    {rank: self.blast_dict[self.query]['d_hits'][rank]})
        except ValueError:
            self.d_hits_array = "N/A"

        # jhits
        try:
            for rank in sorted(self.blast_dict[self.query]['j_hits']):
                self.j_hits_array.append(
                    {rank: self.blast_dict[self.query]['j_hits'][rank]})
        except KeyError:
            self.j_hits_array = "N/A"

        self._return_dict["v_hits"] = self.v_hits_array
        self._return_dict["d_hits"] = self.d_hits_array
        self._return_dict['j_hits'] = self.j_hits_array

                
        for reg in ['j_hits','d_hits','v_hits']:
            try:
                _gene_alignments[reg[0]]={}
                _gene_alignments[reg[0]]['start']=self.blast_dict[self.query][reg]['rank_1']['q_start']
                _gene_alignments[reg[0]]['end']=self.blast_dict[self.query][reg]['rank_1']['q_end']
                _gene_alignments[reg[0]]['seq']=self.blast_dict[self.query][reg]['rank_1']['subject_seq']
            except KeyError:#this means we dont have a reg[0], i.e. for example no D gene was aligned
                del(_gene_alignments[reg[0]])
                continue
        
        """we need this counter, since sometimes igblast aligns a V gene, a J gene, NO D gene, no VD junction, but a DJ junction!
        we work around this by keeping track of how many nts are aligned between V J using vdj_junction-counter
        """
        vdj_junction_counter=0
        try:
            if(self._return_dict['v-d_junction']!="N/A"):
                _gene_alignments['v_d_junction']={}
                _gene_alignments['v_d_junction']['seq']=self._return_dict['v-d_junction']
                _gene_alignments['v_d_junction']['start']=_gene_alignments['v']['end']+1
                _gene_alignments['v_d_junction']['end']=_gene_alignments['v']['end']+len(self._return_dict['v-d_junction'])
                vdj_junction_counter=len(_gene_alignments['v_d_junction']['seq'])

        except:
            pass
        
        try:
            vdj_junction_counter+=len(_gene_alignments['d']['seq'])
        except:
            pass
            
        try:
            if(self._return_dict['d-j_junction']!="N/A"):
                _gene_alignments['d_j_junction']={}
                _gene_alignments['d_j_junction']['seq']=self._return_dict['d-j_junction']
                _gene_alignments['d_j_junction']['start']=_gene_alignments['v']['end']+1+vdj_junction_counter
                _gene_alignments['d_j_junction']['end']=_gene_alignments['v']['end']+vdj_junction_counter+len(self._return_dict['d-j_junction'])
        except:
            pass

        try:
            if(self._return_dict['v-j_junction']!="N/A"):
                _gene_alignments['v_j_junction']={}
                _gene_alignments['v_j_junction']['start']=_gene_alignments['v']['end']+1
                _gene_alignments['v_j_junction']['end']=_gene_alignments['v']['end']+len(self._return_dict['v-j_junction'])
                _gene_alignments['v_j_junction']['seq']=self._return_dict['v-j_junction']
        except:
                pass

        self._return_dict['total_identifiable_cdr3']=self.total_identifiable_cdr3

        try: 
            self._return_dict['cdr3_translated_sequence']=self.cdr3_translated_sequence 
        except AttributeError:
            self._return_dict['cdr3_translated_sequence']="N/A"
        try: 
            self._return_dict['cdr3_sequence']=self.cdr3_sequence 
        except AttributeError:
            self._return_dict['cdr3_sequence']="N/A"

        if self._return_dict["productive"].lower()  == "yes":
        	self._return_dict["partial_cdr3_aa"] = self.cdr3_partial

        self._return_dict['gene_alignments']=_gene_alignments

        
        return self._return_dict
