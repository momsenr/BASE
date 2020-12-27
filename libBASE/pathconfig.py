#!/bin/python
import platform
import os

# copy internal data directory to current location
if(platform.system()=="Linux"):
    igblast_internal_data=os.path.join("igblast","internal_data")   # path to igblast internal data directory
    igblast_auxiliary_data=os.path.join('igblast','optional_file','human_gl.aux') # path to igblast auxiliary data
    constantdb=os.path.join('human_gl_C.fasta') # path to Ig constant part db
    blast_path="blastn" # path to blast binary
    igblast_path="igblastn" # path to igblastn binary
elif(platform.system()=="Windows"):
    #igblast_internal_data="empty"   # path to igblast internal data directory
    igblast_auxiliary_data=os.path.join('ncbi-igblast-1.8.0','optional_file','human_gl.aux') # path to igblast auxiliary data 
    constantdb=os.path.join('germline_c_db',"germline_c_db")  # path to Ig constant part db 
    blast_path=os.path.join('bin','blast-2.9.0+','bin','blastn.exe') # path to blast binary 
    igblast_path="bin\igblastn.exe" # path to igblast binary
    
#os-independent paths
germlinedb_V=os.path.join('germlinedb','db','V_without_orphons')
germlinedb_J=os.path.join('germlinedb','db','J_without_orphons')
germlinedb_D=os.path.join('germlinedb','db','D_without_orphons')
