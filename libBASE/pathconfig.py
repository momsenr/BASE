#!/bin/bash
import platform
import os

# copy internal data directory to current location
if(platform.system()=="Linux"):
    igblast_internal_data="/home/momsen/AGPruess/igblast/internal_data"
    igblast_auxiliary_data=os.path.join('igblast','optional_file','human_gl.aux')
    constantdb=os.path.join('human_gl_C.fasta')
    blast_path="blastn"
    igblast_path="igblastn"
elif(platform.system()=="Windows"):
    igblast_internal_data"\\be-svr-001-vm02-cos.charite.de\ag-pruess\Benutzer Ag Prüß\Momsen\software\internal_data"
    igblast_auxiliary_data=os.path.join('ncbi-igblast-1.8.0','optional_file','human_gl.aux')
    constantdb=os.path.join('germline_c_db',"germline_c_db")
    blast_path=os.path.join('bin','blast-2.9.0+','bin','blastn.exe')
    igblast_path="bin\igblastn.exe"
    
'''os-independent paths'''
germlinedb_V=os.path.join('germlinedb','db','V_without_orphons')
germlinedb_J=os.path.join('germlinedb','db','J_without_orphons')
germlinedb_D=os.path.join('germlinedb','db','D_without_orphons')
