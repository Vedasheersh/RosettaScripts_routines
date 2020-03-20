from Bio.PDB import PDBParser
from Bio.PDB import PDBExceptions
# from Bio.PDB import Polypeptide
import warnings
import sys
import math
from optparse import OptionParser

## 6/1/2015 AMS - This script was adapted from a similar one written by David Nannemann. The purpose is to define residues at the interface of a 
## protein-protein complex. It eliminates the use of a vector cutoff and uses a simpler calculation of cross-interface distances, either by 
## closest side chain atom, Ca distance, or Cb distance.

standard_aa_names={"ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E", "PHE":"F", "GLY":"G", "HIS":"H", "ILE":"I", "LYS":"K", 
                   "LEU":"L", "MET":"M", "ASN":"N", "PRO":"P", "GLN":"Q", "ARG":"R", "SER":"S", "THR":"T", "VAL":"V",
                   "TRP":"W", "TYR":"Y", "TYS":"Y"}
    
def get_antibody_residues( pose, antibody_chains ):
    antibody_residues = []
    for chain in antibody_chains:
		temp = []
        for res in pose[ chain ]:
            seqpos = res.get_id()[1]
            resid  = standard_aa_names[ res.get_resname() ]
            temp.append( ( str(seqpos), chain, resid ) )
        # leave out first and last since they cause problems in backrub
        antibody_residues.extend(temp[1:-1])
    return antibody_residues
    
if __name__ == '__main__':
    usage = "%prog [options] <pdb_file>"
    parser=OptionParser(usage)
    parser.add_option("--antibody",dest="antibody",help="the chains that make up antibody (as a string, e.g. 'AB')", default="")
    parser.add_option("--output",dest="output",default="out", help="Output name for pivots file")
    (options,args)= parser.parse_args()
    
    if len(args) < 1:
        parser.error('specify a pdb file or use -h')
    elif len(args) > 1:
        print('Warning: only the first pdb is considered by this script. The rest will be ignored')

    warnings.simplefilter('ignore',PDBExceptions.PDBConstructionWarning)    

    print('Processing',args[0])
    
    input_pdb = args[0]
    output = options.output.strip()
    antibody_chains=options.antibody
    
    parser = PDBParser()
    structure = parser.get_structure( 'X', input_pdb )

    antibody_residues = get_antibody_residues(structure[0],antibody_chains)
    
    with open(output+'.pivots', 'w')  as out:
        for seqpos, chain, resid in antibody_residues:
            out.write( ''.join([ seqpos.strip(), chain.strip() ])+',' )
