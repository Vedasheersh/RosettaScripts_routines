from Bio.PDB import PDBParser
from Bio.PDB import PDBExceptions
import warnings
import sys
import math
from optparse import OptionParser

standard_aa_names={"ALA":"A", "CYS":"C", "ASP":"D", "GLU":"E", "PHE":"F", "GLY":"G", "HIS":"H", "ILE":"I", "LYS":"K", 
                   "LEU":"L", "MET":"M", "ASN":"N", "PRO":"P", "GLN":"Q", "ARG":"R", "SER":"S", "THR":"T", "VAL":"V",
                   "TRP":"W", "TYR":"Y", "TYS":"Y"}
                   
## Returns the minimal distance between any pair of atoms between two residues
def min_distance( pose, residue1, residue2 ):
    min_dist = -1
    seqpos1, chain1, resid1 = residue1
    seqpos2, chain2, resid2 = residue2
    for atom1 in pose[ chain1 ][ int(seqpos1) ]:
        if 'H' in atom1.get_name():
            continue
        for atom2 in pose[ chain2 ][ int(seqpos2) ]:
            if 'H' in atom2.get_name():
                continue
            if atom1 - atom2 < min_dist or min_dist == -1:
                min_dist = atom1 - atom2

    return min_dist
    
def get_interface_residues( pose, antibody_chains, epitope_residues, nearby_atom_cutoff ):
    antibody_residues = []
    antigen_residues = []
    for chain in antibody_chains:
        for res in pose[ chain ]:
            seqpos = res.get_id()[1]
            resid  = standard_aa_names[ res.get_resname() ]
            antibody_residues.append( ( str(seqpos), chain, resid ) )

    for each in epitope_residues:
        chain = each[-1]
        resn = int(each[:-1])
        res = pose[chain][resn]
        resid  = standard_aa_names[ res.get_resname() ]
        antigen_residues.append( ( str(resn), chain, resid ) )
    
    antigen_interface = []
    antibody_interface = []
    for res1 in antigen_residues:
        for res2 in antibody_residues:
            if min_distance( pose, res1, res2 ) < nearby_atom_cutoff:
                antibody_interface.append( res2 )
            antigen_interface.append( res1 )
    return [ set( antibody_interface ), set( antigen_interface ) ]
    
if __name__ == '__main__':
    usage = "%prog [options] <pdb_file>"
    parser=OptionParser(usage)
    parser.add_option("--antibody",dest="antibody",help="the chains that make up antibody (as a string, e.g. 'AB')", default="")
    parser.add_option("--epitopes",dest="epitopes",help="the residues that make up antigen eptiopes (as a string, e.g. '12A,13A,25B')", default="")
    parser.add_option("--nearby_atom_cutoff",dest="nearby_atom_cutoff",help="SC distance cutoff to define a residue as part of the interface. \
    If any SC atom from a residue on one side is within this cutoff of a residue on the other side it's considered to be in the interface. Default=5.0",default=5.0)
    parser.add_option("--output",dest="output",default="out", help="Output name for resfile")
    parser.add_option("--native",dest="native",help='Just repack the residues on the side flagged "design side"', default=False)
    parser.add_option("--design_side",dest="design_side",default='antibody',help="antibody or antigen needs to be designed. \
    Allowed values:- antibody or antigen. Defaults to antibody. If antigen is specified, only eptiope residues are designed. \
    If antibody is specified all residues of antibody within cutoff from antigen epitope residues are designed")
    parser.add_option("--repack",dest="repack",default='True', help='Repack side of the interface not being designed. By default it is true')
    (options,args)= parser.parse_args()
    
    if len(args) < 1:
        parser.error('specify a pdb file or use -h')
    elif len(args) > 1:
        print('Warning: only the first pdb is considered by this script. The rest will be ignored')

    warnings.simplefilter('ignore',PDBExceptions.PDBConstructionWarning)    

    print('Processing',args[0])
    
    input_pdb = args[0]
    epitopes = options.epitopes.split(',')
    design_side = options.design_side
    native = (options.native.strip())
    output = options.output.strip()
    repack = (options.repack.strip())
    #print(repack)

    antibody_chains=options.antibody
    nearby_atom_cutoff=options.nearby_atom_cutoff
    
    parser = PDBParser()
    structure = parser.get_structure( 'X', input_pdb )

    antibody_int, antigen_int = get_interface_residues(structure[0],antibody_chains,epitopes,5)

    side_dict = {}
    side_dict[ 'antibody' ] = sorted( antibody_int, key=lambda res: ( int(res[0]), res[1] ) )
    side_dict[ 'antigen' ] = sorted( antigen_int, key=lambda res: ( int(res[0]), res[1] ) )

    if design_side=='antibody':
        opposing_side = 'antigen' 
    else:
        opposing_side = 'antibody'
        
    if native == 'True': 
        print('Native is {}\nSo I am True'.format(native))
        packer_aas = 'NATAA'  
    else: 
        packer_aas = 'ALLAA'
    
    #print(side_dict[ design_side ],side_dict[ opposing_side ])
    #print(packer_aas)
    #print(native)

    with open(output+'.resfile', 'w')  as out:
        out.write('NATRO\nstart\n')

        for seqpos, chain, resid in side_dict[ design_side ]:
            out.write( ' '.join([ seqpos, chain, packer_aas ])+'\n' )

        if repack=='True':
            for seqpos, chain, resid in side_dict[ opposing_side ]:
                out.write( ' '.join([ seqpos, chain, 'NATAA' ])+'\n' )
