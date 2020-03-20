import os

PDB = 'WT1.pdb'

EPITOPES = ['346A','347A','348A','351A','352A','354A','355A']

PARTNERS = 'HK_A'

ANTIBODY = 'HK'

ANTIGEN = 'A'

PROTOCOL = 'dock_design.xml'

RESFILE_SCRIPT = 'write_resfile.py'

PIVOTS_SCRIPT = 'write_pivots.py'

DESIGN = 'antibody'

SUBMIT = False

# For each PDB we need 10,000 structures
# Split into 15 individual sets
structs = [700]*15

BASE_DIR = os.path.abspath(os.getcwd())
for set_num in range(1,2):
	# for each set, make a folder and copy files
	os.system('mkdir set_{0}'.format(set_num))
	os.system('cp {0} set_{1}/'.format(PDB,set_num))
	os.system('cp {0} set_{1}/'.format(PROTOCOL,set_num))
	os.system('cp {0} set_{1}/'.format(RESFILE_SCRIPT,set_num))
	os.system('cp {0} set_{1}/'.format(PIVOTS_SCRIPT,set_num))
	
	os.chdir('set_{0}'.format(set_num))
	
	# write and read pivots
	os.system('module load python/3.6.3-anaconda5.0.1')
	print(os.system('/opt/aci/sw/python/3.6.3_anaconda-5.0.1/bin/python write_pivots.py {0} --antibody {1} --output {2}'.format(PDB,ANTIBODY,'set_{0}'.format(set_num))))
	os.system('module unload python/3.6.3-anaconda5.0.1')

	f = open('set_{0}.pivots'.format(set_num))
	pivots = f.read().strip()
	f.close()
	
	# write flags file
	f = open('set_{0}.flags'.format(set_num),'w')
	flags_template = open('../flag_template.txt').read()
	f.write(flags_template.format(PDB,PROTOCOL,PARTNERS,structs[set_num-1],'_set_{0}'.format(set_num),'set_{0}.resfile'.format(set_num),pivots))
	f.close()
	
	# write job file
	f = open('set_{0}.job'.format(set_num),'w')
	job_template = open('../job_template.sh').read()
	f.write(job_template.format(PDB,ANTIBODY,','.join(EPITOPES),5.0,'set_{0}'.format(set_num),DESIGN,'set_{0}.flags'.format(set_num)))
	f.close()
	
	if SUBMIT:
		os.system('qsub set_{0}.job'.format(set_num))
		
	os.chdir(BASE_DIR)
	
	
