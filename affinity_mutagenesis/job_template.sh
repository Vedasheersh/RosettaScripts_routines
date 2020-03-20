#PBS -l nodes=4:ppn=20

#PBS -l walltime=12:00:00

#PBS -l pmem=8gb

#PBS -l mem=8gb

#PBS -A cdm8_b_g_sc_default
#PBS -j oe

set -u

cd $PBS_O_WORKDIR

echo " "

echo " "

echo "JOB Started on $(hostname -s) at $(date)"

module load python/3.6.3-anaconda5.0.1 

python write_resfile.py {0} --antibody {1} --epitopes {2} --nearby_atom_cutoff {3} --output {4} --design_side {5} --native False --repack True

module unload python/3.6.3-anaconda5.0.1 

/opt/aci/sw/mpich/3.2_gcc-5.3.1/bin/mpirun -ppn 20 -np 80 /gpfs/group/cdm8/default/rosetta/new_rosetta_build/rosetta_bin_linux_2018.33.60351_bundle/main/source/bin/rosetta_scripts.mpi.linuxgccrelease @{6}

