#!/usr/bin/env python3
from pyrosetta import *
from pyrosetta.rosetta import *
import os
import argparse


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', '-p', dest='pdb')
    parser.add_argument('--chain_num', '-c', dest='chain', default='2')
    parser.add_argument('--pep_list', '-l', dest='peplist')
    parser.add_argument('--refinement', '-r', dest='refine', action='store_true', default=False)
    parser.add_argument('--bb_min', '-b', dest='bbmin', action='store_true', default=False)
    return parser


def create_init_resfile(pdb_path, chain):
    """Create resfile with appropriate numbering with option PIKAA (sequence of the template)"""
    pose = pose_from_pdb(pdb_path)
    pdb_inf = pose.pdb_info()

    pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]

    first_res_pose = pose.chain_begin(chain)
    last_res_pose = pose.chain_end(chain)

    chain_name = pdb_inf.chain(first_res_pose)
    chain_length = last_res_pose - first_res_pose + 1

    with open('resfile_%s'%pdb_name, 'w') as resf:
        resf.write('NATRO\nstart\n')
        for res in range(first_res_pose, last_res_pose + 1):
            resf.write('{res} {ch} PIKAA {resn}   EX 1 EX ARO 2 USE_INPUT_SC\n'.format(res=str(pdb_inf.number(res)),
                                                                                   ch=chain_name,
                                                                                   resn=pose.residue(res).name1()))
    return chain_length, 'resfile_%s'%pdb_name



def split_pep_to_windows_of_template_length(peplist, window_size):
    """Split the peptide/s to overlapping windows of length n"""
    with open(peplist, 'r') as peplist:
        peps = [p.strip() for p in peplist.readlines()]
    pep_list = ''
    for peptide in peps:
        if len(peptide) - window_size > 0:
            with open (peptide+'.list', 'w') as pep_list:
                for i in range(len(peptide) - window_size):
                    for j in range(i, i + window_size):
                        pep_list.write(peptide[j])
                    pep_list.write('\n')

    if pep_list:
        return True
    else:
        return False


def split_and_run(plist, resfile, pdb, bb_min, protocol, prefix, tasks, run_wins):
    with open(plist, 'r') as peps: # list of different peptides (substrates)
        peptides = [p.strip() for p in peps.readlines()]
    with open(resfile, 'r') as rf:
        resf_lines = rf.readlines()

    for pep in peptides: # for each peptide there is a list of windows of length n (can be n=pep_length )
        if not os.path.isdir(pep):
            os.mkdir(pep)
        if run_wins:
            win_list = pep + '.list'
            os.rename(win_list, os.path.join(pep, win_list))
            with open(os.path.join(pep, win_list), 'r') as p: # list of peptide windows
                wins = [win.strip() for win in p.readlines()]
            for win in wins:
                win_dir = os.path.join(pep, win)
                if not os.path.isdir(win_dir):
                    os.mkdir(win_dir)
                run_fpb(pdb, resf_lines, win_dir, win, bb_min, protocol, prefix, tasks)
        else:
            run_fpb(pdb, resf_lines, pep, pep, bb_min, protocol, prefix, tasks)


def run_fpb(pdb, resf_lines, cur_dir, pepseq, bb_min, protocol, prefix, tasks):
    home_dir = os.getcwd()
    header = resf_lines[:2]
    rf_lines = resf_lines[2:]
    with open(os.path.join(cur_dir, 'resfile.%s' % pepseq), 'w') as new_rf:
        for hline in header: # 2 lines here
            new_rf.write(hline)
        for i, line in enumerate(rf_lines):
            splt_line = line.split()
            splt_line[3] = pepseq[i]
            new_rf.write(' '.join(splt_line) + '\n')
    # write running batch script
    write_batch_script(pdb, 'resfile.%s' % pepseq, pepseq, cur_dir, bb_min, protocol, prefix, tasks)
    # run
    os.chdir(cur_dir)
    os.system('sbatch run_fpbind_picasso')
    os.chdir(home_dir)


def write_batch_script(template_name, resfile_name, pepseq, path, bb_min, protocol, prefix, tasks):
    with open(os.path.join(path, 'run_fpbind_picasso'),'w') as b_file:
        temp_name_noext = os.path.splitext(os.path.basename(template_name))[0]
        # header
        b_file.write('#!/bin/bash\n#SBATCH --ntasks={t}\n#SBATCH --time=50:00:00\n'
                     '#SBATCH --get-user-env\n#SBATCH --mem-per-cpu=1600m\n\n'.format(t=tasks))
        # rosetta directories
        b_file.write('ROSETTA_BIN="/vol/ek/share/rosetta/rosetta_src_2019.14.60699_bundle/main/source/bin"\n'
                     'ROSETTA_DB="/vol/ek/share/rosetta/rosetta_src_2019.14.60699_bundle/main/database/"\n\n')
        # fixbb design
        b_file.write('$ROSETTA_BIN/fixbb.linuxgccrelease -s {temp} -database $ROSETTA_DB '
                     '-flexPepDocking::receptor_chain A -resfile {resf} -scorefile design.score.sc > '
                     'design.log\n'.format(temp=template_name,resf=resfile_name))
        # fpd/minimization run
        b_file.write('$ROSETTA_BIN/FlexPepDocking.linuxgccrelease -database $ROSETTA_DB -s {start} -ex1 -ex2aro '
                     '{protocol} {bbmin} -flexPepDocking:flexpep_score_only -use_input_sc -out:prefix {pseq}_ '
                     '-scorefile {prefix}.score.sc > {prefix}.log'.format(start=temp_name_noext+'_0001.pdb',
                                                                          protocol=protocol, pseq=pepseq, prefix=prefix,
                                                                          bbmin=bb_min))


def main():

    args = arg_parser().parse_args()

    pdb = args.pdb
    pdb_path = os.path.abspath(pdb)

    chain_num = int(args.chain) # Default = 2

    temp_length, resfile_name = create_init_resfile(pdb, chain_num)

    pep_list = args.peplist

    if args.refine:
        protocol = '-flexPepDocking:pep_refine -nstruct 100'
        prefix = 'fpd'
        tasks = '30'
    else:
        protocol = '-flexPepDockingMinimizeOnly'
        prefix = 'min'
        tasks = '1'

    if args.bbmin:
        bb_min = '-min_receptor_bb'
    else:
        bb_min = ''

    mult_windows = split_pep_to_windows_of_template_length(pep_list, temp_length)

    split_and_run(pep_list, resfile_name, pdb_path, bb_min, protocol, prefix, tasks, mult_windows)

if __name__ == "__main__":

    pyrosetta.init()

    main()
