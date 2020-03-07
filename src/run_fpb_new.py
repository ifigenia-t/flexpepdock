#!/usr/bin/env python3
from pyrosetta import *
from pyrosetta.rosetta import *
import os
import argparse

ROSETTA_BIN = (
    "/Users/itsitsa/Downloads/rosetta_bin_mac_2019.35.60890_bundle/main/source/bin"
)
ROSETTA_DB = (
    "/Users/itsitsa/Downloads/rosetta_bin_mac_2019.35.60890_bundle/main/database"
)
AA_LIST = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]
ALA = "A"
OUTPUT_DATA_FOLDER = "output_data"


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", "-p", dest="pdb")
    parser.add_argument("--chain_num", "-c", dest="chain", default=2)
    # Added peptide chain ID option
    parser.add_argument("--chain_id", "-cid", dest="chain_id")
    # Added receptor chain ID option
    parser.add_argument("--receptor_chain_id", "-rcid", dest="receptor_chain_id")
    parser.add_argument("--pep_list", "-l", dest="peplist")
    # Added boolean argument, generate pep list by alanine scanning
    parser.add_argument(
        "--alanine_scanning",
        "-ala",
        dest="alanine_list",
        action="store_true",
        default=False,
    )
    # Added boolean argument, generate peptide list by all by all scanning
    parser.add_argument(
        "--all_by_all_scanning",
        "-all",
        dest="all_by_all_list",
        action="store_true",
        default=False,
    )
    # Added single peptide input option
    parser.add_argument("--peptide", "-sp", dest="single_peptide")
    parser.add_argument(
        "--refinement", "-r", dest="refine", action="store_true", default=False
    )
    parser.add_argument(
        "--bb_min", "-b", dest="bbmin", action="store_true", default=False
    )
    return parser


def create_init_resfile(pdb_path, chain_id, output_path):
    """
    Create resfile with appropriate numbering with option PIKAA (sequence of the template)
    """
    pose = pose_from_pdb(pdb_path)
    pdb_inf = pose.pdb_info()
    chain = 0
    for chain_number in range(1, pose.num_chains() + 1):
        first_res_pose = pose.chain_begin(chain_number)
        chain_name = pdb_inf.chain(first_res_pose)
        print(chain_number, chain_name)
        if chain_id == chain_name:
            print(chain_number, chain_name)
            chain = chain_number
            break
    if chain == 0:
        return
    pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]

    first_res_pose = pose.chain_begin(chain)
    last_res_pose = pose.chain_end(chain)

    peptide_seq = pose.sequence(first_res_pose, last_res_pose)

    chain_name = pdb_inf.chain(first_res_pose)
    chain_length = last_res_pose - first_res_pose + 1
    peptide = []

    # create the output path
    os.makedirs(output_path, exist_ok=True)
    with open(os.path.join(output_path, "resfile_{}s".format(pdb_name)), "w") as resf:
        resf.write("NATRO\nstart\n")
        print("===> writing to results file:{}".format(resf.name))  # --- added
        for res in range(first_res_pose, last_res_pose + 1):
            peptide.append(str(pdb_inf.number(res)))
            resf.write(
                "{res} {ch} PIKAA {resn}   EX 1 EX ARO 2 USE_INPUT_SC\n".format(
                    res=str(pdb_inf.number(res)),
                    ch=chain_name,
                    resn=pose.residue(res).name1(),
                )
            )
        # print("===> peptide striped from structure:{}".format("".join(peptide)))  # --- added
    peptide = "".join(peptide)
    # print(peptide)
    return chain_length, "resfile_{}s".format(pdb_name), peptide, peptide_seq


def read_peptides_file(plist):
    with open(
        plist, "r"
    ) as peps:  # list of different peptides (substrates) from the file
        peptides = [p.strip() for p in peps.readlines()]
    return peptides


# # create the peptide list from one peptide
# def create_pep_list_file(peptide, amino_acid):
#     # Alanine Scanning
#     alanine_list = []
#     alanine_list.append(peptide)
#     for aa in range(0, len(peptide)):
#         if peptide[aa] != amino_acid:
#             new_peptide = peptide[:aa] + amino_acid + peptide[aa + 1 :]
#             alanine_list.append(new_peptide)
#     return alanine_list


def create_pep_list_file(peptide, amino_acids):
    # All by All Scanning
    all_by_all_list = []
    all_by_all_list.append(peptide)
    for i in amino_acids:
        for aa in range(0, len(peptide)):
            if peptide[aa] != i:
                new_peptide_all = peptide[:aa] + i + peptide[aa + 1 :]
                all_by_all_list.append(new_peptide_all)
    return all_by_all_list


def split_and_run(
    peptides,
    resfile,
    pdb,
    output_path,
    bb_min,
    protocol,
    prefix,
    tasks,
    receptor_chain_id,
):
    """

    """
    with open(os.path.join(output_path, resfile), "r") as rf:
        resf_lines = rf.readlines()

    for (
        pep
    ) in (
        peptides
    ):  # for each peptide there is a list of windows of length n (can be n=pep_length)
        print("==> running for peptides: {}".format(pep))

        os.makedirs(os.path.join(output_path, pep), exist_ok=True)
        run_fpb(
            pdb,
            output_path,
            resf_lines,
            pep,
            pep,
            bb_min,
            protocol,
            prefix,
            tasks,
            receptor_chain_id,
        )


def run_fpb(
    pdb,
    output_path,
    resf_lines,
    cur_dir,
    pepseq,
    bb_min,
    protocol,
    prefix,
    tasks,
    receptor_chain_id,
):
    """

    """
    # print("==> pepseq: {}".format(pepseq))  # -- added
    home_dir = os.getcwd()
    header = resf_lines[:2]
    rf_lines = resf_lines[2:]
    with open(os.path.join(output_path, cur_dir, "resfile.%s" % pepseq), "w") as new_rf:
        for hline in header:  # 2 lines here
            new_rf.write(hline)
        for i, line in enumerate(rf_lines):
            splt_line = line.split()
            splt_line[3] = pepseq[i]
            new_rf.write(" ".join(splt_line) + "\n")

    #
    # RUN fixBB - fixbb design
    #
    print(
        pdb,
        "resfile.{}s".format(pepseq),
        pepseq,
        cur_dir,
        bb_min,
        protocol,
        prefix,
        tasks,
    )

    current_dir_path = os.path.join(output_path, cur_dir)

    pdbfile = pdb  # --- added
    basefilename = os.path.basename(pdb)
    threaded_pdbfile_name = "{}_0001.pdb".format(os.path.splitext(basefilename)[0])
    threaded_pdbfile = os.path.join(
        current_dir_path, threaded_pdbfile_name,
    )  # --- added

    resfile = os.path.abspath(
        os.path.join(current_dir_path, "resfile.{}".format(pepseq))
    )  # --- added
    scorefile = os.path.abspath(
        os.path.join(current_dir_path, "scorefile.{}".format(pepseq))
    )  # --- added
    score = os.path.abspath(
        os.path.join(current_dir_path, "{}.{}.score.sc".format(pepseq, prefix))
    )
    fixbb_resfile_log = os.path.abspath(
        os.path.join(current_dir_path, "fixbb.log")
    )  # --- added
    flex_pep_dock_resfile_log = os.path.abspath(
        os.path.join(current_dir_path, "{}.flex_pep_dock.log".format(prefix))
    )

    fixbb_cmd = os.path.join(ROSETTA_BIN, "fixbb.macosclangrelease")
    fixbb_options = [
        "-s " + pdbfile,
        "-database " + ROSETTA_DB,
        "-flexPepDocking::receptor_chain "
        + receptor_chain_id,  # --- If not defined the receptor_chain is A
        "-resfile " + resfile,
        "-scorefile " + scorefile,
        "-overwrite",
    ]
    fixbb_run = fixbb_cmd + " " + " ".join(fixbb_options) + " > " + fixbb_resfile_log
    print("Running fixbb_run - " + os.getcwd())
    print(fixbb_run)
    os.system(fixbb_run)

    # TODO find the proper way to move the file to the right location
    os.rename(os.path.join(os.getcwd(), threaded_pdbfile_name), threaded_pdbfile)

    #
    # RUN FlexPepDocking - Minimization run
    #

    flex_pep_dock_cmd = os.path.join(ROSETTA_BIN, "FlexPepDocking.macosclangrelease ")
    flex_pep_dock_options = [
        "-s " + threaded_pdbfile,
        "-database " + ROSETTA_DB,
        "-ex1",
        "-ex2aro",
        "-flexPepDockingMinimizeOnly",
        "-use_input_sc",
        "-out:prefix {}_".format(os.path.join(current_dir_path, pepseq)),
        "-scorefile {}".format(score),
    ]
    flex_pep_dock_run = (
        flex_pep_dock_cmd
        + " "
        + " ".join(flex_pep_dock_options)
        + " > "
        + flex_pep_dock_resfile_log
    )
    print("Running flex_pep_dock - " + os.getcwd())
    print(flex_pep_dock_run)
    os.system(flex_pep_dock_run)
    os.chdir(home_dir)


def main():
    args = arg_parser().parse_args()
    pdb = args.pdb
    pdb_path = os.path.abspath(pdb)
    pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
    output_path = os.path.join(os.getcwd(), OUTPUT_DATA_FOLDER, pdb_name)

    # chain_num = int(args.chain) # Default = 2
    chain_id = args.chain_id
    receptor_chain_id = args.receptor_chain_id

    temp_length, resfile_name, peptide, peptide_seq = create_init_resfile(
        pdb, chain_id, output_path
    )

    print(
        "==> temp_lenght: {} resfile_name: {} structure_peptide: {}".format(
            temp_length, resfile_name, peptide
        )
    )

    if chain_id is None:
        print("chain ID is not defined")
        sys.exit(-1)

    if args.refine:
        protocol = "-flexPepDocking:pep_refine -nstruct 100"
        prefix = "fpd"
        tasks = "30"
    else:
        protocol = "-flexPepDockingMinimizeOnly"
        prefix = "min"
        tasks = "1"

    if args.bbmin:
        bb_min = "-min_receptor_bb"
    else:
        bb_min = ""

    peptides = []
    amino_acids = []

    if args.alanine_list is True:  # --- added
        # here we call the new function
        # and asign to peptides
        print(
            "===> calculating list from peptide sequence: {} by alanine scanning".format(
                peptide_seq
            )
        )
        amino_acids.append(ALA)

    if args.all_by_all_list is True:  # --- added
        print(
            "===> calculating list from peptide sequence: {} by all by all scanning".format(
                peptide_seq
            )
        )
        amino_acids = AA_LIST

    peptides = create_pep_list_file(peptide_seq, amino_acids)

    print(peptides)

    if args.single_peptide is not None:
        peptide = args.single_peptide
        print("===> using input peptide {}".format(peptide))
        peptides.append(peptide)

    if args.peplist is not None:
        pep_list = args.peplist
        print("===> using peptide file {}".format(pep_list))
        peptides = read_peptides_file(pep_list)

    # split_and_run(
    #     peptides,
    #     resfile_name,
    #     pdb_path,
    #     output_path,
    #     bb_min,
    #     protocol,
    #     prefix,
    #     tasks,
    #     receptor_chain_id,
    # )

    tbt_input_data = combine_score_files(peptides, output_path)

    create_tbt_file(peptide_seq, tbt_input_data, amino_acids, output_path)


# reads a score file and returns the variables we need
def parse_score_file(peptide, output_path):
    scorefile_path = os.path.join(
        output_path, peptide, "{}.min.score.sc".format(peptide)
    )

    with open(scorefile_path, "r") as f:
        line_counter = 0
        scores = {}
        headers = []

        for line in f:
            line_counter += 1
            if line_counter == 1:
                continue
            if line_counter == 2:
                headers = line.split()[1:]
                continue
            values = line.split()[1:]
            for i in range(0, len(values) - 1):
                scores[headers[i]] = values[i]

    return scores


# iterates over all the peptides and appends the output to one combined score file
# +++ and returns the input for the next file
def combine_score_files(peptides, output_path):
    tbt_reference = {}

    # create a new file
    with open(os.path.join(output_path, "finalScores.txt"), "wt") as f:
        f.write("Peptide       total_score     I_sc     pep_sc     reweighted_sc\n")
        # f.write("   {}")
        for peptide in peptides:
            scores = parse_score_file(peptide, output_path)
            f.write(
                "{}    {}    {}    {}   {}\n".format(
                    peptide,
                    scores["total_score"],
                    scores["I_sc"],  # Interface Surface
                    scores["pep_sc"],  # peptide Score
                    scores[
                        "reweighted_sc"
                    ],  # Reweighted Score  (sum(total_score+peptide_score+interface_score))
                )
            )
            tbt_reference[peptide] = scores["reweighted_sc"]

    return tbt_reference


def create_tbt_file(base_peptide, tbt_input_data, amino_acids, output_path):
    results = []
    with open(os.path.join(output_path, "reweightedScores.tbt"), "wt") as f:
        f.write("\t")
        for i in range(0, len(base_peptide)):
            f.write("{}\t".format(base_peptide[i]))
            print("{}\t".format(base_peptide[i]), end="")

        for aa in amino_acids:
            f.write("\n{}\t".format(aa))
            print(aa, end=" ")
            for i in range(0, len(base_peptide)):
                if base_peptide[i] != aa:
                    new_peptide = base_peptide[:i] + aa + base_peptide[i + 1 :]
                    f.write("{}\t".format(tbt_input_data[new_peptide]))
                    print("{}\t".format(tbt_input_data[new_peptide]), end="")
            print()


if __name__ == "__main__":
    pyrosetta.init()
    main()
