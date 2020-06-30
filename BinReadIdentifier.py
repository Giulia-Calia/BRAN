# testing BinReadIdentifier

# This class is aimed to retrieve the ID of reads present in the sig_bins file
# given a certain chromosome and a certain position on the chromosome, together with the bin_size,
# it return a new file in which the ID of reads in that bin are stored
# another file is also created with the IDs of all unmapped reads of the sample

# till now only for one interesting bin per chromosome can be retrieved the read IDs


import os
import pysam
import progressbar
# import pandas as pd
import argparse
import re  # regular expression
import csv


class BinReadIdentifier:
    def __init__(self, bin_size, flag_list, folder, saving_folder, bins, cigar, cigar_filter):
        self.bin_size = int(bin_size)
        self.flags = flag_list
        self.bam_folder = folder  # a list of folders
        self.saving_folder = saving_folder
        self.bins = bins  # a dictionary of chromosome:bin_str_pos pairs
        self.cigar = cigar
        self.cigar_filter = cigar_filter
        self.bam = None

    def set_bam_list(self, bam_list):
        self.bam = bam_list

    def load_bam(self):
        bam = []
        for f in self.bam_folder:
            dir_list = os.listdir(f)
            for el in dir_list:
                if el.endswith(".bam"):
                    bam.append(f + el)
        self.set_bam_list(bam)
        return self.bam

    def get_read_ids(self):
        header = "#  {} bin_size: {}" \
                 "\nID \t\t\t str_pos \t end_pos \t chr \t chr_bin_pos".format(self.bins, self.bin_size)
        for file in self.load_bam():
            bam_file = pysam.AlignmentFile(file, "rb")
            reads_bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
            update_bar = 0
            clone_name = file[file.rfind("/") + 1:file.find(".bam")]
            ids = open(self.saving_folder + "/read_ids_" + str(self.bin_size) + "_" + clone_name + ".tsv", "w")
            ids.write(header + "\n")
            print("\n", clone_name + "\t searching for read ID")
            read_ref_names = []
            # print(bam_file.references)
            for ref in bam_file.references:  # name of chr to which the read belong
                for ch in self.bins.keys():
                    if re.search(r'\b' + ch + r'\b', ref):
                        # the regular expression syntax allows to get only chr1 when specified because
                        # otherwise it would search also for chr10, chr11 ecc
                        read_ref_names.append(ref)
            # print(read_ref_names)
            for ch, pos, ref_names in zip(self.bins.keys(), self.bins.values(), read_ref_names):
                # print(ch)
                # print(pos)
                # print(ref_names)
                # print(len(self.bins[ch]))
                for bin_p in range(len(self.bins[ch])):
                    # print(type(self.bins[ch][bin_p]))
                    # print(pos)
                    for read in bam_file.fetch(contig=ref_names):
                        update_bar += 1
                        reads_bar.update(update_bar)
                        start_p = int(self.bins[ch][bin_p])
                        end_p = int(self.bins[ch][bin_p]) + self.bin_size
                        if str(read.flag) in self.flags and \
                                self.cigar and \
                                read.cigarstring is not None and \
                                start_p <= int(read.reference_start) < end_p:
                            ids.write(read.query_name + " \t " +
                                      str(read.reference_start) + " \t\t " +
                                      str(int(read.reference_start) + len(read.query_sequence) - 1) + "\t" +
                                      ref_names + "\t" +
                                      self.bins[ch][bin_p] + "\n")

                        elif str(read.flag) in self.flags and \
                                not self.cigar and \
                                read.cigarstring is not None and \
                                start_p <= int(read.reference_start) < end_p:
                            ids.write(read.query_name + " \t " +
                                      str(read.reference_start) + " \t " +
                                      str(int(read.reference_start) + len(read.query_sequence) - 1) + "\n")

                        elif read.is_unmapped:
                            # print(read)
                            ids.write(read.query_name + " \t - \t\t - \t\t - \t\t -\n")

                        else:
                            continue
                else:
                    continue

            ids.close()

    # def load_ids(self):
    #     self.load_bam()
    #     ids_files = []
    #     list_header = []
    #     for el in os.listdir(self.saving_folder):
    #         if el.startswith("read_ids"):
    #             ids_files.append(el)
    #     if len(ids_files) == len(self.bam):
    #         for file in ids_files:
    #             with open(self.saving_folder + file) as id_file:
    #                 read_tsv = csv.reader(id_file, delimiter="\t")
    #                 first_line = next(read_tsv)
    #                 # be careful because every time next(read_tsv) appear, the iterator goes ahead of one line
    #                 if str(self.bins) in first_line:
    #                     list_header.append(str(first_line))
    #                     # next return a list with the line as a single string so to be interpreted by the set(),
    #                     # each line-list has to be converted in a string
    #                 else:
    #                     continue
    #         else:
    #             pass
    #
    #     res = len(set(list_header)) == 1
    #     if res == 1:
    #         print("No need to recalculate ids, the files in current directory are what you're searching for")
    #     else:
    #         print("BRAN needs a while to recalculate the read IDs in the required chromosome and bin positions")
    #         self.get_read_ids()

    def old_mapped_ids(self):
        header = "chr \t ID \t\t clone_name \t str_pos \t end_pos  \t type"
        header_unmapped = "chr \t ID \t\t clone_name"
        for file in self.bam:
            bam_file = pysam.AlignmentFile(file, "rb")
            reads_bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
            update_bar = 0
            clone_name = file[file.rfind("/") + 1:file.find(".bam")]
            sig_ids = open(self.saving_folder + "/sig_read_ids_" + clone_name + ".tsv", "w")
            sig_ids.write(header + "\n")
            unmapped_ids = open(self.saving_folder + "/unmapped_ids_" + clone_name + ".tsv", "w")
            unmapped_ids.write((header_unmapped + "\n"))
            print("\n", clone_name)
            for ch, bin_str in zip(self.bins.keys(), self.bins.values()):
                for ref in bam_file.references:
                    if ch in ref:
                        for read in bam_file.fetch(contig=ref):
                            update_bar += 1
                            reads_bar.update(update_bar)
                            if read.is_unmapped:
                                # print(read)
                                unmapped_ids.write(ref + " \t " + read.query_name + " \t " + clone_name + "\n")

                            if str(read.flag) in self.flags and \
                                    self.cigar and \
                                    read.cigarstring is not None and \
                                    bin_str <= int(read.reference_start) < (bin_str + self.bin_size):
                                if not any(filters in read.cigarstring for filters in self.cigar_filter):
                                    sig_ids.write(ref + " \t " +
                                                  read.query_name + " \t " +
                                                  clone_name + " \t " +
                                                  str(read.reference_start) + " \t " +
                                                  str(int(read.reference_start) + len(read.query_sequence) - 1) +
                                                  " \t " + "properly_mapped" + "\n")
                                else:
                                    sig_ids.write(ref + " \t " +
                                                  read.query_name + " \t\t " +
                                                  clone_name + " \t " +
                                                  str(read.reference_start) + " \t " +
                                                  str(int(read.reference_start) + len(read.query_sequence) - 1) +
                                                  " \t " + "clipped" + "\n")
                            elif str(read.flag) in self.flags and \
                                    not self.cigar and \
                                    read.cigarstring is not None and \
                                    bin_str <= int(read.reference_start) < (bin_str + self.bin_size):
                                sig_ids.write(ref + " \t " +
                                              read.query_name + " \t " +
                                              clone_name + " \t " +
                                              str(read.reference_start) + " \t " +
                                              str(int(read.reference_start) + len(read.query_sequence) - 1) +
                                              " \t " + "mapped" + "\n")

            sig_ids.close()
            unmapped_ids.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     # formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="",
                                     epilog="")
    parser.add_argument("-bs", "--bin_size",
                        type=int,
                        default=250000,
                        help="The length of the segments that divide the chromosomes equally")

    parser.add_argument("-fl", "--flag_list",
                        nargs="+",
                        default=["99", "147", "163", "83"],
                        help="""A list of the bitwise-flags in SAM format that identify the reads to be counted 
                           during analyses; if different flags wants to be added, add them as strings 
                           (e.g. "177" "129")""")

    parser.add_argument("-f", "--folder",
                        nargs="+",
                        default=["./"],
                        help="The path to the folder in which are located the files to be analyzed (.bam)")

    parser.add_argument("-sf", "--saving_folder",
                        nargs="+",
                        default="./bins_read_id",
                        help="""The path to the folder in which save the .tsv files of read ids.
                            The program save a file for each .bam containing only one bin per chromosome""")

    parser.add_argument("-ch", "--chromosomes",
                        nargs="+",
                        default=[],
                        type=str,
                        help="""the name of the chromosome for each interesting bin (no repetitions)""")

    parser.add_argument("-bp", "--bin_positions",
                        nargs="+",
                        default=[],
                        type=int,
                        help="""The bin position on the corresponding chromosome, be careful that for each position 
                        there is one and only one chromosome""")

    parser.add_argument("-c", "--cigar",
                        action="store_true",
                        help="If specified, it allows the application of all the filters on cigar_string, per read")

    parser.add_argument("-cf", "--cigar_filter",
                        nargs="+",
                        default=["S", "H"],
                        help="""If specified, the reads mapped with soft and hard clipping (S and H) by default, 
                        are taken out form the read counts; it returns a data frame with same structure of the 
                        default one.\n(Specify other filters like e.g. "I" "D")""")

    parser.add_argument("-id", "--identifier",
                        # type=str,
                        action="store_true",
                        help="if identifier class is needed")

    # print(parser.print_help())
    args = parser.parse_args()

    dict_args = vars(parser.parse_args([]))

    if args.flag_list != dict_args["flag_list"]:
        flags = dict_args["flag_list"] + args.flag_list
    else:
        flags = args.flag_list

    if not os.path.exists(args.saving_folder):
        os.mkdir(args.saving_folder)

    bin_dictionary = dict(zip(args.chromosomes, args.bin_positions))

    ide = TestingBinReadIdentifier(args.bin_size,
                                   flags,
                                   args.folder,
                                   args.saving_folder,
                                   bin_dictionary,
                                   args.cigar,
                                   args.cigar_filter)
    ide.load_bam()
    ide.read_ids()
