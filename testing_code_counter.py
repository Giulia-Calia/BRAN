# Bin Reads Counter Class
#
# This class aims to first divide the chromosomes into small pieces (bins),
# of about 0.25 Mb as default, and then count the number of desired reads
# that fall into each bin. This information is saved in a proper data
# structure that keeps track of the clone-name, the chromosome of interest
# and the bins for each chromosome.
#
# Due to the importance of "N" bases in the genome, also the count for these
# bases could be taken as a record for the counts-normalization, storing each
# count in another data structure that gives information on the chromosome
# and the count of "N" bases in each bin of it.
#
# The data structure as well as their parameters, are saved in a binary file
# created with pickle library, in order to facilitate future works on the
# structures.

import argparse
import os
import gzip
import pickle
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pysam
import progressbar


# ----------------------------------------------
# in readme insert progressbar2 not progressbar
# insert also psutils
# ----------------------------------------------

# ---------------------------------------------------
# see if it is possible to optimize the _load_reads()
# ---------------------------------------------------


class TestingBinReadCounter:
    """This class creates two data structures: read_counts and N_counts

    Attributes:
        folder (str): the folder in which retrieve the alignment files
        bin_size (int): the number of bases that compose the bins
                        (default = 250'000)
        flags (list): list of string; each string is a bitwise-flag in
                    SAM format
        ref (str): the path to the reference file
        out (str): the name of the output file in pickle format
    """

    def __init__(self, bam_folder_path, bin_size, flag_list, reference, out_pickle):
        self.folder = bam_folder_path
        self.bin_size = bin_size
        self.flags = flag_list
        self.ref = reference
        self.out = out_pickle
        self.cigar = None
        # self.read_counts = None
        # self.read_counts_cigar = None

    def get_folder(self):
        """return the folder name"""
        return self.folder

    def set_folder(self, other_folder):
        """set the folder to which retrieve the files"""
        self.folder = other_folder

    def get_bin_size(self):
        """return the number of bases in each bin"""
        return self.bin_size

    def set_bin_size(self, other_size):
        """set the number of bases in each bin"""
        self.bin_size = other_size

    def get_flags(self):
        """return the list of bitwise-flags that are filtered"""
        return self.flags

    def set_flags(self, other_flag):
        """add to the default list of bitwise-flags, the new flags"""
        self.flags = other_flag

    def get_ref_path(self):
        """return the path of the reference file"""
        return self.ref

    def set_ref(self, other_ref):
        """set the file used as a reference"""
        self.ref = other_ref

    def get_ref_name(self):
        """return the name of the reference file"""
        if self.ref is not None and "/" in self.ref:
            ref_name = self.ref[self.ref.rfind("/") + 1:]
            return ref_name
        else:
            return self.ref

    def get_bam_name(self):
        bam_list = []
        dir_list = os.listdir(self.folder)
        for el in dir_list:
            if el.endswith(".bam"):
                bam_list.append(el[:el.find(".bam")])
        return bam_list

    def get_out(self):
        """return the name of the pickle output file"""
        return self.out

    def set_out(self, other_out):
        """set the name of the pickle output file"""
        self.out = other_out

    def set_cigar(self, cigar_input):
        self.cigar = cigar_input

    # def set_cigar_reads(self, cigar_reads):
    #     self.read_counts_cigar = cigar_reads
    #
    # def set_read_counts(self, read_counts):
    #     self.read_counts = read_counts

    def pickle_file_name(self, other_cigar_filters, cigar_filter=False, reference=False, read_info=False):
        # pickle name is a combination of the various parameters:
        # BRAN +
        # bin_size +
        # df (if the param is the standard flag list)/ mf otherwise +
        # the first two letters of the file containing the reference seq
        #   if passed with parameters +
        # id if the read_id_info is required with parameters
        out_name = "BRAN" + str(self.bin_size)
        if self.get_flags() == ["0", "16", "99", "147", "163", "83"]:
            # if self.get_flags() == "3":
            out_name += "_df"
        else:
            out_name += "_mf"

        if cigar_filter and other_cigar_filters:
            out_name += "_cf"
        elif cigar_filter:
            out_name += "_c"

        if reference:
            out_name += "_" + self.get_ref_name()[:2]

        if read_info:
            out_name += "_id"

        out_name += ".p"

        return out_name

    # def load_data(self):
    #     """A simple method to retrieve the data structures"""
    def _load_cigar_read_counts(self, other_cigar_filters, cigar_filter=False):
        # the progress bar is used in this method, as in those after this, just
        # to check if the process is going on well or not
        if cigar_filter is True:
            cigar_clip = ["S", "H"]
            list_chrom = []
            index_column = []
            chrom_column = {}
            bin_column = {}
            read_counts = {}
            read_counts_concat = {}

            dir_list = os.listdir(self.folder)
            for el in dir_list:
                if el.endswith(".bam"):
                    # here pysam is used in order to be able to work on .bam files,
                    # that have a lower weight with respect to .sam files; thus the
                    # process is faster
                    samfile = pysam.AlignmentFile(self.folder + el, "rb")
                    print("\n" + el)
                    bar_reads = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
                    up = 0
                    header = samfile.header["SQ"]
                    clone = el[:el.find(".")]  # name of sample = columns
                    read_counts[clone] = []
                    read_counts_concat[clone] = []
                    chrom_column[clone] = []
                    bin_column[clone] = []

                    for line in header:
                        chr_name = line["SN"]
                        chr_length = line["LN"]
                        # a list of univocal chromosome names is created
                        if chr_name not in list_chrom:
                            list_chrom.append(chr_name)

                        bins = chr_length // self.bin_size + 1  # number of bin in each chromosome
                        read_counts[clone].append([0] * bins)

                        for i in range(bins):
                            chrom_column[clone].append(chr_name)
                            bin_column[clone].append(str(i))

                    # the list of univocal chromosome is used here in order to place the counts
                    # in the right chromosome
                    for chrom in list_chrom:
                        # for each chromosome in the file, at position = bin_location in the list
                        # of the chromosome, the element increments of one
                        for read in samfile.fetch(chrom):
                            # progress bar updating
                            up += 1
                            bar_reads.update(up)
                            if other_cigar_filters:
                                if read.cigarstring is not None:
                                    if not any(clip in read.cigarstring for clip in cigar_clip) and \
                                            not any(el in read.cigarstring for el in other_cigar_filters) and \
                                            str(read.flag) in self.flags:
                                        bit_flag = bin(int(read.flag))  # binary format more easy to check
                                        if bit_flag[-4] == "1":
                                            read_pos = int(read.reference_start) + len(read.query_sequence) - 1
                                        else:
                                            read_pos = int(read.reference_start)
                                        # Place the read in the right bin
                                        #
                                        # Dividing the read position for the length of the bin, the exact bin in
                                        # which the read maps is obtained
                                        #
                                        # The corresponding element of the list is set to 1 if no other reads
                                        # mapped in that bin before and incremented by 1 otherwise
                                        bin_location = read_pos // self.bin_size
                                        read_counts[clone][list_chrom.index(chrom)][int(bin_location)] += 1

                                    else:
                                        continue
                            else:
                                if read.cigarstring is not None:
                                    # if "S" not in read.cigarstring and "H" not in read.cigarstring:
                                    if not any(clip in read.cigarstring for clip in cigar_clip):
                                        if str(read.flag) in self.flags:
                                            bit_flag = bin(int(read.flag))  # binary format more easy to check
                                            if bit_flag[-4] == "1":
                                                read_pos = int(read.reference_start) + len(read.query_sequence) - 1
                                            else:
                                                read_pos = int(read.reference_start)
                                            # Place the read in the right bin
                                            #
                                            # Dividing the read position for the length of the bin, the exact bin in
                                            # which the read maps is obtained
                                            #
                                            # The corresponding element of the list is set to 1 if no other reads
                                            # mapped in that bin before and incremented by 1 otherwise
                                            bin_location = read_pos // self.bin_size
                                            read_counts[clone][list_chrom.index(chrom)][int(bin_location)] += 1

                                        else:
                                            continue

                    # in order to create a DataFrame, the lists of counts in read_counts
                    # have to be merged in a single list
                    for counts in read_counts[clone]:
                        read_counts_concat[clone] += counts

                    samfile.close()

            for j in range(len(bin_column[list(bin_column.keys())[0]])):
                index_column.append(j)
            # preparing for final DataFrame concatenation
            index_column_df = pd.DataFrame({"index": index_column})
            chrom_column_df = pd.DataFrame({"chr": chrom_column[list(chrom_column.keys())[0]]})
            bin_column_df = pd.DataFrame({"bin": bin_column[list(bin_column.keys())[0]]})
            read_counts_concat_df = pd.DataFrame(read_counts_concat)

            read_count_df = pd.concat([index_column_df, chrom_column_df, bin_column_df, read_counts_concat_df],
                                      axis=1)
            # self.set_cigar_reads(read_count_df)
            # return self.read_counts_cigar
            return read_count_df

    def _load_reads(self):
        """Gives a data structure that stores information about the
        clone, the chromosome of each clone and the count of mapping
        reads in each bin of the chromosome (keeping track of the bins)

        Return:
            reads_count_df (DataFrame): pandas DataFrame containing counts of reads
                                for each bin, for each chromosome and for
                                each sample/clone
        """
        # the progress bar is used in this method, as in those after this, just
        # to check if the process is going on well or not
        list_chrom = []
        index_column = []
        chrom_column = {}
        bin_column = {}
        read_counts = {}
        read_counts_concat = {}
        count_none = 0

        dir_list = os.listdir(self.folder)
        for el in dir_list:
            if el.endswith(".bam"):
                # here pysam is used in order to be able to work on .bam files,
                # that have a lower weight with respect to .sam files; thus the
                # process is faster
                samfile = pysam.AlignmentFile(self.folder + el, "rb")
                print("\n" + el)
                bar_reads = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
                up = 0
                header = samfile.header["SQ"]
                clone = el[:el.find(".")]  # name of sample = columns
                read_counts[clone] = []
                read_counts_concat[clone] = []
                chrom_column[clone] = []
                bin_column[clone] = []

                for line in header:
                    chr_name = line["SN"]
                    chr_length = line["LN"]
                    # a list of univocal chromosome names is created
                    if chr_name not in list_chrom:
                        list_chrom.append(chr_name)

                    bins = chr_length // self.bin_size + 1  # number of bin in each chromosome
                    read_counts[clone].append([0] * bins)

                    for i in range(bins):
                        chrom_column[clone].append(chr_name)
                        bin_column[clone].append(str(i))

                # the list of univocal chromosome is used here in order to place the counts
                # in the right chromosome
                for chrom in list_chrom:
                    # for each chromosome in the file, at position = bin_location in the list
                    # of the chromosome, the element increments of one
                    for read in samfile.fetch(chrom):
                        # progress bar updating
                        up += 1
                        bar_reads.update(up)

                        if str(read.flag) in self.flags:
                            bit_flag = bin(int(read.flag))  # binary format more easy to check
                            if bit_flag[-4] == "1":
                                read_pos = int(read.reference_start) + len(read.query_sequence) - 1
                            else:
                                read_pos = int(read.reference_start)
                            # Place the read in the right bin
                            #
                            # Dividing the read position for the length of the bin, the exact bin in
                            # which the read maps is obtained
                            #
                            # The corresponding element of the list is set to 1 if no other reads
                            # mapped in that bin before and incremented by 1 otherwise
                            bin_location = read_pos // self.bin_size
                            read_counts[clone][list_chrom.index(chrom)][int(bin_location)] += 1

                        else:
                            continue
                # in order to create a DataFrame, the lists of counts in read_counts
                # have to be merged in a single list
                for counts in read_counts[clone]:
                    read_counts_concat[clone] += counts

                samfile.close()

        for j in range(len(bin_column[list(bin_column.keys())[0]])):
            index_column.append(j)
        # preparing for final DataFrame concatenation
        index_column_df = pd.DataFrame({"index": index_column})
        chrom_column_df = pd.DataFrame({"chr": chrom_column[list(chrom_column.keys())[0]]})
        bin_column_df = pd.DataFrame({"bin": bin_column[list(bin_column.keys())[0]]})
        read_counts_concat_df = pd.DataFrame(read_counts_concat)

        read_count_df = pd.concat([index_column_df, chrom_column_df, bin_column_df, read_counts_concat_df], axis=1)
        # self.set_read_counts(read_count_df)
        # return self.read_counts
        return read_count_df

    def _load_Ns(self):
        """Gives a data structure that store information about the number
        of 'N' bases in each bin (the same bin of the reads_count), per
        chromosome.

        Return:
            n_per_bin (DataFrame): a pandas DataFrame containing the count of 'N's
                            per each chromosome, per each bin
        """
        if self.get_ref_name().endswith(".fa.gz"):
            with gzip.open(self.get_ref_path(), "rt") as reference:
                print("\nN-counts; used reference is ", self.get_ref_name())
                bar_Ns = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
                i = 0
                index_n = []
                chromosomes = []
                bins = []
                n_per_bin = []
                # SimpleFastaParser is useful for big fasta files, because is
                # a build-in method, faster than SeqIO.parse
                for record in SimpleFastaParser(reference):
                    # record[0] = chromosome name
                    # record[1] = sequence
                    # dividing each sequence in the same bins as the previous method

                    # progress bar updating
                    i += 1
                    bar_Ns.update(i)
                    split_seq = [record[1][x: self.bin_size] if x == 0 else record[1][x: x + self.bin_size]
                                 for x in range(0, len(record[1]), self.bin_size + 1)]
                    for i in range(len(split_seq)):
                        chromosomes.append(record[0])
                        bins.append(i)
                        n_per_bin += [split_seq[i].count("N")]

                for j in range(len(chromosomes)):
                    index_n.append(j)
        n_counts = pd.DataFrame({"index": index_n, "chr": chromosomes, "bin": bins, "N_count": n_per_bin})
        return n_counts

    def _load_read_ID(self):
        """Gives information on the effective presence of the read mate in the
        same bin or not

        Return:
            A pandas DataFrame with three columns, id, chr, mate_in_same_bin
        """
        ids = set()
        chr_location = []
        mate_in_bin = []
        dir_list = os.listdir(self.folder)
        for el in dir_list:
            if el.endswith(".bam"):
                print("\nRead Bin Information for: ", el)
                bar_id = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
                i = 0
                samfile = pysam.AlignmentFile(self.folder + el, "rb")
                for read in samfile.fetch():
                    # progress bar updating
                    i += 1
                    bar_id.update(i)

                    if str(read.flag) in self.flags:
                        read_id = read.query_name
                        read_chr_location = read.reference_name
                        if read_id not in ids:
                            ids.add(read_id)
                            chr_location.append(read_chr_location)
                            # calculate the read position and the beginning of the mate
                            # if both fall into the same bin, an "Y" is assigned and an
                            # "N" otherwise
                            read_pos = int(read.reference_start) // self.bin_size
                            mate_pos = int(read.pnext) // self.bin_size
                            if read_pos == mate_pos:
                                mate_in_bin.append("Y")
                            else:
                                mate_in_bin.append("N")

        read_id_df = pd.DataFrame({"ids": list(ids), "chr": chr_location, "mate_in_bin": mate_in_bin})

        return read_id_df

    def _load_pickle(self, file_name):
        """Retrieve information from the pickle file

        Args:
            file_name (str): the name of the pickle file that has to be imported

        Return:
            out_data (dict): a dictionary containing all the parameters and the data
                            structures required by the user
        """
        with open(self.out + file_name, "rb") as input_param:
            return pickle.load(input_param)

    def _export_pickle(self, other_cigar_filters, cigar=False, reference=False, read_info=False):
        """Export a pickle file containing all the parameters used for
        the counts and the completed data structures

        Args:
            reference (bool): a parameter passed by the user, if true, the n_counts
                              data-structure is added to the dictionary
            read_info (bool): a parameter passed by the user, if true, the n_counts
                              data-structure is added to the dictionary
        Return:
            pickle_file (file): file in pickle format (binary) that can be
                                imported in every script to dealing with the
                                data structures
        """
        with open(self.out + self.pickle_file_name(other_cigar_filters, cigar, reference, read_info), "wb") as exp_file:
            # passing a python object to pickle, it writes the elements in
            # the same order in which they are specified
            out_data = {"bin_size": self.get_bin_size(),
                        "flags": self.get_flags(),
                        "bam": self.get_bam_name(),
                        "cigar_filt": None,
                        "other_cigar_filt": [],
                        "read_counts": None,
                        "ref": None,
                        "n_counts": None,
                        "info": None,
                        "read_id_info": None}

            if reference:
                out_data["ref"] = self.get_ref_name()
                out_data["n_counts"] = self._load_Ns()

            if read_info:
                out_data["info"] = True
                out_data["read_id_info"] = self._load_read_ID()

            if cigar and other_cigar_filters:
                out_data["cigar_filt"] = True
                out_data["other_cigar_filt"] = True
                out_data["read_counts"] = self._load_cigar_read_counts(other_cigar_filters, cigar)
            elif cigar:
                out_data["cigar_filt"] = True
                out_data["read_counts"] = self._load_cigar_read_counts(other_cigar_filters, cigar)
            else:
                out_data["read_counts"] = self._load_reads()

            print(exp_file, "\n", out_data)
            pickle.dump(out_data, exp_file)

    def stats(self):
        """Return:
                a DataFrame containing information on the total number of reads for each sample"""
        # we want verify if the number of reads is correctly calculated
        # number of total reads in the sample
        # number of reads that has flag 99 that has to be equal to reads with flag 147
        # number of reads that has flag 163 that has to be equal to reads with flag 83

        # count_reads_flag = [0] * len(self.flags)  # per each flag in the list the amount of reads in the file is
        # calculated
        read_counts = self._load_reads()
        for col in read_counts.columns:
            if col != "index" and col != "chr" and col != "bin":
                sum_col = sum(read_counts[col])
                print(col, " : ", sum_col, "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="",
                                     epilog="")

    parser.add_argument("-bs", "--bin_size",
                        type=int,
                        default=250000,
                        help="The length of the segments that divide the chromosomes equally")

    parser.add_argument("-f", "--folder",
                        default="./",
                        help="The path to the folder in which are located the files to be analyzed (.bam)")

    parser.add_argument("-fl", "--flag_list",
                        nargs="+",
                        default=["0", "16", "99", "147", "163", "83"],
                        # default="3",
                        help="""A list of the bitwise-flags in SAM format that identify the reads to be counted 
                        during analyses; if different flags wants to be added, add them as single strings
                        \n(Specify it like e.g. "177" "129")""")

    parser.add_argument("-op", "--output_pickle",
                        default="./",
                        help="""Path to the folder where the pickle file containing all the parameters and the 
                        data_structures has to be created; if not specified, the data structure are only displayed in 
                        current terminal and no pickle file is created""")

    parser.add_argument("-r", "--reference",
                        type=str,
                        default=None,
                        help="""The path to the reference file; if specified a data structures containing counts of 
                        N bases in the genome, for each bin, is built-up""")

    parser.add_argument("-i", "--read_info",
                        action="store_true",
                        help="""If specified, a data-frame with information on the read ID, and if the read 
                        and its mate map in the same bin in the same chromosome, is created""")

    parser.add_argument("-c", "--cigar_filter",
                        action="store_true",
                        help="""If specified, the reads mapped with soft and hard clipping (S and H), are taken out 
                        form the read counts; it returns a data frame with same structure of the default one""")

    parser.add_argument("-cf", "--other_cigar_filters",
                        type=str,
                        nargs="+",
                        default=[],
                        help="""An additional parameter to exclude other reads from the count, on the bases of other 
                        information in their cigar, like indels.\n(Specify it like e.g. "I" "D")""")

    args = parser.parse_args()
    dict_args = vars(parser.parse_args([]))

    if args.flag_list != dict_args["flag_list"]:
        flags = dict_args["flag_list"] + args.flag_list
        print(flags)
    else:
        flags = args.flag_list

    counter = TestingBinReadCounter(args.folder, args.bin_size, flags, args.reference, args.output_pickle)

    # if args.cigar:
    #     print(counter._load_cigar_read_counts(args.other_cigar_filters, args.cigar))
    # else:
    #     print(counter._load_reads())

    counter._export_pickle(args.other_cigar_filters, args.cigar_filter, args.reference, args.read_info)
    # print(counter.stats())
    # if args.reference and args.read_info:
    #     counter._export_pickle(reference=True, read_info=True)
    #
    # elif args.reference:
    #     counter._export_pickle(reference=True)
    #
    # elif args.read_info:
    #     counter._export_pickle(read_info=True)
    #
    # else:
    #     counter._export_pickle()

    # ----------------------------------------
    # see if it is possible to optimize the _load_reads()
    # -----------------------------------------
