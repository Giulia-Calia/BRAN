# A raw version of the bin read counter, this script is essentially used
# to verify that the original counter will work correctly
# if it will be possible, it will be truly a testing_code, with test done
# with unittest

import argparse
import os
import gzip
import pickle
import pandas as pd
import plotly
import plotly.graph_objs as go
import plotly.express as px
import plotly.figure_factory as ff
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from pysam import pysam
import time
import progressbar



class BinReadCounter:
    """This class creates two data structures: read_counts and N_counts


    Attributes:
        folder (str): the folder in which retrieve the alignment files
        bin_size (int): the number of bases that compose the bins
                        (default = 250'000)
        flags (list): list of string; each string is a bitwise-flag in 
                    SAM format
        ref (str): the name of the reference file
        out (str): the name of the output file in pickle format 
    """
    def __init__(self, folder, bin_size, flag_list, pck_output, ref):
        self.folder = folder 
        self.bin_size = bin_size
        self.flags = flag_list
        self.out = pck_output
        self.ref = ref
        
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
        """add to the defoult list of bitwise-flags, the new flags"""
        self.flags = self.flags + other_flag

    def get_ref_path(self):
        """return the path of the reference file"""
        return self.ref

    def set_ref(self, other_ref):
        """set the file used as a reference"""
        self.ref = other_ref

    def get_ref_name(self):
        """return the name of the reference file"""
        if "/" in self.ref:
            ref_name = self.ref[self.ref.rfind("/") + 1:]
            return ref_name
        else:
            return self.ref

    def get_out(self):
        """return the name of the pickle output file"""
        return self.out

    def set_out(self, other_out):
        """set the name of the pickle output file"""
        self.out = other_out

    def load_data(self):
        """A simple method to retrieve the data structures"""
        return self._old_load_reads()

    def _old_load_reads(self):
        """Gives a data structure that stores information about the
        clone, the chromosome of each clone and the count of mapping 
        reads in each bin of each chromosome (keeping trak of the bins)


        Args:
            folder (str): folder name 
            bin_size (int): n_bases composing each bin
            flags (list): list of bitwise-flags


        Return:
            reads_count (dict): dictionary containing each count of reads
                                for each bin, for each chromosome and for
                                each sample/clone
            max_bin_number (int): the max number of bin in a chromosome 
        """
        bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
        i = 0
        chrom_list = {}
        id_bin = []
        read_counts = {}
        read_pos = 0
        dir_list = os.listdir(self.folder)
        for el in dir_list:
            if el.endswith(".sam"):  
                with open(el, "r") as current_file:                    
                    clone = el[:el.find(".")]
                    read_counts[clone] = []
                    chrom_list[clone] = []
                    for line in current_file:
                        i += 1
                        bar.update(i)
                        if line.startswith("@"):
                            header = line.replace("\n", "").split("\t")
                            if header[0] == "@SQ":
                                header[1] = header[1].replace("SN:", "")
                                header[2] = header[2].replace("LN:", "")
                                chromosome = header[1]
                                length = header[2]
                                # Here a list of all zeros is created in order to obtain an index on which
                                # iterate subsequently during reads-per-bin counts
                                # reads_count[clone] = []
                                bins = int(length)//self.bin_size + 1
                                for i in range(bins):
                                    chrom_list[clone].append(chromosome)
                                    read_counts[clone].append(0) 
                        else:
                            read_line = line.replace("\n", "").split("\t")
                            # Controll that the flag is in the list specified by the user 
                            if read_line[1] in self.flags:
                                # If so, the flag is transformed in a bit-format
                                #
                                # Once the bit-format is determined, the 4th bit starting from the end of
                                # the string is checked --> this correspond to the information of reverse/
                                # forward strand --> if "1", the read is in reverse strand and this imply
                                # that the position has to be calculated, otherwise the information available
                                # at the 3rd element in the file line, remains unchanged 
                                bit_flag = bin(int(read_line[1]))
                                if bit_flag[-4] == "1":
                                    read_pos = int(read_line[3]) + len(read_line[9]) - 1
                                else:
                                    read_pos = int(read_line[3])
                                if read_line[2] in chrom_list[clone]:
                                    # Place the read in the right bin 
                                    #
                                    # Dividing the read position for the length of the bin, the exact bin in
                                    # which the read maps is obtained
                                    # 
                                    # The corresponding element of the list is set to 1 if no other reads
                                    # mapped in that bin before and incremented by 1 otherwise 
                                    bin_location = read_pos//self.bin_size
                                    if read_counts[clone][bin_location] == 0:
                                        read_counts[clone][bin_location] = 1
                                    else:
                                        read_counts[clone][bin_location] += 1
                                else:
                                    continue
        # in order to retrieve only one time the bins that are equal for all
        # the clones
   
        for i in range(len(chrom_list[list(chrom_list.keys())[0]])):
            id_bin.append(i)
        df_bins = pd.DataFrame({"bin" : id_bin})
        df_chrom = pd.DataFrame({"chr" : chrom_list[list(chrom_list.keys())[0]]})
        df_counts = pd.DataFrame(read_counts)

        read_count_df = pd.concat([df_chrom, df_bins, df_counts], axis=1)    
        return read_count_df

    def _old_get_read_ID(self):
        # ok for small files but it take an infinity of time for bigger files
        ids = set()
        chr_location = []
        mate_in_bin = []
        read_pos = 0
        dir_list = os.listdir(self.folder)
        for el in dir_list:
            if el.endswith(".sam"):        
                with open(el, "r") as current_file:
                    for line in current_file:
                        if not line.startswith("@"):
                            read = line.replace("\n", "").split("\t")
                            if read[1] in self.flags:
                                read_id = read[0]
                                read_chr_location = read[2]
                                if read_id not in ids:
                                    ids.add(read_id)
                                    chr_location.append(read_chr_location)
                                    read_pos = int(read[3])//self.bin_size
                                    mate_pos = int(read[7])//self.bin_size
                                    if read_pos == mate_pos:
                                        mate_in_bin.append("Y")
                                    else:
                                        mate_in_bin.append("N")                        
                                                            
        read_id_df = pd.DataFrame({"ids" : list(ids), "chr" : chr_location, "mate_in_bin" : mate_in_bin}) 
        # print(len(ids), len(chr_location), len(mate_in_bin))
        # print(len(read_id_df[read_id_df["mate_in_bin"] == "N"]))
        return read_id_df

    def _old_load_Ns(self,  ref=None):
            """Gives a data structure that store information about the number
            of 'N' bases in each bin (the same bin of the reads_count), per
            chromosome.


            Args:
                ref (str): the name of the reference file
                bin_size (int): number of bases per each bin


            Return:
                n_per_bin (dict): a dictionary containing the count of 'N's
                                per each chromosome, per each bin
            """
            if ref != None:
                with gzip.open(ref, "rt") as reference:
                    chromosomes = []
                    bins = []
                    n_per_bin = []
                    # SimpleFastaParser is usefull for big fata files, because is
                    # a build-in method, faster than SeqIO.parse
                    for record in SimpleFastaParser(reference):
                        # record[0] = chrmosome name
                        # record[1] = sequence
                        # dividing each sequence in the same bins as the previous method
                        split_seq = [record[1][x : self.bin_size] if x == 0 else record[1][x : x + self.bin_size] for x in range(0, len(record[1]), self.bin_size + 1)]
                        for i in range(len(split_seq)):
                            chromosomes.append(record[0])
                            bins.append(i)
                            n_per_bin += [split_seq[i].count("N")]

                n_count = pd.DataFrame({"chr" : chromosomes, "bin" : bins, "N_count" : n_per_bin})

                return n_count

    def _old_export_pickle(self, out=None):
        """Export a pickle file containing all the parameters used for
        the counts and the completed data structures

        Return:
            pickle_file (file): file in pickle format (binary) that can be
                                imported in every script to dealing with the
                                data structures
        """
        if out != None:
            with open(self.out, "wb") as exp_file:
                # passing a list of elements to pickle, it writes them in
                # the same order in which they are specified
                pickle.dump([self.get_bin_size(), self.get_flags(), self.get_ref_name(), self.get_out(), self._load_reads(),
                             self._load_Ns()], exp_file)

    def _load_reads(self):
        """Gives a data structure that stores information about the
        clone, the chromosome of each clone and the count of mapping
        reads in each bin of each chromosome (keeping trak of the bins)

        Return:
            reads_count (dict): dictionary containing each count of reads
                                for each bin, for each chromosome and for
                                each sample/clone
            max_bin_number (int): the max number of bin in a chromosome
        """
        # the progress bar is used in this method, as in those after this, just
        # to check if the process is going on well or not
        bar_reads = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
        i = 0
        list_chrom = []
        chrom_column = {}
        bin_column = {}
        read_counts = {}
        read_counts_concat = {}
        read_pos = 0
        dir_list = os.listdir(self.folder)
        for el in dir_list:
            if el.endswith(".bam"):
                # here pysam is used in order to be able to work on .bam files,
                # that have a lower weight with respect to .sam files; thus the
                # process is faster
                samfile = pysam.AlignmentFile(el, "rb")
                clone = el[:el.find(".")]  # name of sample columns
                # clone_chrom_list["chr"] = []
                read_counts[clone] = []
                read_counts_concat[clone] = []
                chrom_column[clone] = []
                bin_column[clone] = []
                header = samfile.header["SQ"]

                for line in header:
                    chr_name = line["SN"]
                    chr_length = line["LN"]
                    if chr_name not in list_chrom:
                        list_chrom.append(chr_name)
                    bins = chr_length // self.bin_size + 1
                    # read_counts[clone].append([0] * bins)
                    read_counts[clone].append([0] * bins)
                # print(read_counts)
                    for i in range(bins):
                        chrom_column[clone].append(chr_name)
                        bin_column[clone].append(str(i))
                for chrom in list_chrom:
                    # for each chromosome in the file = each chromosome in the dictionary,
                    # at position = bin_location in the list of the chromosome, the element
                    # increments of one
                    for read in samfile.fetch(chrom):
                        i += 1
                        bar_reads.update(i)
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

                            # here i have to check for the chromosome, otherwise the count goes always to the first bins
                            bin_location = read_pos // self.bin_size
                            # print(bin_location)
                            read_counts[clone][list_chrom.index(chrom)][int(bin_location)] += 1

                        else:
                            continue

                for counts in read_counts[clone]:
                    read_counts_concat[clone] += counts
                samfile.close()
        # print(read_counts_concat)
        # print(clone_chrom_list)
        # print(read_counts)

        read_counts_concat_df = pd.DataFrame(read_counts_concat)
        chrom_column_df = pd.DataFrame({"chr": chrom_column[list(chrom_column.keys())[0]]})
        bin_column_df = pd.DataFrame({"bin": bin_column[list(bin_column.keys())[0]]})
        # df_chr = pd.DataFrame(clone_chrom_list)
        # df_bin = pd.DataFrame(id_bin)
        # df_chr_bin = pd.concat([df_chr, df_bin], axis=1)
        # print(id_bin)
        #
        # for i in range(len(clone_chrom_list[list(clone_chrom_list.keys())[0]])):
        #     # to append to the columns a progressive number of bins
        #     id_bin.append(i)
        # df_bins = pd.DataFrame({"bin": id_bin})
        # df_chrom = pd.DataFrame({"chr": clone_chrom_list[list(clone_chrom_list.keys())[0]]})
        # df_counts = pd.DataFrame(read_counts)
        #
        read_count_df = pd.concat([chrom_column_df, bin_column_df, read_counts_concat_df], axis=1)
        print(read_count_df)

        # return read_count_df

    def _load_Ns(self):
        """Gives a data structure that store information about the number
        of 'N' bases in each bin (the same bin of the reads_count), per
        chromosome.

        Return:
            n_per_bin (dict): a dictionary containing the count of 'N's
                            per each chromosome, per each bin
        """
        bar_Ns = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
        i = 0
        ref_name = self.get_ref_name()
        if ref_name.endswith(".fa.gz"):
            with gzip.open(self.ref, "rt") as reference:
                chromosomes = []
                bins = []
                n_per_bin = []
                # SimpleFastaParser is useful for big fasta files, because is
                # a build-in method, faster than SeqIO.parse
                for record in SimpleFastaParser(reference):
                    # record[0] = chromosome name
                    # record[1] = sequence
                    # dividing each sequence in the same bins as the previous method
                    i += 1
                    bar_Ns.update(i)
                    split_seq = [record[1][x: self.bin_size] if x == 0 else record[1][x: x + self.bin_size]
                                 for x in range(0, len(record[1]), self.bin_size + 1)]
                    for i in range(len(split_seq)):
                        chromosomes.append(record[0])
                        bins.append(i)
                        n_per_bin += [split_seq[i].count("N")]

        n_count = pd.DataFrame({"chr": chromosomes, "bin": bins, "N_count": n_per_bin})

        return n_count

    def _export_pickle(self):
        """Export a pickle file containing all the parameters used for
        the counts and the completed data structures

        Return:
            pickle_file (file): file in pickle format (binary) that can be
                                imported in every script to dealing with the
                                data structures
        """
        # pickle name is a combination of the various parameters:
        # BRAN +
        # bin_size +
        # df (if the param is the standard flag list)/ mf otherwise +
        # the first two letters of the file containing the reference seq
        out_name = "BRAN" + str(self.bin_size)
        if self.get_flags() == ["0", "16", "99", "147", "163", "83"]:
            out_name += "_df_" + self.get_ref_name()[:2] + ".p"
        else:
            out_name += "_mf_" + self.get_ref_name()[:2] + ".p"
        self.set_out(self.out + out_name)
        with open(self.out, "wb") as exp_file:
            # passing a list of elements to pickle, it writes them in
            # the same order in which they are specified
            pickle.dump([self.get_bin_size(), self.get_flags(), self.get_ref_name(), self._load_reads(), self._load_Ns()], exp_file)


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
                        type=str,
                        default="./",
                        help="The name of the path to the folder in which the files to be analyzed are situated")

    parser.add_argument("-fl", "--flag_list",
                        nargs="+",
                        default=["0","16","99","147","163","83"],
                        help="""A list of the bitwise-flags in SAM format that identify the reads to be caunted during analyses""")

    parser.add_argument("-op", "--output_pickle", 
                        type=str,
                        default="./",
                        help="""If specified creates the pickle file cointanining all the parameters and the 
                        data_structures""")

    parser.add_argument("-r", "--reference",
                        type=str,
                        default="/",
                        help="""The path of the reference file, if not specified, the coverage is calculated with the 
                        mean of the mapping reads for all the samples""")
                
    parser.add_argument("-ir", "--info_reads",
                        type=str,
                        default=False,
                        help="""If specified creates a data structure explaining if the read (read_ID) and its mate, 
                        are in the same bin or not""")

    args = parser.parse_args()
    dict_args = vars(parser.parse_args([]))
    if args.flag_list != dict_args["flag_list"]:
        flags = dict_args["flag_list"] + args.flag_list
    else:
        flags = args.flag_list
    
    counter = BinReadCounter(args.folder, args.bin_size, flags, args.output_pickle, args.reference)

    # with ChargingBar('Processing') as cbar:
    #     counter._old_load_reads()
    #     cbar.next()

    print(counter._load_reads())
  
    # counter._export_pickle(out="asd.p")
    # print(counter._get_read_ID())
    # counter.plot_hist()
    # print(counter._load_Ns(ref=args.reference))
