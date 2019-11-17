# Bin Reads Counter Class
# 
# This class aims to first divide the chromosomes into small pieces (bins),
# of about 0.25 Mb as default, and then count the number of desired reads
# that fall into each bin. This information is saved in an appropiate data
# structure that keeps trak of the clone-name, the chromosome of interest
# and the bins for each chromosome.
# 
# Due to the importance of "N" bases in the genome, also the count for these
# bases is take as a record for the counts-normalization, storing each count
# in another data structure that gives inforamtion on the chromosome and the 
# count of "N" bases in each bin of it.
# 
# The data structure as well as their paramerers, are saved in a binary file
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
    def __init__(self, folder, bin_size, flag_list, reference, out_pickle):
        self.folder = folder 
        self.bin_size = bin_size
        self.flags = flag_list
        self.ref = reference 
        self.out = out_pickle

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
        self.flags = self.flags + (other_flag)


    def get_ref(self):
        """return the name of the reference file"""
        return self.ref


    def set_ref(self, other_ref):
        """set the file used as a reference"""
        self.ref = other_ref
    

    def get_out(self):
        """return the name of the pickle output file"""
        return self.out


    def set_out(self, other_out):
        """set the name of the pickle output file"""
        self.out = other_out


    def load_data(self):
        """A simple method to retrieve the data structures"""
        return self._load_reads(), self._load_Ns()


    def _load_reads(self):
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
        # the progress bar is used in this method as in those after this, just
        # to check if the process is going on well or not 
        bar_reads = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
        i = 0
        clone_chrom_list = {}
        id_bin = []
        read_counts = {}
        read_pos = 0
        dir_list = os.listdir(self.folder)
        for el in dir_list:
            if el.endswith(".bam"):
                # here pysam is used in order to be able to work on .bam files,
                # that have a lower weight with respect to .sam files; thus the
                # process is faster 
                samfile = pysam.AlignmentFile(el, "rb")
                clone = el[:el.find(".")]
                read_counts[clone] = []
                clone_chrom_list[clone] = []
                header = samfile.header["SQ"]
                    
                for el in header:
                    chr_name = el["SN"]
                    chr_length = el["LN"]
                    bins = chr_length//self.bin_size + 1
                    for i in range(bins):
                        clone_chrom_list[clone].append(chr_name)
                        read_counts[clone].append(0)

                for read in samfile.fetch():
                    i += 1
                    bar_reads.update(i)
                    if str(read.flag) in self.flags:
                        bit_flag = bin(int(read.flag))
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
                        bin_location = read_pos//self.bin_size
                        if read_counts[clone][bin_location] == 0:
                            read_counts[clone][bin_location] = 1
                        else:
                            read_counts[clone][bin_location] += 1
                    else:
                        continue
                    
                samfile.close()
            
        for i in range(len(clone_chrom_list[list(clone_chrom_list.keys())[0]])):
            id_bin.append(i)
        df_bins = pd.DataFrame({"bin" : id_bin})
        df_chrom = pd.DataFrame({"chr" : clone_chrom_list[list(clone_chrom_list.keys())[0]]})
        df_counts = pd.DataFrame(read_counts)

        read_count_df = pd.concat([df_chrom, df_bins, df_counts], axis=1)    
        return read_count_df
        
    
    def _load_Ns(self):
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
        bar_Ns = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
        i = 0 
        with gzip.open(self.ref, "rt") as reference:
            chromosomes = []
            bins = []
            n_per_bin = []
            # SimpleFastaParser is usefull for big fata files, because is
            # a build-in method, faster than SeqIO.parse 
            for record in SimpleFastaParser(reference):
                # record[0] = chrmosome name
                # record[1] = sequence   
                # dividing each sequence in the same bins as the previous method
                i += 1
                bar_Ns.update(i)
                split_seq = [record[1][x : self.bin_size] if x == 0 else record[1][x : x + self.bin_size] for x in range(0, len(record[1]), self.bin_size + 1)]
                for i in range(len(split_seq)):
                    chromosomes.append(record[0])
                    bins.append(i)
                    n_per_bin += [split_seq[i].count("N")]

        n_count = pd.DataFrame({"chr" : chromosomes, "bin" : bins, "N_count" : n_per_bin})
    
        return n_count
         

    def load_read_ID(self):  
        """Gives information on the effective presence of the read mate in the
        same bin or not 

        Args:
            folder (str): the name of the folder in which the .bam files are saved
            flags (list): list of string, where the string are the bitwise flags of
                        the reads of interest
            bin_size (int): the length of the bin

        Return:
            A pandas DataFrame with three columns, id, chr, mate_in_same_bin 
        """
        bar_id = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
        i = 0
        ids = set()
        chr_location = []
        mate_in_bin = []
        read_pos = 0
        dir_list = os.listdir(self.folder)
        for el in dir_list:
            if el.endswith(".bam"):   
                samfile = pysam.AlignmentFile(el, "rb")
                for read in samfile.fetch():
                    i += 1
                    bar_id.update(i)
                    if str(read.flag) in self.flags:
                        read_id = read.query_name
                        read_chr_location = read.reference_name
                        if read_id not in ids:
                            ids.add(read_id)
                            chr_location.append(read_chr_location)
                            # calculate the read position and the beginning of the mate 
                            # if both fall into the same bin, an "Y" is assigned 
                            read_pos = int(read.reference_start)//self.bin_size
                            mate_pos = int(read.pnext)//self.bin_size
                            if read_pos == mate_pos:
                                mate_in_bin.append("Y")
                            else:
                                mate_in_bin.append("N")
                                                                                    
        read_id_df = pd.DataFrame({"ids" : list(ids), "chr" : chr_location, "mate_in_bin" : mate_in_bin}) 
        print(len(ids), len(chr_location), len(mate_in_bin))
        print(len(read_id_df[read_id_df["mate_in_bin"] == "N"]))
        return read_id_df
            

    def _load_pickle(self):
        """Retrieve information from the pickle file


        Args:
            out (str): the name of the output file, in pickle format


        Return:
            parameters (list): list of all the parameter's value
            data_structures (list): list of the data structure created
                                    with this class 
        """
        with open(self.out, "rb") as input_param:
            # all the elements excluded the last two are the parameters used 
            # to build-up the data structures
            parameters = pickle.load(input_param)[:-2]
        with open(self.out, "rb") as input_data:
            # the last two elements are the two data structures 
            data_structures = pickle.load(input_data)[-2:]
        return parameters, data_structures


    def _export_pickle(self):
        """Export a pickle file containing all the parameters used for
        the counts and the completed data structures


        Args:
            out (str): the name of the output file, in pickle format
            _load_reads() (method): reads_count per bin, per chromosome,
                                    per clone
            _load_Ns() (method):'N' count per bin, per chromosome of the
                                reference genome
            bin_size (int): the number of bases that compose the bins
                            (default = 250'000)
            flags (list): list of string; each string is a bitwise-flag in 
                        SAM format
            folder (str): the folder in which retrieve the alignment files
            ref (str): the name of the reference file


        Return:
            pickle_file (file): file in pickle format (binary) that can be
                                imported in every script to dealing with the
                                data structures
        """                                     
        out_name = "BRAN" + str(self.bin_size) 
        if self.get_flags() == ["0","16","99","147","163","83"]:
            out_name += "_df_" + self.get_ref()[:2] + ".p"
        else:
            out_name += "_mf_" + self.get_ref()[:2] + ".p"
        self.set_out(out_name)
        with open(self.out, "wb") as exp_file:
            # passing a list of elements to pickle, it writes them in
            # the same order in which they are specified 
            pickle.dump([self.get_bin_size(), self.get_flags(), self.get_ref(), self._load_reads(), self._load_Ns()], exp_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                description="",
                                epilog="")

    parser.add_argument("-bs", "--bin_size",
                        type=int,
                        default=250000,
                        help="The lenght of the segments that divide the chromosomes equally")

    parser.add_argument("-f", "--folder",
                        type=str,
                        default="./",
                        help="The name of the path to the folder in which the files to be analyzed are situated")

    parser.add_argument("-fl", "--flag_list",
                        nargs="+",
                        default=["0","16","99","147","163","83"],
                        help="""A list of the bitwise-flags in SAM format that identify the reads to be caunted during analyses""")

    parser.add_argument("-op", "--output_pickle", 
                        default=None,
                        help="""If specified creates the pickle file cointanining all the parameters and the data_structures""")

    parser.add_argument("-r", "--reference",
                        type=str,
                        default=None,
                        help="""The path of the reference file, if not specified, the coverage is calculated with the mean of the mapping reads for all the samples""")

    parser.add_argument("-i", "--read_info",
                        type=str,
                        default=None,
                        help="""If specified, it returns a dataframe with information on the read ID, and if the read and its mate map in the same bin in the same chromosome""")

    args = parser.parse_args()
    dict_args = vars(parser.parse_args([]))
    
    if args.flag_list != dict_args["flag_list"]:
        flags = dict_args["flag_list"] + args.flag_list
        print(flags)
    else:
        flags = args.flag_list
    
    counter = BinReadCounter(args.folder, args.bin_size, flags, args.reference, args.output_pickle)
    counter.set_ref("chardonnay_primary_contigs_chromosome-order.fa.gz")

    if args.output_pickle:
        counter._export_pickle()
    else: 
        print(counter.load_data())
    
    if args.read_info:
        print(counter.load_read_ID())
    
    # print(counter._load_reads())
    # print(counter._load_Ns())
    # print(counter._load_pickle()[0])