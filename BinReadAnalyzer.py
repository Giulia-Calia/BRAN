""" Bin read analyzer class"""

import pickle
import os
import sys 
import argparse
import pandas as pd
import numpy as np
import gzip
import plotly
import plotly.graph_objects as go
from Bio.SeqIO.FastaIO import SimpleFastaParser
from BinReadCounter import BinReadCounter
import rpy2.robjects as robjects
from rpy2.robjects.vectors import FactorVector, StrVector, IntVector
from rpy2.robjects.packages import importr

class BinReadAnalyzer:
    def __init__(self, folder_name, bin_size, ref_genome, flags, out_pickle):
        self.folder = folder_name
        self.bin_size = bin_size
        self.flags = flags 
        self.ref = ref_genome
        self.out = out_pickle
        # self.control_file = None

    def load_data(self):
        """If the pickle file is present and the parameters are not 
        changed, the data structures were kept otherwise it runs the 
        methods load_reads and loads_Ns of the imported class BinReadCounter
        """
        actual_params = [self.bin_size, self.flags, self.ref]
        equals = []
       
        list_dir = os.listdir(self.folder)
        # verify if a pickle file is present in the desired folder
        for file in list_dir:
            if file.endswith(".p"):
                # print(file)
                # if so, checks if the parameters for the building-up of
                # the previous pickle file, are the same as the actual ones
                counter = BinReadCounter(self.folder, self.bin_size, self.flags, self.ref, file)
                counter.set_out(file)
                old_params = counter._load_pickle()[0]
                # print(old_params)
                if actual_params == old_params:
                    # if the parameters are the same, the data strucutures were imported
                    equals.append(file)
                    print("No difference in parameters")
                    return counter._load_pickle()[1]

                else:
                    continue
            else:
                continue
        
        if not equals:
            # if the list is empty, none of the pickle files in the current directory
            # have the same parameters of the actual running thus the algorithm uses 
            # the modules of the BinReadCounter to calculate the new data structures 
            # and it saves a new pickle file 
            print("Parameters are changed\nBinReadCounter is running with actual parameters\nIT COULD TAKE A WHILE")
            counter._export_pickle()
            return counter.load_data()


    def plot_chrom_sample(self, chrom, sample, ns, fig=go.Figure()):
        read_df = self.load_data()[0]

        if not os.path.exists("plots"):
            os.mkdir("plots")
        
        fig.update_xaxes(title_text="Cromosomes_Bins")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        single_chrom = read_df[read_df["chr"] == "CH.chr" + str(chrom)]
        col_list = list(single_chrom.columns)
        for i in range(len(col_list[:col_list.index(sample)]) + 1):
            if col_list[i] == sample:
                fig.add_trace(go.Scatter(x=single_chrom["bin"], 
                                        y=single_chrom[sample], 
                                        mode="markers",
                                        name=str(col_list[i])))

        if ns:
            count_N_df = self.load_data()[1]
            single_chrom_N = count_N_df[count_N_df["chr"] == "CH.chr" + str(chrom)]
            fig.add_trace(go.Scatter(x=single_chrom_N["bin"],
                                    y=single_chrom_N["N_count"],
                                    mode="markers",
                                    name="N_counts"))

        fig.update_layout(title="Read Counts - Clone: " + sample + " - Chr: " + str(chrom) + " - Bin Size: " + str(self.bin_size))
        
        fig.show()
        

    def plot_chromosome(self, chrom, ns, fig=go.Figure()):
        read_df = self.load_data()[0]
        
        if not os.path.exists("plots"):
            os.mkdir("plots")
        
        fig.update_xaxes(title_text="Cromosomes_Bins")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        single_chrom = read_df[read_df["chr"] == "CH.chr" + str(chrom)]
        col_list = list(single_chrom.columns)
        for i in range(len(col_list)):
                if col_list[i] != "bin" and col_list[i] != "chr":
                    fig.add_trace(go.Scatter(x=single_chrom["bin"], 
                                            y=single_chrom[col_list[i]], 
                                            mode="markers",
                                            name=str(col_list[i])))

        if ns:
            count_N_df = self.load_data()[1]
            single_chrom_N = count_N_df[count_N_df["chr"] == "CH.chr" + str(chrom)]
            fig.add_trace(go.Scatter(x=single_chrom_N["bin"],
                                    y=single_chrom_N["N_count"],
                                    mode="markers",
                                    name="N_counts"))


        fig.update_layout(title="Read Counts - All Clones - Chr: " + str(chrom) + " - Bin Size: " + str(self.bin_size))
        fig.write_image("plots/counts_chr" + str(chrom) + "_all_"  + str(self.bin_size) + ".pdf")
        
        fig.show()
   

    def plot_sample(self, sample, ns, fig=go.Figure()):
        read_df = self.load_data()[0]
        
        if not os.path.exists("plots"):
            os.mkdir("plots")
        
        fig.update_xaxes(title_text="Cromosomes_Bins")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        col_read_df = list(read_df.columns)
        for i in range(len(col_read_df[:col_read_df.index(sample)]) + 1):
            if col_read_df[i] == sample:
                fig.add_trace(go.Scatter(x=read_df["bin"], 
                                        y=read_df[sample], 
                                        mode="markers",
                                        name=str(col_read_df[i])))

        if ns:
            count_N_df = self.load_data()[1]
            fig.add_trace(go.Scatter(x=count_N_df["bin"],
                                    y=count_N_df["N_count"],
                                    mode="markers",
                                    name="N_counts"))

        
        
        fig.update_layout(title="Read Counts - Clone: " + sample + " - Bin Size: " + str(self.bin_size))
        
        fig.show()


    def plot_all(self, ns):
        read_df = self.load_data()[0]
        
        if not os.path.exists("plots"):
            os.mkdir("plots")

        fig = go.Figure()
        fig.update_xaxes(title_text="Cromosomes_Bins")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        col_list = list(read_df.columns)
        for i in range(len(col_list)):
                if col_list[i] != "bin" and col_list[i] != "chr":
                    fig.add_trace(go.Scatter(x=read_df["bin"], 
                                            y=read_df[col_list[i]], 
                                            mode="markers",
                                            name=str(col_list[i])))

        if ns:
            count_N_df = self.load_data()[1]
            fig.add_trace(go.Scatter(x=count_N_df["bin"],
                                    y=count_N_df["N_count"],
                                    mode="markers",
                                    name="N_counts"))

        fig.update_layout(title="Read Counts - All Clones - All Chromosomes - Bin Size: " + str(self.bin_size))
        fig.write_image("plots/counts_all_" + str(self.bin_size) + ".pdf")

    
    def normalize_bins(self):
        # if self.control == True:
        # robjects.r("library(edgeR)")
        read_counts = self.load_data()[0]
        print(read_counts.sum(axis=0)[2])
        col_list = list(read_counts.columns)
        print(col_list)
        pass



        
        


    
    def plot_norm_data(self):
        pass


    def get_significant_diff(self, chr = None, sample = None):
        pass


    def plot_sig_data(self):
        pass


if __name__ == "__main__":
    # ad an argument to argparse, with the name of the file as control for the normalization

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
                        default="chardonnay_primary_contigs_chromosome-order.fa.gz",
                        help="""The path of the reference file, if not specified, the coverage is calculated with the mean of the mapping reads for all the samples""")

    parser.add_argument("-ch", "--chromosome",
                        default=None,
                        type=int,
                        help="""The number of the chromosome of interest for the plot of counts""")
                    
    parser.add_argument("-s", "--sample", 
                        default=None,
                        type=str,
                        help="""The name of the clone of interest for the plot of counts""")

    parser.add_argument("-N", "--Ns_count",
                        default=None,
                        type=str,
                        help="""Specify if the Ns counts has to be included in the plot of the read counts""")

    args = parser.parse_args()
    dict_args = vars(parser.parse_args([]))
    
    if args.flag_list != dict_args["flag_list"]:
        flags = dict_args["flag_list"] + args.flag_list
    else:
        flags = args.flag_list
    
    
    analyzer = BinReadAnalyzer(args.folder, args.bin_size, args.reference, flags, args.output_pickle)
   
    # analyzer.load_data()
    
    if args.chromosome and args.sample:
        if args.Ns_count:
            analyzer.plot_chrom_sample(args.chromosome, args.sample, ns=args.Ns_count)
        else:
            analyzer.plot_chrom_sample(args.chromosome, args.sample, ns=False)

    elif args.chromosome:
        if args.Ns_count:
            analyzer.plot_chromosome(args.chromosome, ns=args.Ns_count)
        else:
            analyzer.plot_chromosome(args.chromosome, ns=False)

    elif args.sample:
        if args.Ns_count:
            analyzer.plot_sample(args.sample, ns=args.Ns_count)
        else:
            analyzer.plot_sample(args.sample, ns=False)
    
    else:
        if args.Ns_count:
            analyzer.plot_all(ns=args.Ns_count)
        else:
            analyzer.plot_all(ns=False)
    

    analyzer.normalize_bins()