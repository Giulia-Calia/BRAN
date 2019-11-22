# BinReadAnalyzer class
#
# This class aims to analyze the data structure created with the imported
# class: BinReadCounter
# One can consider this class as the very effector of the analysis on the
# genome of the plant of interest. For analysis is intended, first of all,
# the normalization of the raw_read_counts and than the detection of
# statistically significant difference in terms of number of reads in bin,
# that could lead to the identification of part of chromosome or even an
# entire chromosome that have a significantly different count of reads;
# proving that something in that region or chromosome is happened.
#
# The script first checks if already exist, in the folder of interest, a
# pickle file; if so, the actual parameters and the parameters of the file
# are compared --> if they are the same, the method _load_pickle() from
# BinReadCounter is called and the data structures are available for the
# analysis, otherwise the _export_pickle() method of BinReadCounter is
# imported and the Counter starts to produce the data structures and the
# new pickle file
#
# After that, on the basis of the parameters given by the user, visualization
# plots are build up and shows in an as an interactive way via a web page;
# also an offline version of all the plots is saved in an appropriate folder

import os
import argparse
import plotly.graph_objects as go
from BinReadCounter import BinReadCounter
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
base = importr('base')
utils = importr('utils')
stats = importr('stats')


class BinReadAnalyzer:
    """This class provides different analysis of the data coming form
    BinReadCounter class; here imported.

    Options:
    -   plot read_counts for all chromosomes, for a specific one or for
        a specific sample
    -   normalize read_counts
    -   plot normalized_read_counts for all chromosomes, for a specific
        one or for a specific sample

    Attributes:
            folder_name: the folder in which the necessary files are stored
            bin_size: the length of the bin in which each chromosome is
                      divided (in bp)
            ref_genome: a string being the name of the reference_genome file
            flags: a list of string, each corresponding to a particular
                   pairwise-flag in the .bam/.sam format
            out_pickle: default is None because it takes the name of an
                        already existing file or the name of a new file
                        created during the running of the code
    """
    def __init__(self, folder_name, bin_size, ref_genome, flags, out_pickle):
        self.folder = folder_name
        self.bin_size = bin_size
        self.flags = flags 
        self.ref = ref_genome
        self.out = out_pickle
        # self.control_file = None

    def load_data(self, verbose=False):
        """If the pickle file is already present in the folder and the parameters
        are not changed, the data structures were kept, otherwise it runs the
        methods _export_pickle() and load_data() of the imported class BinReadCounter
        in order to create a new pickle file and to import the corresponding
        data-structures

        Args:
            verbose (bool): default is False, but if set to True, it allows to print
                            and see if a pickle file is already present or if the
                            Analyzer is running to create a new one

        Returns:
            A list of data structures from the pickle file or form the BinReadCounter;
            the first data structure is the read_counts and the second is the N bases
            count
        """
        counter = BinReadCounter(self.folder, self.bin_size, self.flags, self.ref, self.out)
        actual_params = [self.bin_size, self.flags, self.ref]
        equals = []
        new_file = ""
        list_dir = os.listdir(self.folder)
        # verify if a pickle file is present in the desired folder
        for file in list_dir:
            if file.endswith(".p"):
                # if so, checks if the parameters for the building-up of
                # the previous pickle file, are the same as the actual ones
                counter.set_out(file)
                old_params = counter._load_pickle()[0]
                # print(old_params)
                if actual_params == old_params:
                    # if the parameters are the same, the data structures were imported
                    equals.append(file)
                    if verbose:
                        print("Same parameters\nThe pickle file is directly imported from: ", file)
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
            if verbose:
                print("Parameters are changed\n",
                      "BinReadCounter is running with actual parameters\n",
                      "IT COULD TAKE A WHILE to create and import the pickle file")
            counter._export_pickle()
            return counter._load_pickle()[1]

    def add_ns_trace(self, fig, chrom=None):
        """This method is used in other plotting methods, in order to add the
        representation of N bases in the same plot

        Args:
            fig (obj): is a go.Figure() object, taken from the different plot
                       methods, it represent the figure to which add the Ns
                       representation
            chrom (int): is a parameter of the script; it is different from None
                        in those methods that uses this information to build up
                        the plot

        Returns:
            A trace to add
        """
        count_N_df = self.load_data()[1]  # second data structure
        # if the chromosome is specified, the trace that could be added
        # have to contain only N counts for that specific chromosome
        if chrom != None:
            single_chrom_N = count_N_df[count_N_df["chr"] == "CH.chr" + str(chrom)]
            return fig.add_trace(go.Scatter(x=single_chrom_N["bin"],
                                            y=single_chrom_N["N_count"],
                                            mode="markers",
                                            name="N_counts"))

        else:
            return fig.add_trace(go.Scatter(x=count_N_df["bin"],
                                            y=count_N_df["N_count"],
                                            mode="markers",
                                            name="N_counts"))

    def plot_chrom_sample(self, chrom, sample, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of raw read_counts
        for a specific chromosome and a specific sample

        Args:
            chrom (int): a number representing the chromosome of interest
            sample (str): the name of the sample of interest (it corresponds
                          to the name of the column in the data structure)
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of counts
        """
        read_df = self.load_data()[0]

        fig.update_xaxes(title_text="Chromosomes_Bins")
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
            self.add_ns_trace(fig, chrom=chrom)

        fig.update_layout(
            title="Read Counts - Clone: " + sample + " - Chr: " + str(chrom) + " - Bin Size: " + str(self.bin_size))

        fig.show()

    def plot_chromosome(self, chrom, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of raw_read_counts
        of all samples, but for a specific chromosome of interest

        Args:
            chrom (int): a number representing the chromosome of interest
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of counts
        """
        read_df = self.load_data()[0]

        fig.update_xaxes(title_text="Chromosomes_Bins")
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
            self.add_ns_trace(fig, chrom=chrom)

        fig.update_layout(title="Read Counts - All Clones - Chr: " + str(chrom) + " - Bin Size: " + str(self.bin_size))
        # fig.write_image("plots/counts_chr" + str(chrom) + "_all_" + str(self.bin_size) + ".pdf")

        fig.show()

    def plot_sample(self, sample, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of raw_read_counts
        of all chromosomes, but for a specific sample of interest

        Args:
            sample (str): the name of the sample of interest (it corresponds
                          to the name of the column in the data structure)
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of counts
        """
        read_df = self.load_data()[0]

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        col_read_df = list(read_df.columns)
        for i in range(len(col_read_df[:col_read_df.index(sample)]) + 1):
            if col_read_df[i] == sample:
                fig.add_trace(go.Scatter(x=read_df["bin"],
                                         y=read_df[sample],
                                         mode="markers",
                                         name=str(col_read_df[i])))

        if ns:
            self.add_ns_trace(fig)

        fig.update_layout(title="Read Counts - Clone: " + sample + " - Bin Size: " + str(self.bin_size))

        fig.show()

    def plot_all(self, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of raw_read_counts
        of all chromosomes and all samples

        Args:
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of counts
        """
        read_df = self.load_data()[0]

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        col_list = list(read_df.columns)
        for i in range(len(col_list)):
            if col_list[i] != "bin" and col_list[i] != "chr":
                fig.add_trace(go.Scatter(x=read_df["bin"],
                                         y=read_df[col_list[i]],
                                         mode="markers",
                                         name=str(col_list[i])))

        if ns:
            self.add_ns_trace(fig)

        fig.update_layout(title="Read Counts - All Clones - All Chromosomes - Bin Size: " + str(self.bin_size))
#         # fig.write_image("plots/counts_all_" + str(self.bin_size) + ".pdf")

        fig.show()

    def normalize_bins(self):
        """This method handles the normalization of the raw read_counts
        using an R package, edgeR, imported thanks to rpy2 that provides
        an already implemented function, cpm, for the normalization of a
        table of counts, as well as a series of other function specifically
        implemented for RNA-Seq"""
        read_counts = self.load_data()[0]
        # the following line has to be deleted when the class will work on the entire .bam files
        chr1 = read_counts[read_counts['chr'] == 'CH.chr1']
        # the edgeR package is imported using rpy2 syntax to access to all its built-in functions
        edger = rpackages.importr('edgeR')
        read_counts_edger = {}  # a dictionary of sample: vector_of_counts to work with edger
        col_list = list(chr1.columns)
        for i in range(len(col_list)):
            if col_list[i] != "chr" and col_list[i] != 'bin':
                read_counts_edger[col_list[i]] = robjects.IntVector(chr1[col_list[i]])

        read_counts_edger_df = robjects.DataFrame(read_counts_edger)
        norm_counts = edger.cpm(read_counts_edger_df, normalized_lib_sizes=True)

        return norm_counts  # an edger object (table) of counts

    def plot_norm_data_chr_sample(self, chrom, sample, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of normalized_read_counts
        of a specific chromosome of a specific sample

        Args:
            chrom (int): a number representing the chromosome of interest
            sample (str): the name of the sample of interest (it corresponds
                          to the name of the column in the data structure)
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of normalized counts
        """
        read_counts = self.load_data()[0]
        norm_counts = self.normalize_bins()

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

        single_chrom = read_counts[read_counts["chr"] == "CH.chr" + str(chrom)]
        col_list = list(single_chrom.columns)
        for i in range(1, len(col_list[1:])):
            if col_list[i] == sample:
                fig.add_trace(go.Scatter(x=single_chrom["bin"],
                                         y=list(norm_counts.rx(True, i)),
                                         mode="markers",
                                         name=str(col_list[i + 1])))

        if ns:
            self.add_ns_trace(fig, chrom=chrom)

        fig.update_layout(title="Normalized Read Counts - Clone: " + sample + " - Chr: " + str(chrom) + " - Bin Size: "
                                + str(self.bin_size))

        fig.show()

    def plot_norm_data_chr(self, chrom, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of normalized_read_counts
        of a specific chromosome for all samples

        Args:
            chrom (int): a number representing the chromosome of interest
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of normalized counts
        """
        read_counts = self.load_data()[0]
        norm_counts = self.normalize_bins()

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

        single_chrom = read_counts[read_counts["chr"] == "CH.chr" + str(chrom)]
        col_list = list(single_chrom.columns)
        for i in range(1, len(col_list[1:])):
            fig.add_trace(go.Scatter(x=single_chrom["bin"],
                                     y=list(norm_counts.rx(True, i)),
                                     mode="markers",
                                     name=str(col_list[i + 1])))
        if ns:
            self.add_ns_trace(fig, chrom=chrom)

        fig.update_layout(
            title="Normalized Read Counts - Clone: all - Chr: " + str(chrom) + " - Bin Size: " + str(self.bin_size))

        fig.show()

    def plot_norm_data_sample(self, sample, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of normalized_read_counts
        in all chromosomes of a specific sample

        Args:
            sample (str): the name of the sample of interest (it corresponds
                          to the name of the column in the data structure)
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of normalized counts
        """
        read_counts = self.load_data()[0]
        norm_counts = self.normalize_bins()
        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

        single_sample = read_counts[["chr", "bin", sample]]
        print(single_sample)
        col_list = list(single_sample.columns)
        for i in range(len(col_list)):
            if col_list[i] == sample:
                fig.add_trace(go.Scatter(x=single_sample["bin"],
                                         y=list(norm_counts.rx(True, i)),
                                         mode="markers",
                                         name=str(col_list[i])))

        if ns:
            self.add_ns_trace(fig)

        fig.update_layout(
            title="Normalized Read Counts - Clone: " + sample + " - Chr: all - Bin Size: " + str(self.bin_size))

        fig.show()

    def plot_norm_data_all(self, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of normalized_read_counts
        in all chromosomes and for all samples

        Args:
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of normalized counts
        """
        read_counts = self.load_data()[0]
        norm_counts = self.normalize_bins()
        col_list = list(read_counts.columns)
        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

        for i in range(1, len(col_list[1:])):
            fig.add_trace(go.Scatter(x=read_counts["bin"],
                                     y=list(norm_counts.rx(True, i)),
                                     mode="markers",
                                     name=col_list[i + 1]))

        if ns:
            self.add_ns_trace(fig)

        fig.update_layout(title="Normalized Read Counts - Clone: all - Chr: all - Bin Size: " + str(self.bin_size))

        fig.show()

    def get_significant_diff(self, chr=None, sample=None):
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
                        help="The length of the segments that divide the chromosomes equally")

    parser.add_argument("-f", "--folder",
                        type=str,
                        default="./",
                        help="The path to the folder in which are situated the files to be analyzed (.bam)")

    parser.add_argument("-fl", "--flag_list",
                        nargs="+",
                        default=["0", "16", "99", "147", "163", "83"],
                        help="""A list of the bitwise-flags in SAM format that identify the reads to be counted 
                        during analyses; if different flags wants to be added, add them as strings 
                        (e.g. "177" "129")""")

    parser.add_argument("-op", "--output_pickle", 
                        default='./',
                        help="The folder path to which search the pickle file already created")

    parser.add_argument("-r", "--reference",
                        type=str,
                        default="chardonnay_primary_contigs_chromosome-order.fa.gz",
                        help="""The path of the reference file, if not specified, the coverage is calculated with the 
                        mean of the mapping reads for all the samples""")

    parser.add_argument("-ch", "--chromosome",
                        default=None,
                        type=int,
                        help="""The number of the chromosome of interest for the plot of counts""")
                    
    parser.add_argument("-s", "--sample", 
                        default=None,
                        type=str,
                        help="""[optional] - The name of the clone of interest for the plot of counts""")

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

    analyzer.load_data(verbose=True)

    if not os.path.exists("plots"):
        os.mkdir("plots")

    if args.chromosome and args.sample:
        if args.Ns_count:
            analyzer.plot_chrom_sample(args.chromosome, args.sample, args.Ns_count)
            analyzer.plot_norm_data_chr_sample(args.chromosome, args.sample, args.Ns_count)
        else:
            analyzer.plot_chrom_sample(args.chromosome, args.sample)
            analyzer.plot_norm_data_chr_sample(args.chromosome, args.sample)

    elif args.chromosome:
        if args.Ns_count:
            analyzer.plot_chromosome(args.chromosome, args.Ns_count)
            analyzer.plot_norm_data_chr(args.chromosome, args.Ns_count)
        else:
            analyzer.plot_chromosome(args.chromosome)
            analyzer.plot_norm_data_chr(args.chromosome)

    elif args.sample:
        if args.Ns_count:
            analyzer.plot_sample(args.sample, args.Ns_count)
            analyzer.plot_norm_data_sample(args.sample, args.Ns_count)
        else:
            analyzer.plot_sample(args.sample)
            analyzer.plot_norm_data_sample(args.sample)

    else:
        if args.Ns_count:
            analyzer.plot_all(args.Ns_count)
            analyzer.plot_norm_data_all(args.Ns_count)
        else:
            analyzer.plot_all()
            analyzer.plot_norm_data_all()