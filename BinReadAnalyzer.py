""" Bin read analyzer class"""

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
    def __init__(self, folder_name, bin_size, ref_genome, flags, out_pickle):
        self.folder = folder_name
        self.bin_size = bin_size
        self.flags = flags 
        self.ref = ref_genome
        self.out = out_pickle
        # self.control_file = None

    def load_data(self, verbose=False):
        """If the pickle file is present and the parameters are not
        changed, the data structures were kept otherwise it runs the
        methods load_reads and loads_Ns of the imported class BinReadCounter
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
        count_N_df = self.load_data()[1]
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
        read_df = self.load_data()[0]

        if not os.path.exists("plots"):
            os.mkdir("plots")

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
        read_df = self.load_data()[0]

        if not os.path.exists("plots"):
            os.mkdir("plots")

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
        fig.write_image("plots/counts_chr" + str(chrom) + "_all_" + str(self.bin_size) + ".pdf")

        fig.show()

    def plot_sample(self, sample, ns=False, fig=go.Figure()):
        read_df = self.load_data()[0]

        if not os.path.exists("plots"):
            os.mkdir("plots")

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

    def plot_all(self, ns=False):
        read_df = self.load_data()[0]

        if not os.path.exists("plots"):
            os.mkdir("plots")

        fig = go.Figure()
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
        fig.write_image("plots/counts_all_" + str(self.bin_size) + ".pdf")

        fig.show()

    def normalize_bins(self):

        read_counts = self.load_data()[0]
        # the following line has to be deleted when the class will work on the entire .bam files
        chr1 = read_counts[read_counts['chr'] == 'CH.chr1']
        edger = rpackages.importr('edgeR')
        read_counts_edger = {}
        col_list = list(chr1.columns)
        for i in range(len(col_list)):
            if col_list[i] != "chr" and col_list[i] != 'bin':
                read_counts_edger[col_list[i]] = robjects.IntVector(chr1[col_list[i]])  # include after testing

        read_counts_edger_df = robjects.DataFrame(read_counts_edger)

        # in order to obtain a normalized counts table, is not necessary
        # have a DGEList object, it is sufficient to have a table of raw
        # count as in this case
        # summarized = edger.DGEList(counts=read_counts_edger_df)
        # norm_counts = edger.cpm(summarized, normalized_lib_sizes=False)

        norm_counts = edger.cpm(read_counts_edger_df, normalized_lib_sizes=True)

        # log_cpm = edger.cpm(summarized, log=True) NOT NECESSARY
        # norm_factors = edger.calcNormFactors(summarized) NOT NECESSARY
        # pseudo_counts = edger.estimateCommonDisp(norm_factors) NOT USEFUL NOT NORM_COUNTS
        # print(norm_counts.rx(True, 1)) in other cases useful in selecting an entire column
        return norm_counts

        # d = {'bin': robjects.IntVector(chr1['bin']), 'counts': robjects.IntVector(chr1['Ch15_3_1_ref_IlluminaPE_aligns
        # _Primary_chr1'])}

    def plot_norm_data_chr_sample(self, chrom, sample, ns=False, fig=go.Figure()):
        if not os.path.exists("plots"):
            os.mkdir("plots")
        read_counts = self.load_data()[0]
        norm_counts = self.normalize_bins()

        if not os.path.exists("plots"):
            os.mkdir("plots")

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
        if not os.path.exists("plots"):
            os.mkdir("plots")

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
        if not os.path.exists("plots"):
            os.mkdir("plots")
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
        if not os.path.exists("plots"):
            os.mkdir("plots")

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

    analyzer.load_data(verbose=True)

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