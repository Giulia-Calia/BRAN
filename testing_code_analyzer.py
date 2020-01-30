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
import pandas as pd
import numpy as np
from testing_code_counter import TestingBinReadCounter
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

pandas2ri.activate()
base = importr('base')
utils = importr('utils')
stats = importr('stats')
edger = rpackages.importr('edgeR')


# --------------------------
# install tzlocal for pandas2ri
# --------------------------


class TestingBinReadAnalyzer:
    """This class provides different analysis of the data coming form
    BinReadCounter class; here imported.

    Options:
    -   plot read_counts for all chromosomes, for a specific one or for
        a specific sample
    -   normalize read_counts
    -   plot normalized_read_counts for all chromosomes, for a specific
        one or for a specific sample

    Attributes:
            folder_path: the folder in which the necessary files are stored
            bin_size: the length of the bin in which each chromosome is
                      divided (in bp)
            ref_genome: a string being the name of the reference_genome file
            flags: a list of string, each corresponding to a particular
                   pairwise-flag in the .bam/.sam format
            out_pickle: default is None because it takes the name of an
                        already existing file or the name of a new file
                        created during the running of the code
    """

    def __init__(self, bam_folder_path, bin_size, ref_genome, flags, out_pickle):
        self.folder = bam_folder_path
        self.bin_size = bin_size
        self.flags = flags
        self.ref = ref_genome
        self.out = out_pickle
        self.parameters = None
        self.norm = None
        self.fold_change = None
        self.sig_fc = None

    def set_parameters(self, param):
        self.parameters = param

    def set_norm(self, norm):
        self.norm = norm

    def set_fold_change(self, diff):
        self.fold_change = diff

    def set_sig_fc(self, fc):
        self.sig_fc = fc

    def load_data(self, other_cigar_filters, cigar_filter=None, reference=None, read_info=None, unmapped=None,
                  verbose=False):
        """If the pickle file is already present in the folder and the parameters
        are not changed, the data structures were kept, otherwise it runs the
        methods _export_pickle() and load_data() of the imported class BinReadCounter
        in order to create a new pickle file and to import the corresponding
        data-structures

        Args:
            cigar_filter (bool)
            other_cigar_filters (bool)
            reference (bool): true if the reference is necessary
            read_info (bool): true if read info are necessary
            verbose (bool): default is False, but if set to True, it allows to print
                            and see if a pickle file is already present or if the
                            Analyzer is running to create a new one

        Returns:
            A list of data structures from the pickle file or form the BinReadCounter;
            the first data structure is the read_counts and the second is the N bases
            count
        """
        counter = TestingBinReadCounter(self.folder, self.bin_size, self.flags, self.ref, self.out)
        found = None
        bam_files = []
        i = 0
        list_bam_dir = os.listdir(self.folder)
        # print("bam folder: ", self.folder)
        list_pick_dir = os.listdir(self.out)
        # print("pickle folder: ", self.out)
        for file in list_bam_dir:
            if file.endswith(".bam"):
                bam_files.append(file[:file.find(".bam")])
        # print("bam files:", bam_files)

        # verify if a pickle file is present in the desired folder
        for file in list_pick_dir:
            if file.endswith(".p"):
                # if so, checks if the parameters for the building-up of
                # the previous pickle file, are the same as the actual ones
                # with open(file, "rb") as input_pickle:
                #     print("ok")
                #     counter = pickle.load(input_pickle)
                try:
                    parameters = counter._load_pickle(file)
                    if parameters["bin_size"] == self.bin_size and \
                            parameters["flags"] == self.flags and \
                            all(el in parameters["bam"] for el in bam_files) and \
                            parameters["cigar_filt"] == cigar_filter and \
                            parameters["other_cigar_filt"] == other_cigar_filters and \
                            parameters["unmapped"] == unmapped and \
                            parameters["ref"] == counter.get_ref_name() and \
                            parameters["info"] == read_info:
                        # print(all(file in parameters["bam"] for file in bam_files))
                        # print(parameters["bam"])
                        print(parameters)
                        found = True
                        if verbose:
                            print("Same parameters; import from: ", file)

                        self.set_parameters(parameters)
                        return self.parameters

                except:
                    i += 1
                    continue
        if verbose:
            print("Warning:\nException occurred - other processes running at this moment: ", i)
        if not found:
            # if not found, none of the pickle files in the current directory
            # have the same parameters of the actual running or a file pickle is not p
            # resent at all in the directory thus the algorithm uses the modules of the
            # BinReadCounter to calculate the new data structures and it saves a new pickle file
            if verbose:
                print("Parameters are changed or no pickle file exists",
                      "\nBinReadCounter is running with actual parameters",
                      "\nIT COULD TAKE A WHILE to create and import the pickle file")

            counter._export_pickle(other_cigar_filters, cigar_filter, reference, read_info, unmapped)
            name_pickle = counter.pickle_file_name(other_cigar_filters, cigar_filter, reference, read_info, unmapped)
            parameters = counter._load_pickle(name_pickle)

            self.set_parameters(parameters)
            print(self.parameters)
            return self.parameters

    def normalize_bins(self):  # , df_counts, bin_size
        """This method handles the normalization of the raw read_counts
        using an R package, edgeR, imported thanks to rpy2 that provides
        an already implemented function, cpm, for the normalization of a
        table of counts, as well as a series of other function specifically
        implemented for RNA-Seq"""
        read_counts = self.parameters["read_counts"]
        # read_counts = df_counts
        # the edgeR package is imported using rpy2 syntax to access to all its built-in functions
        read_counts_edger = {}  # a dictionary of sample: vector_of_counts to work with edger
        col_list = list(read_counts.columns)
        group = []
        for i in range(len(col_list)):
            if col_list[i] != "index" and col_list[i] != "chr" and col_list[i] != 'bin':
                # creates an element of the dictionary [sample_name]: r_vector_of_counts
                read_counts_edger[col_list[i]] = robjects.IntVector(read_counts[col_list[i]])
                if "ref_ref" in col_list[i]:
                    group.append(1)
                else:
                    group.append(2)
        read_counts_edger_df = robjects.DataFrame(read_counts_edger)  # R data frame

        # --create a DGEList obj--
        dge_list = edger.DGEList(counts=read_counts_edger_df, group=robjects.IntVector(group))
        # --calc norm factors--
        norm_fact = edger.calcNormFactors(dge_list)
        # --normalization function--
        # norm_counts = edger.cpm(read_counts_edger_df, log=True)
        norm_counts = edger.cpm(norm_fact, normalized_lib_sizes=True)  # the log is calculated with fold-change function

        # norm_counts_dict = {}
        # for i in range(1, norm_counts.ncol + 1):
        #     norm_counts_dict[norm_counts.colnames[i - 1]] = list(norm_counts.rx(True, i))
        # norm_counts_df = pd.DataFrame(norm_counts_dict)
        # print(norm_counts_df)

        self.set_norm(norm_counts)
        return self.norm  # an edger object (table) of counts

    def get_fold_change(self, saving_folder, control_name=None):  # df_counts, bin_size,
        """This method calculates the pairwise logFoldChange for the different clones
        against the control plant counts, in order to see, if are present, significant
        differences in terms of mapping reads"""
        if control_name is not None:
            read_counts = self.parameters["read_counts"]
            # read_counts = df_counts

            control = control_name
            read_counts_edger = {}  # {}  # a dictionary of sample: vector_of_counts to work with edger
            col_list = list(read_counts.columns)
            fold_change = {}

            bcv = 0.1  # reasonable dispersion value; typical values for common BCV (square-root dispersion) of data-set
            # arising from well-controlled experiments are 0.4 for human data, 0.1 for data of genetically identical
            # model organisms or 0.01 fo technical replicates

            for i in range(len(col_list)):
                if col_list[i] != "index" and col_list[i] != "chr" and col_list[i] != 'bin':
                    read_counts_edger[col_list[i]] = robjects.r.matrix(robjects.IntVector(read_counts[col_list[i]]))

            # --with normalized data--
            for i in range(1, self.norm.ncol + 1):
                col = self.norm.rx(True, i)
                col_name = self.norm.colnames[i - 1]
                if col_name != control:
                    print("\n", col_name)
                    tmp_dict = {col_name: col, control: self.norm.rx(True, control)}
                    tmp_df = robjects.DataFrame(tmp_dict)
                    dge_list = edger.DGEList(robjects.DataFrame(tmp_df), group=robjects.IntVector([2, 1]))
                    # print(dge_list)
                    dge_list_ed = edger.estimateDisp(dge_list)
                    exact_test = edger.exactTest(dge_list_ed, dispersion=bcv ** 2)
                    # print(edger.topTags(exact_test))
                    exact_test_df = robjects.DataFrame(exact_test[0])
                    # print(type(exact_test_df))
                    exact_test_pd_df = pandas2ri.ri2py_dataframe(exact_test_df)
                    # print(exact_test_pd_df)
                    fold_change[col_name] = exact_test_pd_df
                    fold_change_df = fold_change[col_name]
                    fold_change_df.to_csv(saving_folder +
                                          "fold_change.txt",
                                          sep="\t")

            self.set_fold_change(fold_change)
            return self.fold_change

        else:
            print("No control group is specified, please try again, specifying the column name of control file with "
                  "parameter '-co'")

            # --with raw data--
            # for el in read_counts_edger:
            #     if el != control:
            #         # print(el)
            #         tmp_dict = {el: read_counts_edger[el], control: read_counts_edger[control]}
            #         # print(tmp_dict)
            #         tmp_df = robjects.DataFrame(tmp_dict)
            #         dge_list = edger.DGEList(robjects.DataFrame(tmp_df), group=robjects.IntVector([2, 1]))
            #         # print(dge_list)
            #         dge_list_ed = edger.estimateDisp(dge_list)
            #         # print(dge_list_ed)
            #         exact_test = edger.exactTest(dge_list_ed, dispersion=bcv ** 2)
            #         # print(exact_test)
            #         # print(edger.topTags(exact_test))
            #         exact_test_df = robjects.DataFrame(exact_test[0])
            #         # print(type(exact_test_df))
            #         exact_test_pd_df = pandas2ri.ri2py_dataframe(exact_test_df)
            #         # print(exact_test_pd_df)
            #         fold_change[el] = exact_test_pd_df

    def get_sig_pos(self, fc, p_value):  # , df_counts, bin_size
        """This method return the chromosome position in terms of basis, of bins
        retained significantly different from the control"""
        fold_change = self.fold_change
        read_counts = self.parameters["read_counts"]
        # read_counts = df_counts
        sig_data = {"clone": [], "bin": [], "chr": [], "genome_position": [], "index": [], "logFC": [], "PValue": []}
        for el in fold_change:
            sig_fc_values = fold_change[el][fold_change[el]["logFC"] >= fc]
            sig_fc_neg_values = fold_change[el][fold_change[el]["logFC"] <= -fc]
            sig_fc = pd.concat([sig_fc_values, sig_fc_neg_values])
            sig_fc_pvalues = sig_fc[sig_fc["PValue"] <= p_value]

            for i in sig_fc_pvalues.index.values:
                effective_bin = read_counts["bin"].iloc[i]
                sig_data["clone"].append(el)
                sig_data["bin"].append(effective_bin)  # not index but bin!!!
                sig_data["chr"].append(read_counts["chr"].iloc[i])
                sig_data["genome_position"].append(effective_bin * self.bin_size)
                sig_data["index"].append(i)
                sig_data["logFC"].append(sig_fc_pvalues["logFC"][i])
                sig_data["PValue"].append(sig_fc_pvalues["PValue"][i])

        sig_data_df = pd.DataFrame(sig_data)
        self.set_sig_fc(sig_data_df)
        return self.sig_fc

    def get_no_sig_pos(self, fc, pvalue):
        fold_change = self.fold_change
        read_counts = self.parameters["read_counts"]
        no_sig_data = {"clone": [], "bin": [], "chr": [], "genome_position": [], "logFC": [], "PValue": []}

        # for el in fold_change:
        pass

    def plot_background(self, fig):  # , df_counts
        read_counts = self.parameters["read_counts"]
        # read_counts = df_counts
        coordinates_x = []
        coordinates_y = []
        chromosomes = []
        start = 0

        for ch in read_counts["chr"]:
            if ch[ch.find("c"):] not in chromosomes:
                chromosomes.append(ch[ch.find("c"):])
        print(chromosomes)
        for chrom in chromosomes:
            single_df = read_counts[read_counts["chr"] == "CH." + chrom]
            # print(single_df)
            # here start and end are updated every time, to allow concatenation of chromosomes on the plot
            # the average position within the interval of each chromosome take in account only the length of the
            # interval itself
            length = (single_df["bin"].iloc[-1] + 1) * self.parameters["bin_size"]
            tmp_end = start + length
            # print("start", start)
            # print("end", tmp_end)
            avg = start + length / 2
            coordinates_x.append(avg)
            coordinates_y.append(-2.6)
            if int(chrom[chrom.find("r") + 1:]) % 2 == 0:
                fig.add_shape(go.layout.Shape(type="rect",
                                              xref="x",
                                              yref="paper",
                                              x0=start,
                                              y0=0,
                                              x1=tmp_end,
                                              y1=1,
                                              fillcolor="rgb(230, 230, 250)",  # lavande
                                              opacity=0.5,
                                              layer="below",
                                              line_width=0))
            else:
                fig.add_shape(go.layout.Shape(type="rect",
                                              xref="x",
                                              yref="paper",
                                              x0=start,
                                              y0=0,
                                              x1=tmp_end,
                                              y1=1,
                                              fillcolor="rgb(240, 248, 255)",  # ~light_mint_green
                                              opacity=0.5,
                                              layer="below",
                                              line_width=0))

            start = tmp_end

        fig.add_trace(go.Scatter(x=coordinates_x,
                                 y=coordinates_y,
                                 text=chromosomes,
                                 mode="text",
                                 showlegend=False,
                                 ))

        return fig

    def add_threshold_fc(self, fig, fc):
        read_counts = self.parameters["read_counts"]

        fig.add_shape(go.layout.Shape(type="line",
                                      x0=0,
                                      y0=fc,
                                      x1=len(read_counts) * self.bin_size,
                                      y1=fc))

        fig.add_shape(go.layout.Shape(type="line",
                                      x0=0,
                                      y0=-fc,
                                      x1=len(read_counts) * self.bin_size,
                                      y1=-fc))

        fig.update_shapes(dict(xref="x",
                               yref="y",
                               line=dict(color="crimson",  # dark red
                                         width=1)))

        return fig

    def add_ns_trace(self, fig, reference=None, chrom=None):
        """This method is used in other plotting methods, in order to add the
        representation of N bases in the same plot

        Args:
            reference (bool): true if the reference is declared
            fig (obj): is a go.Figure() object, taken from the different plot
                       methods, it represent the figure to which add the Ns
                       representation
            chrom (int): is a parameter of the script; it is different from None
                        in those methods that uses this information to build up
                        the plot

        Returns:
            A trace to add
        """
        if reference is not None:
            count_N_df = self.parameters["n_counts"]  # second data structure
            # if the chromosome is specified, the trace that could be added
            # have to contain only N counts for that specific chromosome
            if chrom is not None:
                single_chrom_N = count_N_df[count_N_df["chr"] == "CH.chr" + str(chrom)]
                return fig.add_trace(go.Scatter(x=single_chrom_N["index"],
                                                y=single_chrom_N["N_count"],
                                                mode="markers",
                                                name="N_counts"))

            else:
                return fig.add_trace(go.Scatter(x=count_N_df["index"],
                                                y=count_N_df["N_count"],
                                                mode="markers",
                                                name="N_counts"))
        else:
            print("Sorry, to add this trace to the plot, you have to specify the reference file name\nPlease try again")

    def plot_counts_distributions(self, saving_folder, cigar):
        read_counts = self.parameters["read_counts"]
        counts_to_plot = []
        col_labels = []

        # --plot_raw_counts_distribution--
        fig_all = go.Figure()
        fig_all.update_xaxes(title_text="Raw_Counts")
        fig_all.update_yaxes(title_text="Frequency")

        for col in read_counts.columns:
            if col != "index" and col != "chr" and col != 'bin':
                counts_to_plot.append(read_counts[col])
                col_labels.append(col)
                fig_all.add_trace(go.Histogram(x=read_counts[col],
                                               name=col))

        # --plot_norm_counts_distribution--

        fig_all_norm = go.Figure()
        fig_all_norm.update_traces(opacity=0.75)
        fig_all_norm.update_xaxes(title_text="Normalized_Counts")
        fig_all_norm.update_yaxes(title_text="Frequency")

        for i in range(0, self.norm.ncol):
            fig_all_norm.add_trace(go.Histogram(x=list(self.norm.rx(True, i + 1)),
                                                name=str(self.norm.colnames[i])))

        # --plot_log_scaled_norm_counts_distribution--
        fig_log_norm = go.Figure()
        fig_log_norm.update_xaxes(title_text="Log-scaled_Normalized_Counts")
        fig_log_norm.update_yaxes(title_text="Frequency")

        for key in self.fold_change:
            fig_log_norm.add_trace(go.Histogram(x=self.fold_change[key]["logCPM"],
                                                name=key))

        # --Titles--
        if cigar:

            fig_all.update_layout(title_text="Cigar Filter - Histogram RAW data - BRAN" +
                                             str(self.bin_size) + " - all samples",
                                  legend_orientation="h",
                                  barmode="overlay")
            fig_all_norm.update_layout(title_text="Cigar Filter - Histogram NORMALIZED data - BRAN" +
                                                  str(self.bin_size) + " - all samples",
                                       legend_orientation="h",
                                       barmode="overlay")

            fig_log_norm.update_layout(title_text="Cigar Filter - Histogram LOG-SCALE NORMALIZED data - BRAN" +
                                                  str(self.bin_size) + " - all samples",
                                       legend_orientation="h",
                                       barmode="overlay")

        else:
            fig_all.update_layout(title_text="Histogram RAW data - BRAN" + str(self.bin_size) + " - all samples",
                                  legend_orientation="h",
                                  barmode="overlay")

            fig_all_norm.update_layout(title_text="Histogram NORMALIZED data - BRAN" +
                                                  str(self.bin_size) + " - all samples",
                                       legend_orientation="h",
                                       barmode="overlay")

            fig_log_norm.update_layout(title_text="Histogram LOG-SCALE NORMALIZED data - BRAN" +
                                                  str(self.bin_size) + " - all samples",
                                       legend_orientation="h",
                                       barmode="overlay")

        fig_all.show()
        fig_all_norm.show()
        fig_log_norm.show()
        save_fig_all = fig_all.write_image(saving_folder + "all_sample_dist_" + str(self.bin_size) + ".pdf",
                                           width=1280,
                                           height=1024)
        save_fig_all_norm = fig_all_norm.write_image(saving_folder + "all_sample_norm_dist_" +
                                                     str(self.bin_size) + ".pdf",
                                                     width=1280,
                                                     height=1024)
        save_fig_log_norm = fig_log_norm.write_image(saving_folder + "all_sample_norm_log_dist_" +
                                                     str(self.bin_size) + ".pdf",
                                                     width=1280,
                                                     height=1024)

    def plot_chrom_sample(self, saving_folder, reference, chrom, sample, cigar, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of raw read_counts
        for a specific chromosome and a specific sample

        Args:
            reference (bool): true if the reference is declared
            chrom (int): a number representing the chromosome of interest
            sample (str): the name of the sample of interest (it corresponds
                          to the name of the column in the data structure)
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of counts
        """
        read_counts = self.parameters["read_counts"]

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        c_name = None
        for c in list(read_counts["chr"]):
            if c.endswith(str(chrom)):
                c_name = c

        single_chrom = read_counts[read_counts["chr"] == c_name]
        col_list = list(single_chrom.columns)
        for i in range(len(col_list[:col_list.index(sample)]) + 1):
            if col_list[i] == sample:
                fig.add_trace(go.Scatter(x=single_chrom["index"],
                                         y=single_chrom[sample],
                                         mode="markers",
                                         name=str(col_list[i])))

        if ns:
            self.add_ns_trace(fig, reference=reference, chrom=chrom)

        if cigar:
            fig.update_layout(title="Cigar Filter - Read Counts - Clone: " + sample +
                                    " - Chr: " + str(chrom) + " - Bin Size: " + str(self.bin_size),
                              legend_orientation="h")
        else:
            fig.update_layout(title="Read Counts - Clone: " + sample +
                                    " - Chr: " + str(chrom) + " - Bin Size: " + str(self.bin_size),
                              legend_orientation="h")

        fig.show()
        save_fig = fig.write_image(saving_folder + "scatter_counts_chr" +
                                   str(chrom) + sample + str(self.bin_size) + ".pdf",
                                   width=1280,
                                   height=1024)

    def plot_chromosome(self, saving_folder, reference, chrom, cigar, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of raw_read_counts
        of all samples, but for a specific chromosome of interest

        Args:
            reference (bool): true if the reference is declared
            chrom (int): a number representing the chromosome of interest
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of counts
        """
        read_counts = self.parameters["read_counts"]

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        c_name = None
        for c in list(read_counts["chr"]):
            if c.endswith(str(chrom)):
                c_name = c

        single_chrom = read_counts[read_counts["chr"] == c_name]
        col_list = list(single_chrom.columns)
        for i in range(len(col_list)):
            if col_list[i] != "index" and col_list[i] != "chr" and col_list[i] != "bin":
                fig.add_trace(go.Scatter(x=single_chrom["index"],
                                         y=single_chrom[col_list[i]],
                                         mode="markers",
                                         name=str(col_list[i])))
        if ns:
            self.add_ns_trace(fig, reference=reference, chrom=chrom)

        if cigar:
            fig.update_layout(title="Cigar Filter - Read Counts - All Clones - Chr: " + str(chrom) +
                                    " - Bin Size: " + str(self.bin_size),
                              legend_orientation="h")

        else:
            fig.update_layout(title="Read Counts - All Clones - Chr: " + str(chrom) +
                                    " - Bin Size: " + str(self.bin_size),
                              legend_orientation="h")

        fig.show()
        save_fig = fig.write_image(saving_folder + "scatter_counts_chr" + str(chrom) + str(self.bin_size) + ".pdf",
                                   width=1280,
                                   height=1024)

    def plot_sample(self, saving_folder, reference, sample, cigar, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of raw_read_counts
        of all chromosomes, but for a specific sample of interest

        Args:
            reference (bool): true if the reference is declared
            sample (str): the name of the sample of interest (it corresponds
                          to the name of the column in the data structure)
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of counts
        """
        read_counts = self.parameters["read_counts"]

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        col_read_counts = list(read_counts.columns)
        for i in range(len(col_read_counts[:col_read_counts.index(sample)]) + 1):
            if col_read_counts[i] == sample:
                fig.add_trace(go.Scatter(x=read_counts["index"],
                                         y=read_counts[sample],
                                         mode="markers",
                                         name=str(col_read_counts[i])))

        if ns:
            self.add_ns_trace(fig, reference=reference)

        self.plot_background(fig)

        if cigar:
            fig.update_layout(title="Cigar Filter - Read Counts - Clone: " + sample +
                                    " - Bin Size: " + str(self.bin_size),
                              legend_orientation="h")

        else:
            fig.update_layout(title="Read Counts - Clone: " + sample + " - Bin Size: " + str(self.bin_size),
                              legend_orientation="h")

        fig.show()
        save_fig = fig.write_image(saving_folder + "scatter_counts_" + sample + str(self.bin_size) + ".pdf",
                                   width=1280,
                                   height=1024)

    def plot_all(self, saving_folder, reference, cigar, ns=False, fig=go.Figure()):  # df_counts, bin_size
        """This method allows to obtain a scatter-plot of raw_read_counts
        of all chromosomes and all samples

        Args:
            reference (bool): true if the reference is declared
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of counts
        """
        read_counts = self.parameters["read_counts"]
        # read_counts = df_counts
        col_list = list(read_counts.columns)

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        for i in range(len(col_list)):
            if col_list[i] != "index" and col_list[i] != "chr" and col_list[i] != "bin":
                fig.add_trace(go.Scatter(x=read_counts["index"],
                                         y=read_counts[col_list[i]],
                                         mode="markers",
                                         name=str(col_list[i])))

        if ns:
            self.add_ns_trace(fig, reference=reference)

        self.plot_background(fig)

        if cigar:
            fig.update_layout(title="Cigar Filter - Read Counts - All Clones - All Chromosomes - Bin Size: " +
                                    str(self.bin_size),
                              legend_orientation="h")

        else:
            fig.update_layout(title="Read Counts - All Clones - All Chromosomes - Bin Size: " + str(self.bin_size),
                              legend_orientation="h")

        fig.show()
        save_fig = fig.write_image(saving_folder + "scatter_all_counts" + str(self.bin_size) + ".pdf",
                                   width=1280,
                                   height=1024)

    def plot_norm_data_chr_sample(self, saving_folder, reference, chrom, sample, cigar, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of normalized_read_counts
        of a specific chromosome of a specific sample

        Args:
            reference (bool): true if the reference is declared
            chrom (int): a number representing the chromosome of interest
            sample (str): the name of the sample of interest (it corresponds
                          to the name of the column in the data structure)
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of normalized counts
        """
        read_counts = self.parameters["read_counts"]
        norm_counts = self.norm

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

        c_name = None
        for c in list(read_counts["chr"]):
            if c.endswith(str(chrom)):
                c_name = c

        single_chrom = read_counts[read_counts["chr"] == c_name]
        for i in range(0, norm_counts.ncol):
            if norm_counts.colnames[i] == sample:
                fig.add_trace(go.Scatter(x=single_chrom["index"],
                                         y=list(norm_counts.rx(True, i + 1)[single_chrom["index"].iloc[0]:
                                                                            single_chrom["index"].iloc[-1] + 1]),
                                         mode="markers",
                                         name=str(norm_counts.colnames[i])))
        if ns:
            self.add_ns_trace(fig, reference=reference, chrom=chrom)

        if cigar:
            fig.update_layout(title="Cigar Filter - Normalized Read Counts - Clone: " +
                                    sample +
                                    " - Chr: " +
                                    str(chrom) +
                                    " - Bin Size: " +
                                    str(self.bin_size),
                              legend_orientation="h")

        else:
            fig.update_layout(title="Normalized Read Counts - Clone: " +
                                    sample +
                                    " - Chr: " +
                                    str(chrom) +
                                    " - Bin Size: " +
                                    str(self.bin_size),
                              legend_orientation="h")

        fig.show()
        save_fig = fig.write_image(saving_folder + "scatter__norm_counts_chr" +
                                   str(chrom) + sample + str(self.bin_size) + ".pdf",
                                   width=1280,
                                   height=1024)

    def plot_norm_data_chr(self, saving_folder, reference, chrom, cigar, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of normalized_read_counts
        of a specific chromosome for all samples

        Args:
            reference (bool): true if the reference is declared
            chrom (int): a number representing the chromosome of interest
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of normalized counts
        """
        read_counts = self.parameters["read_counts"]
        norm_counts = self.norm

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

        c_name = None
        for c in list(read_counts["chr"]):
            if c.endswith(str(chrom)):
                c_name = c

        single_chrom = read_counts[read_counts["chr"] == c_name]
        for i in range(0, norm_counts.ncol):
            fig.add_trace(go.Scatter(x=single_chrom["index"],
                                     y=list(norm_counts.rx(True, i + 1)[single_chrom["index"].iloc[0]:
                                                                        single_chrom["index"].iloc[-1] + 1]),
                                     mode="markers",
                                     name=norm_counts.colnames[i]))
        if ns:
            self.add_ns_trace(fig, reference=reference, chrom=chrom)

        if cigar:
            fig.update_layout(title="Cigar Filter - Normalized Read Counts - Clone: all - Chr: " +
                                    str(chrom) +
                                    " - Bin Size: " +
                                    str(self.bin_size),
                              legend_orientation="h")

        else:
            fig.update_layout(title="Normalized Read Counts - Clone: all - Chr: " +
                                    str(chrom) +
                                    " - Bin Size: " +
                                    str(self.bin_size),
                              legend_orientation="h")

        fig.show()
        save_fig = fig.write_image(saving_folder + "scatter_norm_counts_chr" + str(chrom) + str(self.bin_size) + ".pdf",
                                   width=1280,
                                   height=1024)

    def plot_norm_data_sample(self, saving_folder, reference, sample, cigar, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of normalized_read_counts
        in all chromosomes of a specific sample

        Args:
            reference (bool): true if the reference is declared
            sample (str): the name of the sample of interest (it corresponds
                          to the name of the column in the data structure)
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of normalized counts
        """
        read_counts = self.parameters["read_counts"]
        norm_counts = self.norm

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

        single_sample = read_counts[["index", "chr", "bin", sample]]
        col_list = list(single_sample.columns)
        for i in range(0, norm_counts.ncol):
            if norm_counts.colnames[i] == sample:
                fig.add_trace(go.Scatter(x=single_sample["index"],
                                         y=list(norm_counts.rx(True, i + 1)),
                                         mode="markers",
                                         name=str(norm_counts.colnames[i])))

        if ns:
            self.add_ns_trace(fig, reference=reference)

        self.plot_background(fig)

        if cigar:
            fig.update_layout(title="Cigar Filter - Normalized Read Counts - Clone: " +
                                    sample +
                                    " - Chr: all - Bin Size: " +
                                    str(self.bin_size),
                              legend_orientation="h")

        else:
            fig.update_layout(title="Normalized Read Counts - Clone: " +
                                    sample +
                                    " - Chr: all - Bin Size: " +
                                    str(self.bin_size),
                              legend_orientation="h")

        fig.show()
        save_fig = fig.write_image(saving_folder + "scatter_norm_counts_" + sample + str(self.bin_size) + ".pdf",
                                   width=1280,
                                   height=1024)

    def plot_norm_data_all(self, saving_folder, reference, cigar, ns=False, fig=go.Figure()):  # df_counts, bin_size
        """This method allows to obtain a scatter-plot of normalized_read_counts
        in all chromosomes and for all samples

        Args:
            reference (bool): true if the reference is declared
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of normalized counts
        """
        read_counts = self.parameters["read_counts"]
        # read_counts = df_counts
        norm_counts = self.norm

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

        for i in range(0, norm_counts.ncol):
            fig.add_trace(go.Scatter(x=read_counts["index"],
                                     y=list(norm_counts.rx(True, i + 1)),
                                     mode="markers",
                                     name=norm_counts.colnames[i]))

        if ns:
            self.add_ns_trace(fig, reference=reference)

        self.plot_background(fig)

        if cigar:
            fig.update_layout(title="Cigar Filter - Normalized Read Counts - Clone: all - Chr: all - Bin Size: " +
                                    str(self.bin_size),
                              legend_orientation="h")

        else:
            fig.update_layout(title="Normalized Read Counts - Clone: all - Chr: all - Bin Size: " + str(self.bin_size),
                              legend_orientation="h")

        fig.show()
        save_fig = fig.write_image(saving_folder + "scatter_norm_all_counts" + str(self.bin_size) + ".pdf",
                                   width=1280,
                                   height=1024)

    def plot_sig_data_chr_sample(self, saving_folder, fc, p_val, cigar, chrom, sample):
        read_counts = self.parameters["read_counts"]
        # read_counts = df_counts
        c_name = None
        fold_change = self.fold_change
        sig_fc = self.sig_fc
        fig = go.Figure()

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Fold-Change")

        for c in list(read_counts["chr"]):
            if c.endswith("chr" + chrom):
                c_name = c

        # print(sample)
        # print(sig_fc.clone.unique())
        if sample in sig_fc.clone.unique():
            sub_df = sig_fc[sig_fc["clone"] == sample]  # significant data
            # print(sub_df)
            sub_df_chr = sub_df[sub_df["chr"] == c_name]
            # print(sub_df_chr)
            no_sig_fc = fold_change[sample].drop(list(sub_df["bin"]))  # not significant data
            print(no_sig_fc)
            no_sig_chr = no_sig_fc[no_sig_fc["chr"] == c_name]

            # fig.add_trace(go.Scatter(x=sub_df_chr["bin"],
            #                          y=sub_df_chr["logFC"],
            #                          mode="markers",
            #                          name=sample))
            #
            # fig.add_trace(go.Scatter(x=no_sig_chr.index.values,
            #                          y=no_sig_chr["logFC"],
            #                          mode="markers",
            #                          marker=dict(color="rgb(176, 196, 222)"),  # silver
            #                          showlegend=False))

        # self.add_threshold_fc(fig, fc)
        #
        # self.plot_background(fig)

        # if cigar:
        #     fig.update_layout(title="Cigar Filter - Significant Fold Change Counts - Clone: all - Chr: all - "
        #                             "Bin Size: " + str(self.bin_size) + " fc_threshold: " + str(fc) + " p_value: " +
        #                             str(p_val),
        #                       legend_orientation="h")
        # else:
        #     fig.update_layout(title="Significant Fold Change Counts - Clone: all - Chr: all - "
        #                             "Bin Size: " + str(self.bin_size) + " fc_threshold: " + str(fc) + " p_value: " +
        #                             str(p_val),
        #                       legend_orientation="h")
        #
        # fig.show()
        # save_fig = fig.write_image(saving_folder + "counts_fold_change" + str(self.bin_size) + ".pdf",
        #                            width=1280,
        #                            height=1024)

    def plot_sig_data_chr(self):
        pass

    def plot_sig_data_sample(self):
        pass

    def plot_sig_data(self, saving_folder, fc, p_value, cigar):  # , df_counts, bin_size
        read_counts = self.parameters["read_counts"]
        # read_counts = df_counts

        fold_change = self.fold_change
        # print(fold_change['test1_alignSort'].index.values * self.bin_size)
        sig_fc = self.sig_fc
        fig = go.Figure()

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Fold-Change")

        for clone in fold_change.keys():
            if sig_fc.empty:
                fig.add_trace(go.Scatter(x=fold_change[clone].index.values * self.bin_size,
                                         y=fold_change[clone]["logFC"],
                                         mode="markers",
                                         marker=dict(color="rgb(176, 196, 222)"),  # silver
                                         name=clone))

            else:
                sub_df = sig_fc[sig_fc["clone"] == clone]
                no_sig_fc = fold_change[clone].drop(list(sub_df["index"]), axis=0)

                fig.add_trace(go.Scatter(x=sub_df["index"] * self.bin_size,
                                         y=sub_df["logFC"],
                                         mode="markers",
                                         name=clone))

                fig.add_trace(go.Scatter(x=no_sig_fc.index.values * self.bin_size,
                                         y=no_sig_fc["logFC"],
                                         mode="markers",
                                         marker=dict(color="rgb(176, 196, 222)"),  # silver
                                         showlegend=False))

        # fig.add_shape(go.layout.Shape(type="line",
        #                               x0=0,
        #                               y0=fc,
        #                               x1=len(read_counts),
        #                               y1=fc))
        #
        # fig.add_shape(go.layout.Shape(type="line",
        #                               x0=0,
        #                               y0=-fc,
        #                               x1=len(read_counts),
        #                               y1=-fc))
        #
        # fig.update_shapes(dict(xref="x",
        #                        yref="y",
        #                        line=dict(color="crimson",
        #                                  width=1)))

        self.add_threshold_fc(fig, fc)
        self.plot_background(fig)

        if cigar:
            fig.update_layout(title="Cigar Filter - Significant Fold Change Counts - Clone: all - Chr: all - "
                                    "Bin Size: " + str(self.bin_size) + " - Threshold_FC: " + str(fc) +
                                    " - p-value: " + str(p_value),
                              legend_orientation="h")
        else:
            fig.update_layout(title="Significant Fold Change Counts - Clone: all - Chr: all - "
                                    "Bin Size: " + str(self.bin_size) + " - Threshold_FC: " + str(fc) +
                                    " - p-value: " + str(p_value),
                              legend_orientation="h")

        fig.show()
        save_fig = fig.write_image(saving_folder + "counts_fold_change" + str(self.bin_size) + ".pdf",
                                   width=1280,
                                   height=1024)

    def raw_plot_sh_clipping(self, control_ref):
        """"""
        sh_clipping = self.parameters["summary_clip_file"]
        read_counts = self.parameters["read_counts"]
        # print(read_counts["chr"].value_counts())
        clone_df = sh_clipping[sh_clipping["clone"] == "test1_alignSort"]
        chr_col = []

        for el in clone_df["read_id"]:
            for chrom in read_counts["chr"].value_counts().index:
                if str(chrom) in str(el):
                    chr_col.append(chrom)

        # chr_col = pd.Series(chr_col)
        # chr_col.to_csv("chr_col.txt", sep="\t")

        clone_df["chr"] = chr_col
        chr1 = clone_df[clone_df["chr"] == "CH.chr1"]
        chr1_pos = list(chr1["chrom_pos"])
        # print(chr1_pos)
        chr2 = clone_df[clone_df["chr"] == "CH.chr2"]
        chr2_pos = list(chr2["chrom_pos"])
        # print(chr2_pos)
        chr3 = clone_df[clone_df["chr"] == "CH.chr3"]
        chr3_pos = list(chr3["chrom_pos"])
        # print(chr3_pos)

        augmented_pos2 = []
        augmented_pos3 = []
        for pos in chr2_pos:
            augmented_pos2.append(int(pos) + 21500000)

        for pos in chr3_pos:
            augmented_pos3.append(int(pos) + 39000000)

        plot_pos = chr1_pos + augmented_pos2 + augmented_pos3
        clone_df["plot_positions"] = plot_pos

        clone_df.to_csv("chr1_augmented_counts_chr2_3.txt", sep="\t")
        # print(clone_df)

        occur_count = clone_df["plot_positions"].value_counts()  # counts how many times that position is occurred
        # for different reads
        index_count = list(occur_count.index)
        count_df = pd.DataFrame({"plot_positions": index_count, "counts": occur_count})
        clone_df_counts = pd.merge(clone_df, count_df, on="plot_positions")
        print(clone_df_counts)
        # occur_count.to_csv("counts_repeat_pos.txt", sep="\t")
        clone_df_counts.to_csv("correct_counts_pos.txt", sep="\t")

        # print(occur_count)
        pos = clone_df_counts["plot_positions"]
        hover_pos = clone_df_counts["chrom_pos"]
        fig = go.Figure()
        # Here I changed the pop-up information on the plot, writing the exact position of the read on the chromosome,
        # instead the position in the entire genome = x-axis
        fig.add_trace(go.Scatter(x=pos,
                                 y=clone_df_counts["counts"],
                                 hovertext=hover_pos,
                                 hovertemplate=
                                 "<b>Chrom_position</b>: %{hovertext:,}" + "<br><b>Genome_position</b>: %{x:,}",
                                 hoverinfo="text",
                                 mode="markers",
                                 name="test_1"
                                 ))

        self.plot_background(fig)
        fig.update_xaxes(title_text="Genome Length")
        fig.update_yaxes(title_text="Count Repeated Clip Positions")
        fig.update_layout(title="Soft_Hard Clipped Read Positions - "
                                "Bin Size: " + str(self.bin_size))
        fig.show()

    def plot_sh_clipping(self, control_ref, saving_folder):
        """"""
        sh_clipping = self.parameters["summary_clip_file"]
        read_counts = self.parameters["read_counts"]
        chromosomes = sh_clipping["chr"].value_counts().index
        genome_pos = []
        start = 0

        clone_clipping = sh_clipping[sh_clipping["clone"] != control_ref]
        for chrom in chromosomes:
            bin_number = read_counts["chr"].value_counts()[chrom]
            single_df = sh_clipping[sh_clipping["chr"] == chrom]
            single_df_clone = single_df[single_df["clone"] != control_ref]
            aug_length = bin_number * self.parameters["bin_size"]
            # aug_length = self.parameters["chrom_length"][chrom]
            print(aug_length)
            tmp_end = start + aug_length
            print("tmp_end", tmp_end)
            if chrom == chromosomes[0]:
                genome_pos = list(single_df_clone["chrom_pos"])
            else:
                genome_pos += list(single_df_clone["chrom_pos"] + start)
            start = tmp_end

        # print(genome_pos)

        clone_clipping["genome_pos"] = genome_pos
        counts = pd.DataFrame({"genome_pos": clone_clipping["genome_pos"].value_counts().index,
                               "read_rep": clone_clipping["genome_pos"].value_counts()})

        clone_clipping_counts = pd.merge(clone_clipping, counts, on="genome_pos")
        clone_clipping_counts.to_csv(saving_folder +
                                     "clipped_reads_count_" +
                                     clone_clipping["clone"].iloc[0] +
                                     ".txt",
                                     sep="\t")
        # print(clone_clipping_counts)

        pos = clone_clipping["genome_pos"]
        rep = clone_clipping_counts["read_rep"]
        hover_pos = clone_clipping["chrom_pos"]
        fig = go.Figure()
        # Here I changed the pop-up information on the plot, writing the exact position of the read on the chromosome,
        # instead the position in the entire genome = x-axis
        fig.add_trace(go.Scatter(x=pos,
                                 y=rep,
                                 hovertext=hover_pos,
                                 # for CONTROL
                                 # hovertemplate=
                                 # "<b>Chrom_position</b>: %{hovertext:,}" + "<br><b>Genome_position</b>: %{x:,}",
                                 hovertemplate=
                                 "<b>Chrom_position</b>: %{hovertext:,}" + "<br>%{y}",
                                 hoverinfo="text",
                                 mode="markers",
                                 name="test_1"
                                 ))

        self.plot_background(fig)
        fig.update_xaxes(title_text="Genome Length")
        fig.update_yaxes(title_text="Count Repeated Clip Positions")
        fig.update_layout(title=clone_clipping["clone"].iloc[0] + " - Soft_Hard Clipped Read Positions")
        fig.show()
        save_fig = fig.write_image(saving_folder + "clipped_reads_counts" + str(self.bin_size) + ".pdf",
                                   width=1280,
                                   height=1024)


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
                        help="The path to the folder in which are situated the files to be analyzed (.bam)")

    parser.add_argument("-fl", "--flag_list",
                        nargs="+",
                        default=["99", "147", "163", "83"],
                        help="""A list of the bitwise-flags in SAM format that identify the reads to be counted 
                        during analyses; if different flags wants to be added, add them as strings 
                        (e.g. "177" "129")""")

    parser.add_argument("-op", "--output_pickle",
                        default='./',
                        help="The folder path to which search the pickle file already created")

    parser.add_argument("-r", "--reference",
                        type=str,
                        default=None,
                        help="""The path of the reference file, if not specified, the coverage is calculated with the 
                        mean of the mapping reads for all the samples""")

    parser.add_argument("-i", "--read_info",
                        action="store_true",
                        default=None,
                        help="""If specified, a data-frame with information on the read ID, and if the read 
                            and its mate map in the same bin in the same chromosome, is created""")

    parser.add_argument("-ch", "--chromosome",
                        default=None,
                        type=str,
                        help="""The number of the chromosome of interest for the plot of counts""")

    parser.add_argument("-s", "--sample",
                        default=None,
                        type=str,
                        help="""[optional] - The name of the clone of interest for the plot of counts""")

    parser.add_argument("-N", "--Ns_count",
                        action="store_true",
                        help="""Specify if the Ns counts has to be included in the plot of the read counts""")

    parser.add_argument("-co", "--control_name",
                        type=str,
                        default=None,
                        help="""The name of the control group for the fold_change analysis, it should be the same 
                        name as the column name in the read_counts data structure, thus the name of the alignment 
                        file used as baseline without the ".bam" string""")

    parser.add_argument("-c", "--cigar_filter",
                        action="store_true",
                        default=None,
                        help="""If specified, the reads mapped with soft and hard clipping (S and H), are taken out 
                            form the read counts; it returns a data frame with same structure of the default one""")

    parser.add_argument("-cf", "--other_cigar_filters",
                        type=str,
                        nargs="+",
                        default=[],
                        help="""An additional parameter to exclude other reads from the count, on the bases of other 
                            information in their cigar, like indels.\n(Specify it like e.g. "I" "D")""")

    parser.add_argument("-fc", "--fold_change",
                        type=float,
                        default=1.5,
                        help="An integer to set the fold-change cut-off")

    parser.add_argument("-p", "--p_value",
                        type=float,
                        default=0.05,
                        help="An integer to set the p-value in order to retain the fold-change value"
                             "not picked up by chance")

    parser.add_argument("-u", "--unmapped",
                        action="store_true",
                        help="""If specified, also a .txt file is  created, with all the unmapped reads and, as 
                        last raw, the counts for each sample""")

    parser.add_argument("-sf", "--saving_folder",
                        type=str,
                        default="./plots/",
                        help="""Path to the directory in which create the plots folder and save all images; 
                        if not specified, a directory 'plots' will be created in the current one""")

    args = parser.parse_args()
    dict_args = vars(parser.parse_args([]))

    if args.flag_list != dict_args["flag_list"]:
        flags = dict_args["flag_list"] + args.flag_list
    else:
        flags = args.flag_list

    analyzer = TestingBinReadAnalyzer(args.folder, args.bin_size, args.reference, flags, args.output_pickle)
    analyzer.load_data(other_cigar_filters=args.other_cigar_filters,
                       cigar_filter=args.cigar_filter,
                       reference=args.reference,
                       read_info=args.read_info,
                       unmapped=args.unmapped,
                       verbose=True)

    if not os.path.exists(args.saving_folder):
        os.mkdir(args.saving_folder)

    analyzer.normalize_bins()
    analyzer.get_fold_change(args.saving_folder, args.control_name)
    analyzer.get_sig_pos(args.fold_change, args.p_value)
    analyzer.plot_sig_data(args.saving_folder, args.fold_change, args.p_value, args.cigar_filter)
    analyzer.plot_sh_clipping(args.control_name, args.saving_folder)

    exit(1)

    analyzer.normalize_bins()
    # analyzer.get_fold_change(args.control_name)
    fc = analyzer.get_fold_change(args.control_name)
    analyzer.get_sig_pos(args.fold_change, args.p_value)
    analyzer.get_no_sig_pos(args.fold_change, args.p_value)

    counts = params["read_counts"]
    start = 15309000 // args.bin_size
    # print(start)
    end = 15309706 // args.bin_size
    # print(end)

    data = counts[counts["chr"] == "CH.chr1"]
    data = data[data['bin'] >= start]
    data = data[data['bin'] <= end]
    print(data)
    # fc_start_filt = data["index"].iloc[0]
    # fc_end_filt = data["index"].iloc[-1]

    # for clone in fc:
    #     print(fc[clone].iloc[fc_start_filt:fc_end_filt + 1])

    # exit(1)



    analyzer.plot_counts_distributions(args.saving_folder, args.cigar_filter)

    if args.chromosome and args.sample:
        analyzer.plot_chrom_sample(args.saving_folder, args.reference, args.chromosome, args.sample, args.cigar_filter,
                                   args.Ns_count)
        analyzer.plot_norm_data_chr_sample(args.saving_folder, args.reference, args.chromosome, args.sample,
                                           args.cigar_filter,
                                           args.Ns_count)

    elif args.chromosome:
        analyzer.plot_chromosome(args.saving_folder, args.reference, args.chromosome, args.cigar_filter, args.Ns_count)
        analyzer.plot_norm_data_chr(args.saving_folder, args.reference, args.chromosome, args.cigar_filter,
                                    args.Ns_count)

    elif args.sample:
        analyzer.plot_sample(args.saving_folder, args.reference, args.sample, args.cigar_filter, args.Ns_count)
        analyzer.plot_norm_data_sample(args.saving_folder, args.reference, args.sample, args.cigar_filter,
                                       args.Ns_count)

    else:
        analyzer.plot_all(args.saving_folder, args.reference, args.cigar_filter, args.Ns_count)
        analyzer.plot_norm_data_all(args.saving_folder, args.reference, args.cigar_filter, args.Ns_count)

    # analyzer.plot_sig_data_chr_sample(args.saving_folder, args.fold_change, args.p_value, args.cigar_filter, args.chromosome, args.sample)

    # ---------------------------complete_file_pickle--------------------------------
    # with open("../all_samples_pickles/BRAN250000_df.p", "rb") as input_param:
    #     param = pickle.load(input_param)
    #     chr6 = param["read_counts"][param["read_counts"]["chr"] == "CH.chr6"]
    #     chr1 = param["read_counts"][param["read_counts"]["chr"] == "CH.chr1"]
    #     chr19 = param["read_counts"][param["read_counts"]["chr"] == "CH.chr19"]
    #     chr8 = param["read_counts"][param["read_counts"]["chr"] == "CH.chr8"]
    #     chr13 = param["read_counts"][param["read_counts"]["chr"] == "CH.chr13"]
    #
    #     # print(chr6[["Ch15_3_1_ref_IlluminaPE_aligns_Primary", "bin"]])
    #     print("all\n", param["read_counts"][["Ch15_3_1_ref_IlluminaPE_aligns_Primary", "Ch15_2_1_ref_IlluminaPE_aligns_Primary"]])
    #     print("chr6\n", chr6[["Ch15_3_1_ref_IlluminaPE_aligns_Primary", "Ch15_2_1_ref_IlluminaPE_aligns_Primary"]])
    #     print("chr1\n", chr1[["Ch15_3_1_ref_IlluminaPE_aligns_Primary", "Ch15_2_1_ref_IlluminaPE_aligns_Primary"]])
    #     print("chr19\n", chr19[["Ch15_3_1_ref_IlluminaPE_aligns_Primary", "Ch15_2_1_ref_IlluminaPE_aligns_Primary"]])
    #     print("chr8\n", chr8[["Ch15_3_1_ref_IlluminaPE_aligns_Primary", "Ch15_2_1_ref_IlluminaPE_aligns_Primary"]])
    #     print("chr13\n", chr13[["Ch15_3_1_ref_IlluminaPE_aligns_Primary", "Ch15_2_1_ref_IlluminaPE_aligns_Primary"]])
    #
    #     analyzer.normalize_bins(param["read_counts"], param["bin_size"])
    #     analyzer.get_fold_change(param["read_counts"], param["bin_size"], args.control_name)
    #     # analyzer.plot_all(param["read_counts"], param["bin_size"], args.reference)
    #     analyzer.plot_norm_data_all(param["read_counts"], param["bin_size"], args.reference)
    #     analyzer.get_sig_pos(param["read_counts"], param["bin_size"])
    #     analyzer.plot_sig_data(param["read_counts"], param["bin_size"])
    # analyzer.plot_background(param["read_counts"])
    # ------------------------------------------------------------------------------
