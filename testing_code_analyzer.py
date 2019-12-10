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
import pickle
import plotly.graph_objects as go
# import plotly.express as px  # plot norm_counts only one sample
# import plotly.figure_factory as ff  # plot norm_counts multiple sample
import pandas as pd
# import numpy as np
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

    def set_parameters(self, param):
        self.parameters = param

    def set_norm(self, norm):
        self.norm = norm

    def set_fold_change(self, diff):
        self.fold_change = diff

    def load_data(self, other_cigar_filters, cigar_filter=None, reference=None, read_info=None, verbose=False):
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
        found = False
        bam_files = []
        list_bam_dir = os.listdir(self.folder)
        list_pick_dir = os.listdir(self.out)
        # verify if a pickle file is present in the desired folder
        for file in list_bam_dir:
            if file.endswith(".bam"):
                bam_files.append(file[:file.find(".bam")])

        for file in list_pick_dir:
            if file.endswith(".p"):
                # if so, checks if the parameters for the building-up of
                # the previous pickle file, are the same as the actual ones
                # with open(file, "rb") as input_pickle:
                #     print("ok")
                #     counter = pickle.load(input_pickle)
                parameters = counter._load_pickle(file)
                # print(parameters)
                if parameters["bin_size"] == self.bin_size and \
                        parameters["flags"] == self.flags and \
                        parameters["bam"] == bam_files and \
                        parameters["cigar_filt"] == cigar_filter and \
                        parameters["other_cigar_filt"] == other_cigar_filters and \
                        parameters["ref"] == counter.get_ref_name() and \
                        parameters["info"] == read_info:
                    found = True
                    if verbose:
                        print("Same parameters; import from: ", file)
                    self.set_parameters(parameters)
                    return self.parameters
        if not found:
            # if not found, none of the pickle files in the current directory
            # have the same parameters of the actual running or a file pickle is not p
            # resent at all in the directory thus the algorithm uses the modules of the
            # BinReadCounter to calculate the new data structures and it saves a new pickle file
            if verbose:
                print(" Parameters are changed or no pickle file exists\n",
                      "BinReadCounter is running with actual parameters\n",
                      "IT COULD TAKE A WHILE to create and import the pickle file")

            counter._export_pickle(other_cigar_filters, cigar_filter, reference, read_info)
            name_pickle = counter.pickle_file_name(other_cigar_filters, cigar_filter, reference, read_info)
            parameters = counter._load_pickle(name_pickle)
            self.set_parameters(parameters)
            return self.parameters

    def normalize_bins(self):  # , df_counts
        """This method handles the normalization of the raw read_counts
        using an R package, edgeR, imported thanks to rpy2 that provides
        an already implemented function, cpm, for the normalization of a
        table of counts, as well as a series of other function specifically
        implemented for RNA-Seq"""
        read_counts = self.parameters["read_counts"]
        # read_counts = df_counts

        # fig = px.histogram(read_counts, x="FEM_ref_IlluminaPE_aligns_Primary")
        # fig.show()
        # counts_to_plot = []
        # col_labels = []
        # for col in read_counts.columns:
        #     if col != "index" and col != "chr" and col != 'bin':
        #         counts_to_plot.append(read_counts[col])
        #         col_labels.append(col)
        # fig = ff.create_distplot(counts_to_plot, group_labels=col_labels)
        # fig.show()

        # the following line has to be deleted when the class will work on the entire .bam files
        # chr1 = read_counts[read_counts['chr'] == 'CH.chr1']
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

        # create a DGEList
        dge_list = edger.DGEList(counts=read_counts_edger_df, group=robjects.IntVector(group))
        # print(dge_list)

        # calc norm factors
        norm_fact = edger.calcNormFactors(dge_list)
        # print(norm_fact)

        # normalization function
        # norm_counts = edger.cpm(read_counts_edger_df, log=True)
        norm_counts = edger.cpm(norm_fact, normalized_lib_sizes=True)
        print(norm_counts.colnames)
        # print(norm_counts)
        # norm_counts_rdf = robjects.DataFrame(norm_counts)
        # norm_counts_df = pandas2ri.ri2py_dataframe(norm_counts_rdf)
        norm_counts_dict = {}
        for i in range(1, norm_counts.ncol + 1):
            norm_counts_dict[norm_counts.colnames[i-1]] = list(norm_counts.rx(True, i))
        norm_counts_df = pd.DataFrame(norm_counts_dict)
        print(norm_counts_df)

        # plot only one sample

        # norm_counts_df["index"] = read_counts["index"]
        # norm_counts_df["Ch15_3_1_ref_IlluminaPE_aligns_Primary"] =
        # np.log(norm_counts_df["Ch15_3_1_ref_IlluminaPE_aligns_Primary"]) fig = px.histogram(norm_counts_df,
        # x="Ch15_3_1_ref_IlluminaPE_aligns_Primary") fig.show()

        # plot all samples
        # norm_counts_to_plot = []
        # for col in norm_counts_df.columns:
        #     norm_counts_to_plot.append(norm_counts_df[col])
        # fig_norm = ff.create_distplot(norm_counts_to_plot, group_labels=norm_counts_df.columns)
        # fig_norm = ff.create_distplot(norm_counts_to_plot, group_labels=norm_counts_df.columns, bin_size=0.05)
        # fig_norm.show()

        # print(norm_counts_df)
        # print(list(norm_counts.rx(True, 1)))
        self.set_norm(norm_counts)
        return self.norm  # an edger object (table) of counts

    def get_fold_change(self, control_name=None):  # df_counts,
        """This method calculates the pairwise logFoldChange for the different clones
        against the control plant counts, in order to see, if are present, significant
        differences in terms of mapping reads"""
        if control_name is not None:
            read_counts = self.parameters["read_counts"]
            # read_counts = df_counts
            control = control_name
            read_counts_edger = {}  # {}  # a dictionary of sample: vector_of_counts to work with edger
            bcv = 0.1  # reasonable dispersion value; typical values for common BCV (square-root dispersion) of data-set
            # arising from well-controlled experiments are 0.4 for human data, 0.1 for data of genetically identical
            # model organisms or 0.01 fo technical replicates
            col_list = list(read_counts.columns)
            fold_change = {}
            for i in range(len(col_list)):
                if col_list[i] != "index" and col_list[i] != "chr" and col_list[i] != 'bin':
                    read_counts_edger[col_list[i]] = robjects.r.matrix(robjects.IntVector(read_counts[col_list[i]]))

            # with raw data
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

            # with normalized data
            for i in range(1, self.norm.ncol + 1):
                col = self.norm.rx(True, i)
                col_name = self.norm.colnames[i-1]
                if col_name != control:
                    print(col_name)
                    tmp_dict = {col_name: col, control: self.norm.rx(True, control)}
                    tmp_df = robjects.DataFrame(tmp_dict)
                    dge_list = edger.DGEList(robjects.DataFrame(tmp_df), group=robjects.IntVector([2, 1]))
                    # print(dge_list)
                    dge_list_ed = edger.estimateDisp(dge_list)
                    exact_test = edger.exactTest(dge_list_ed, dispersion=bcv**2)
                    exact_test_df = robjects.DataFrame(exact_test[0])
                    # print(type(exact_test_df))
                    exact_test_pd_df = pandas2ri.ri2py_dataframe(exact_test_df)
                    # print(exact_test_pd_df)
                    fold_change[col_name] = exact_test_pd_df

                    # print(edger.topTags(exact_test))
            # print(fold_change)
            self.set_fold_change(fold_change)
            return self.fold_change

        else:
            print("No control group is specified, please try again, specifying the column name of control file as "
                  "parameter '-co'")

    def get_sig_pos(self):  # , df_counts
        """This method return the chromosome position in terms of basis, of bins
        retained significantly different from the control"""
        fold_change = self.fold_change
        read_count = self.parameters["read_counts"]
        # read_counts = df_counts
        sig_data = {"clone": [], "bin": [], "chr": [], "genome_position": []}

        print(fold_change)
        for el in fold_change:
            sig_fc = fold_change[el][fold_change[el]["logFC"] > 1.5]
            sig_row_index = sig_fc.index.values
            for i in range(len(sig_row_index)):
                sig_data["clone"].append(el)
                sig_data["bin"].append(read_count["bin"].iloc[sig_row_index[i]])
                sig_data["chr"].append(read_count["chr"].iloc[sig_row_index[i]])
                sig_data["genome_position"].append(sig_row_index[i] * self.parameters["bin_size"])
        sig_data_df = pd.DataFrame(sig_data)
        print(sig_data_df)
        #     print(sig_fc.index.values)
        #     print(fc_data)
        #     if fc_data[fc_data > 1.5]:
        #         print(fc_data[fc_data["logFC"] > 1.5])
        # here I have to retrieve information about th position of the bin in the chromosome
        return sig_data_df

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

    def plot_chrom_sample(self, reference, chrom, sample, ns=False, fig=go.Figure()):
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

        single_chrom = read_counts[read_counts["chr"] == "CH.chr" + str(chrom)]
        col_list = list(single_chrom.columns)
        for i in range(len(col_list[:col_list.index(sample)]) + 1):
            if col_list[i] == sample:
                fig.add_trace(go.Scatter(x=single_chrom["index"],
                                         y=single_chrom[sample],
                                         mode="markers",
                                         name=str(col_list[i])))

        if ns:
            self.add_ns_trace(fig, reference=reference, chrom=chrom)

        fig.update_layout(
            title="Read Counts - Clone: " + sample + " - Chr: " + str(chrom) + " - Bin Size: " + str(self.bin_size))

        fig.show()

    def plot_chromosome(self, reference, chrom, ns=False, fig=go.Figure()):
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

        single_chrom = read_counts[read_counts["chr"] == "CH.chr" + str(chrom)]
        col_list = list(single_chrom.columns)
        for i in range(len(col_list)):
            if col_list[i] != "index" and col_list[i] != "chr" and col_list[i] != "bin":
                fig.add_trace(go.Scatter(x=single_chrom["index"],
                                         y=single_chrom[col_list[i]],
                                         mode="markers",
                                         name=str(col_list[i])))
        if ns:
            self.add_ns_trace(fig, reference=reference, chrom=chrom)

        fig.update_layout(title="Read Counts - All Clones - Chr: " + str(chrom) + " - Bin Size: " + str(self.bin_size))

        fig.show()

    def plot_sample(self, reference, sample, ns=False, fig=go.Figure()):
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

        fig.update_layout(title="Read Counts - Clone: " + sample + " - Bin Size: " + str(self.bin_size))

        fig.show()

    def plot_all(self, reference, ns=False, fig=go.Figure()):  # df_counts,
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
        coordinates_x = []
        coordinates_y = []

        chromosomes = []
        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        col_list = list(read_counts.columns)
        for i in range(len(col_list)):
            if col_list[i] != "index" and col_list[i] != "chr" and col_list[i] != "bin":
                fig.add_trace(go.Scatter(x=read_counts["index"],
                                         y=read_counts[col_list[i]],
                                         mode="markers",
                                         name=str(col_list[i])))

        for ch in read_counts["chr"]:
            if ch[ch.find("c"):] not in chromosomes:
                chromosomes.append(ch[ch.find("c"):])

        for chrom in chromosomes:
            single_df = read_counts[read_counts["chr"] == "CH." + chrom]
            coordinates_x.append(single_df["index"].mean())
            coordinates_y.append(-1)
            if int(chrom[chrom.find("r") + 1:]) % 2 == 0:
                fig.add_shape(go.layout.Shape(type="rect",
                                              xref="x",
                                              yref="paper",
                                              x0=single_df["index"].iloc[0],
                                              y0=0,
                                              x1=single_df["index"].iloc[-1],
                                              y1=1,
                                              fillcolor="rgb(230, 230, 250)",  # lavande
                                              opacity=0.5,
                                              layer="below",
                                              line_width=0))
            else:
                fig.add_shape(go.layout.Shape(type="rect",
                                              xref="x",
                                              yref="paper",
                                              x0=single_df["index"].iloc[0],
                                              y0=0,
                                              x1=single_df["index"].iloc[-1],
                                              y1=1,
                                              fillcolor="rgb(240, 248, 255)",  # ~light_green
                                              opacity=0.5,
                                              layer="below",
                                              line_width=0))

        fig.add_trace(go.Scatter(x=coordinates_x,
                                 y=coordinates_y,
                                 text=chromosomes,
                                 mode="text",
                                 showlegend=False,
                                 ))

        if ns:
            self.add_ns_trace(fig, reference=reference)

        fig.update_layout(title="Read Counts - All Clones - All Chromosomes - Bin Size: " + str(self.bin_size))

        fig.show()

    def plot_norm_data_chr_sample(self, reference, chrom, sample, ns=False, fig=go.Figure()):
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

        single_chrom = read_counts[read_counts["chr"] == "CH.chr" + str(chrom)]
        for i in range(0, norm_counts.ncol):
            if norm_counts.colnames[i] == sample:
                fig.add_trace(go.Scatter(x=single_chrom["index"],
                                         y=list(norm_counts.rx(True, i + 1)[single_chrom["index"].iloc[0]:
                                                                            single_chrom["index"].iloc[-1] + 1]),
                                         mode="markers",
                                         name=str(norm_counts.colnames[i])))
        if ns:
            self.add_ns_trace(fig, reference=reference, chrom=chrom)

        fig.update_layout(title="Normalized Read Counts - Clone: " +
                                sample +
                                " - Chr: " +
                                str(chrom) +
                                " - Bin Size: " +
                                str(self.bin_size))

        fig.show()

    def plot_norm_data_chr(self, reference, chrom, ns=False, fig=go.Figure()):
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

        single_chrom = read_counts[read_counts["chr"] == "CH.chr" + str(chrom)]
        for i in range(0, norm_counts.ncol):
            fig.add_trace(go.Scatter(x=single_chrom["index"],
                                     y=list(norm_counts.rx(True, i + 1)[single_chrom["index"].iloc[0]:
                                                                        single_chrom["index"].iloc[-1] + 1]),
                                     mode="markers",
                                     name=norm_counts.colnames[i]))
        if ns:
            self.add_ns_trace(fig, reference=reference, chrom=chrom)

        fig.update_layout(title="Normalized Read Counts - Clone: all - Chr: " +
                                str(chrom) +
                                " - Bin Size: " +
                                str(self.bin_size))

        fig.show()

    def plot_norm_data_sample(self, reference, sample, ns=False, fig=go.Figure()):
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

        fig.update_layout(title="Normalized Read Counts - Clone: " +
                                sample +
                                " - Chr: all - Bin Size: " +
                                str(self.bin_size))

        fig.show()

    def plot_norm_data_all(self, reference, ns=False, fig=go.Figure()):  # df_counts,
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
        coordinates_x = []
        coordinates_y = []
        chromosomes = []
        norm_counts = self.norm
        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")
        for i in range(0, norm_counts.ncol):
            fig.add_trace(go.Scatter(x=read_counts["index"],
                                     y=list(norm_counts.rx(True, i + 1)),
                                     mode="markers",
                                     name=norm_counts.colnames[i]))
        for ch in read_counts["chr"]:
            if ch[ch.find("c"):] not in chromosomes:
                chromosomes.append(ch[ch.find("c"):])

        for chrom in chromosomes:
            single_df = read_counts[read_counts["chr"] == "CH." + chrom]
            coordinates_x.append(single_df["index"].mean())
            coordinates_y.append(-2)
            if int(chrom[chrom.find("r") + 1:]) % 2 == 0:
                fig.add_shape(go.layout.Shape(type="rect",
                                              xref="x",
                                              yref="paper",
                                              x0=single_df["index"].iloc[0],
                                              y0=0,
                                              x1=single_df["index"].iloc[-1],
                                              y1=1,
                                              fillcolor="rgb(230, 230, 250)",
                                              opacity=0.5,
                                              layer="below",
                                              line_width=0))
            else:
                fig.add_shape(go.layout.Shape(type="rect",
                                              xref="x",
                                              yref="paper",
                                              x0=single_df["index"].iloc[0],
                                              y0=0,
                                              x1=single_df["index"].iloc[-1],
                                              y1=1,
                                              fillcolor="rgb(240, 248, 255)",
                                              opacity=0.5,
                                              layer="below",
                                              line_width=0))

        fig.add_trace(go.Scatter(x=coordinates_x,
                                 y=coordinates_y,
                                 text=chromosomes,
                                 mode="text",
                                 showlegend=False,
                                 ))

        if ns:
            self.add_ns_trace(fig, reference=reference)

        fig.update_layout(title="Normalized Read Counts - Clone: all - Chr: all - Bin Size: " + str(self.bin_size))

        fig.show()

    def plot_sig_data(self, df_counts, bs, fig=go.Figure()):  #  df_counts,
        read_counts = df_counts
        # read_counts = self.parameters["read_counts"]
        fold_change = self.fold_change
        coordinates_x = []
        coordinates_y = []
        chromosomes = []
        bin_size = bs
        for el in fold_change:
            fold_change[el]["index"] = read_counts["index"]
            # threshold = fold_change[el][fold_change[el]["logFC"] > 1.5]
            # print(threshold)
            fig.update_xaxes(title_text="Chromosomes_Bins")
            fig.update_yaxes(title_text="Read_Count_Per_Bin")
            fig.add_trace(go.Scatter(x=fold_change[el]["index"],
                                     y=fold_change[el]["logFC"],
                                     mode="markers",
                                     name=el))

            # for every bin with significant difference form the control, plot a dashed-dot line
            # fig.add_shape(go.layout.Shape(type="line",
            #                               x0=threshold["index"],
            #                               y0=threshold["logFC"],
            #                               x1=threshold["index"],
            #                               y1=threshold["logFC"] + 10,
            #                               line=dict(color="black",
            #                                         width=3,
            #                                         dash="dashdot")))

        for ch in read_counts["chr"]:
            if ch[ch.find("c"):] not in chromosomes:
                chromosomes.append(ch[ch.find("c"):])

        for chrom in chromosomes:
            single_df = read_counts[read_counts["chr"] == "CH." + chrom]
            coordinates_x.append(single_df["index"].mean())
            coordinates_y.append(-1)
            if int(chrom[chrom.find("r") + 1:]) % 2 == 0:
                fig.add_shape(go.layout.Shape(type="rect",
                                              xref="x",
                                              yref="paper",
                                              x0=single_df["index"].iloc[0],
                                              y0=0,
                                              x1=single_df["index"].iloc[-1],
                                              y1=1,
                                              fillcolor="rgb(230, 230, 250)",
                                              opacity=0.5,
                                              layer="below",
                                              line_width=0))
            else:
                fig.add_shape(go.layout.Shape(type="rect",
                                              xref="x",
                                              yref="paper",
                                              x0=single_df["index"].iloc[0],
                                              y0=0,
                                              x1=single_df["index"].iloc[-1],
                                              y1=1,
                                              fillcolor="rgb(240, 248, 255)",
                                              opacity=0.5,
                                              layer="below",
                                              line_width=0))

        fig.add_trace(go.Scatter(x=coordinates_x,
                                 y=coordinates_y,
                                 text=chromosomes,
                                 mode="text",
                                 showlegend=False,
                                 ))

        # fig.update_layout(title="Read Counts - All Clones - All Chromosomes - Bin Size: " + str(self.bin_size))
        fig.update_layout(title="Read Counts - All Clones - All Chromosomes - Bin Size: " + str(bin_size))

        fig.show()


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
                        type=int,
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
                       verbose=True)

    # if not os.path.exists("plots"):
    #     os.mkdir("plots")
    analyzer.normalize_bins()
    analyzer.get_fold_change(args.control_name)
    analyzer.get_sig_pos()

    # if args.chromosome and args.sample:
    #     if args.Ns_count:
    #         analyzer.plot_chrom_sample(args.reference, args.chromosome, args.sample, args.Ns_count)
    #         analyzer.plot_norm_data_chr_sample(args.reference, args.chromosome, args.sample, args.Ns_count)
    #     else:
    #         analyzer.plot_chrom_sample(args.reference, args.chromosome, args.sample)
    #         analyzer.plot_norm_data_chr_sample(args.reference, args.chromosome, args.sample)
    #
    # elif args.chromosome:
    #     if args.Ns_count:
    #         analyzer.plot_chromosome(args.reference, args.chromosome, args.Ns_count)
    #         analyzer.plot_norm_data_chr(args.reference, args.chromosome, args.Ns_count)
    #     else:
    #         analyzer.plot_chromosome(args.reference, args.chromosome)
    #         analyzer.plot_norm_data_chr(args.reference, args.chromosome)
    #
    # elif args.sample:
    #     if args.Ns_count:
    #         analyzer.plot_sample(args.reference, args.sample, args.Ns_count)
    #         analyzer.plot_norm_data_sample(args.reference, args.sample, args.Ns_count)
    #     else:
    #         analyzer.plot_sample(args.reference, args.sample)
    #         analyzer.plot_norm_data_sample(args.reference, args.sample)
    #
    # else:
    #     if args.Ns_count:
    #         analyzer.plot_all(args.reference, args.Ns_count)
    #         analyzer.plot_norm_data_all(args.reference, args.Ns_count)
    #     else:
    #         analyzer.plot_all(args.reference)
    #         analyzer.plot_norm_data_all(args.reference)
    #         analyzer.plot_sig_data(args.bin_size)

    # ---------------------------complete_file_pickle--------------------------------
    # with open("/home/giulia/Downloads/BRAN30000_df_c.p", "rb") as input_param:
    #     param = pickle.load(input_param)
    #     print(param)
    #     analyzer.normalize_bins(param["read_counts"])
    #     analyzer.get_fold_change(param["read_counts"], args.control_name)
    #     # # analyzer.plot_all(param["read_counts"], args.reference)
    #     # # analyzer.plot_norm_data_all(param["read_counts"], args.reference)
    #     analyzer.plot_sig_data(param["read_counts"], param["bin_size"])
    #     # # analyzer.get_sig_pos(param["read_counts"])
    # ------------------------------------------------------------------------------
