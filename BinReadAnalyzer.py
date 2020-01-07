# BinReadAnalyzer class
#
# This class aims to analyze the data structure created with the imported
# class: BinReadCounter
# One can consider this class as the very effector of analyses on plant of
# interest' genome. For analysis is intended, first of all,
# the normalization of the raw_read_counts and than the detection of
# significant differences in terms of number of reads per bin, against
# reference genome; this could lead to the identification of part of chromosome
# or even an entire chromosome that have a significantly different count of reads;
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
# plots are build up and shows as an interactive way via a web page;
# also an offline version of all the plots is saved in an appropriate folder

import os
import argparse
import plotly.graph_objects as go
import pandas as pd
from BinReadCounter import BinReadCounter
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

pandas2ri.activate()
base = importr('base')
utils = importr('utils')
stats = importr('stats')
edger = rpackages.importr('edgeR')


class BinReadAnalyzer:
    """This class provides different analysis of the data coming form
    BinReadCounter class; here imported.

    Options:
    -   plot counts distributions, normalized and log scaled
    -   plot read_counts for all chromosomes, for a specific one or for
        a specific sample
    -   normalize read_counts
    -   plot normalized_read_counts for all chromosomes, for a specific
        one or for a specific sample
    -   fold-change calculation
    -   plot of significant values in fold-change

    Attributes:
            bam_folder_path: the folder in which the necessary files are stored
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
        """set the parameters of pickle file"""
        self.parameters = param

    def set_norm(self, norm):
        """set the data structure of normalization counts"""
        self.norm = norm

    def set_fold_change(self, diff):
        """set the data structure of differences in counts"""
        self.fold_change = diff

    def set_sig_fc(self, fc):
        """set the data structure of significat differences in counts"""
        self.sig_fc = fc

    def load_data(self, other_cigar_filters, cigar_filter=None, reference=None,
                  read_info=None, unmapped=None, verbose=False):
        """If the pickle file is already present in the folder and the parameters
        are not changed, the data structures were kept, otherwise it runs the
        methods _export_pickle() and load_data() of the imported class BinReadCounter
        in order to create a new pickle file and to import the corresponding
        data-structures

        Args:
            other_cigar_filters (list): list of default filters + the optional ones
            cigar_filter (bool): if specified = True, the default filters are applied
            reference (bool): true if the reference is necessary
            read_info (bool): true if read info are necessary
            unmapped (bool): if specified != None, counts for unmapped reads are calculated
                            together with a file containing all unmapped reads
            verbose (bool): default is False, but if set to True, it allows to print
                            and see if a pickle file is already present or if the
                            Analyzer is running to create a new one

        Returns:
            A list of parameters and data structures from the pickle file or from the BinReadCounter;
        """
        counter = BinReadCounter(self.folder, self.bin_size, self.flags, self.ref, self.out)
        found = None
        bam_files = []
        i = 0
        list_bam_dir = os.listdir(self.folder)
        list_pick_dir = os.listdir(self.out)
        for file in list_bam_dir:
            if file.endswith(".bam"):
                bam_files.append(file[:file.find(".bam")])

        # verify if a pickle file is present in the desired folder
        for file in list_pick_dir:
            if file.endswith(".p"):
                # if so, checks if the parameters for the building-up of
                # the previous pickle file, are the same as the actual ones
                # with open(file, "rb") as input_pickle:
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
                        print(parameters)
                        found = True
                        if verbose:
                            print("Same parameters; import from: ", file)

                        self.set_parameters(parameters)
                        return self.parameters

                except:  # this allows to create/import different pickle files at the same time
                    i += 1
                    continue
        if verbose:
            print("Warning:\nException occurred - other processes running at this moment: ", i)

        if not found:
            # a file pickle is not present at all in the directory thus the algorithm
            # uses the modules of the BinReadCounter to calculate the new data structures and it saves a new pickle file
            if verbose:
                print("Parameters are changed or no pickle file exists",
                      "\nBinReadCounter is running with actual parameters",
                      "\nIT COULD TAKE A WHILE to create and import the pickle file")

            counter._export_pickle(other_cigar_filters, cigar_filter, reference, read_info, unmapped)
            name_pickle = counter.pickle_file_name(other_cigar_filters, cigar_filter, reference, read_info, unmapped)
            parameters = counter._load_pickle(name_pickle)

            self.set_parameters(parameters)
            return self.parameters

    def normalize_bins(self):
        """This method handles the normalization of the raw read_counts
        using an R package, edgeR, imported thanks to rpy2 that provides
        an already implemented function, cpm, for the normalization of a
        table of counts, as well as a series of other function specifically
        implemented for RNA-Seq"""
        read_counts = self.parameters["read_counts"]
        # the edgeR package is imported using rpy2 syntax to access to all its built-in functions
        read_counts_edger = {}  # a dictionary = sample: vector_of_counts to work with edger
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
        norm_counts = edger.cpm(norm_fact, normalized_lib_sizes=True)  # the log is calculated with fold-change method

        self.set_norm(norm_counts)
        return self.norm  # an edger object (table) of counts object (table) of counts

    def get_fold_change(self, control_name=None):
        """This method calculates the pairwise logFoldChange for the different clones
        against the control plant counts, in order to see, if are present, significant
        differences in terms of mapping reads"""
        if control_name is not None:
            read_counts = self.parameters["read_counts"]

            control = control_name
            read_counts_edger = {}  # a dictionary of sample: vector_of_counts to work with edger
            col_list = list(read_counts.columns)
            fold_change = {}

            bcv = 0.1  # reasonable dispersion value; typical values for common BCV (square-root dispersion) of data-set
            # arising from well-controlled experiments are 0.4 for human data, 0.1 for data of genetically identical
            # model organisms or 0.01 fo technical replicates

            for i in range(len(col_list)):
                if col_list[i] != "index" and col_list[i] != "chr" and col_list[i] != 'bin':
                    read_counts_edger[col_list[i]] = robjects.r.matrix(robjects.IntVector(read_counts[col_list[i]]))

            for i in range(1, self.norm.ncol + 1):
                col = self.norm.rx(True, i)
                col_name = self.norm.colnames[i - 1]
                if col_name != control:
                    print("\n", col_name)
                    tmp_dict = {col_name: col, control: self.norm.rx(True, control)}
                    tmp_df = robjects.DataFrame(tmp_dict)
                    dge_list = edger.DGEList(robjects.DataFrame(tmp_df), group=robjects.IntVector([2, 1]))

                    dge_list_ed = edger.estimateDisp(dge_list)
                    exact_test = edger.exactTest(dge_list_ed, dispersion=bcv ** 2)

                    exact_test_df = robjects.DataFrame(exact_test[0])

                    exact_test_pd_df = pandas2ri.ri2py_dataframe(exact_test_df)

                    fold_change[col_name] = exact_test_pd_df

            print(fold_change)
            self.set_fold_change(fold_change)
            return self.fold_change

        else:
            print("No control group is specified, please try again, specifying the column name of control file with "
                  "parameter '-co'")

    def get_sig_pos(self, fc, p_value):
        """This method return the chromosome position in terms of basis, of bins
        retained significantly different from the control = fold-change >= 1.5 and p-value <= 0.05"""
        fold_change = self.fold_change
        read_counts = self.parameters["read_counts"]
        sig_data = {"clone": [], "bin": [], "chr": [], "genome_position": [], "logFC": [], "PValue": []}

        for el in fold_change:
            sig_fc_values = fold_change[el][fold_change[el]["logFC"] >= fc]
            sig_fc_neg_values = fold_change[el][fold_change[el]["logFC"] <= -fc]
            sig_fc = pd.concat([sig_fc_values, sig_fc_neg_values])
            sig_fc_pvalues = sig_fc[sig_fc["PValue"] <= p_value]
            for i in sig_fc_pvalues.index.values:
                sig_data["clone"].append(el)
                sig_data["bin"].append(i)
                sig_data["chr"].append(read_counts["chr"].iloc[i])
                sig_data["genome_position"].append(i * self.bin_size)
                sig_data["logFC"].append(sig_fc_pvalues["logFC"][i])
                sig_data["PValue"].append(sig_fc_pvalues["PValue"][i])

        sig_data_df = pd.DataFrame(sig_data)
        self.set_sig_fc(sig_data_df)
        print(self.sig_fc)
        return self.sig_fc

    def plot_background(self, fig):
        """A trace tobe added to all scatter plots that allows the visualization of
        chromosomes on the image, in background"""
        read_counts = self.parameters["read_counts"]
        coordinates_x = []
        coordinates_y = []
        chromosomes = []

        for ch in read_counts["chr"]:
            if ch[ch.find("c"):] not in chromosomes:
                chromosomes.append(ch[ch.find("c"):])

        for chrom in chromosomes:
            single_df = read_counts[read_counts["chr"] == "CH." + chrom]
            coordinates_x.append(single_df["index"].mean())
            coordinates_y.append(-2.6)
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
                                              fillcolor="rgb(240, 248, 255)",  # ~light_mint_green
                                              opacity=0.5,
                                              layer="below",
                                              line_width=0))

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
                                      x1=len(read_counts),
                                      y1=fc))

        fig.add_shape(go.layout.Shape(type="line",
                                      x0=0,
                                      y0=-fc,
                                      x1=len(read_counts),
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
        """Histogram plots of counts, in order to see if the normalization was successful"""
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

        # fig_all.show()
        # fig_all_norm.show()
        # fig_log_norm.show()
        save_fig_all = fig_all.write_image(saving_folder + "/all_sample_dist_" + str(self.bin_size) + ".pdf",
                                           width=1280,
                                           height=1024)
        save_fig_all_norm = fig_all_norm.write_image(saving_folder + "/all_sample_norm_dist_" +
                                                     str(self.bin_size) + ".pdf",
                                                     width=1280,
                                                     height=1024)
        save_fig_log_norm = fig_log_norm.write_image(saving_folder + "/all_sample_norm_log_dist_" +
                                                     str(self.bin_size) + ".pdf",
                                                     width=1280,
                                                     height=1024)

    def plot_chrom_sample(self, saving_folder, reference, chrom, sample, cigar, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of raw read_counts
        for a specific chromosome and a specific sample

        Args:
            saving_folder (str): the path to where save the images
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
        c_name = None
        for c in list(read_counts["chr"]):
            if chrom in c:
                c_name = c
        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

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

        # fig.show()
        save_fig = fig.write_image(saving_folder + "/scatter_counts_chr" +
                                   str(chrom) + sample + str(self.bin_size) + ".pdf",
                                   width=1280,
                                   height=1024)

    def plot_chromosome(self, saving_folder, reference, chrom, cigar, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of raw_read_counts
        of all samples, but for a specific chromosome of interest

        Args:
            saving_folder (str): the path to where save the images
            reference (bool): true if the reference is declared
            chrom (int): a number representing the chromosome of interest
            cigar (list)
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of counts
        """
        read_counts = self.parameters["read_counts"]
        c_name = None
        for c in list(read_counts["chr"]):
            if chrom in c:
                c_name = c

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

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

        # fig.show()
        save_fig = fig.write_image(saving_folder + "/scatter_counts_chr" + str(chrom) + str(self.bin_size) + ".pdf",
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

        # fig.show()
        save_fig = fig.write_image(saving_folder + "/scatter_counts_" + sample + str(self.bin_size) + ".pdf",
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

        # fig.show()
        save_fig = fig.write_image(saving_folder + "/scatter_all_counts" + str(self.bin_size) + ".pdf",
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
        c_name = None

        for c in list(read_counts["chr"]):
            if chrom in c:
                c_name = c

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

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

        # fig.show()
        save_fig = fig.write_image(saving_folder + "/scatter__norm_counts_chr" +
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
        c_name = None
        for c in list(read_counts["chr"]):
            if str(chrom) in c:
                c_name = c

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

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

        # fig.show()
        save_fig = fig.write_image(
            saving_folder + "/scatter_norm_counts_chr" + str(chrom) + str(self.bin_size) + ".pdf",
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

        # fig.show()
        save_fig = fig.write_image(saving_folder + "/scatter_norm_counts_" + sample + str(self.bin_size) + ".pdf",
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
            fig.update_layout(
                title="Normalized Read Counts - Clone: all - Chr: all - Bin Size: " + str(self.bin_size),
                legend_orientation="h")

        # fig.show()
        save_fig = fig.write_image(saving_folder + "/scatter_norm_all_counts" + str(self.bin_size) + ".pdf",
                                   width=1280,
                                   height=1024)

    def plot_sig_data_chr_sample(self, saving_folder, fc, p_val, cigar, chrom, sample):
        pass

    def plot_sig_data_chr(self):
        pass

    def plot_sig_data_sample(self):
        pass

    def plot_sig_data(self, saving_folder, fc, p_val, cigar):  # , df_counts, bin_size
        read_counts = self.parameters["read_counts"]
        # read_counts = df_counts

        fold_change = self.fold_change
        sig_fc = self.sig_fc
        fig = go.Figure()

        fig.update_xaxes(title_text="Chromosomes_Bins")
        fig.update_yaxes(title_text="Fold-Change")

        for clone in sig_fc.clone.unique():
            sub_df = sig_fc[sig_fc["clone"] == clone]
            no_sig_fc = fold_change[clone].drop(list(sub_df["bin"]))

            fig.add_trace(go.Scatter(x=sub_df["bin"],
                                     y=sub_df["logFC"],
                                     mode="markers",
                                     name=clone))

            fig.add_trace(go.Scatter(x=no_sig_fc.index.values,
                                     y=no_sig_fc["logFC"],
                                     mode="markers",
                                     marker=dict(color="rgb(176, 196, 222)"),  # silver
                                     showlegend=False))

        self.add_threshold_fc(fig, fc)

        self.plot_background(fig)

        if cigar:
            fig.update_layout(title="Cigar Filter - Significant Fold Change Counts - Clone: all - Chr: all - "
                                    "Bin Size: " + str(self.bin_size) + " fc_threshold: " + str(fc) + " p_value: " +
                                    str(p_val),
                              legend_orientation="h")
        else:
            fig.update_layout(title="Significant Fold Change Counts - Clone: all - Chr: all - "
                                    "Bin Size: " + str(self.bin_size) + " fc_threshold: " + str(fc) + " p_value: " +
                                    str(p_val),
                              legend_orientation="h")

        # fig.show()
        save_fig = fig.write_image(saving_folder + "/counts_fold_change" + str(self.bin_size) + ".pdf",
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

    analyzer = BinReadAnalyzer(args.folder, args.bin_size, args.reference, flags, args.output_pickle)

    analyzer.load_data(other_cigar_filters=args.other_cigar_filters,
                       cigar_filter=args.cigar_filter,
                       reference=args.reference,
                       read_info=args.read_info,
                       unmapped=args.unmapped,
                       verbose=True)

    analyzer.normalize_bins()
    analyzer.get_fold_change(args.control_name)
    analyzer.get_sig_pos(args.fold_change, args.p_value)

    if not os.path.exists(args.saving_folder):
        os.mkdir(args.saving_folder)

    analyzer.plot_counts_distributions(args.saving_folder, args.cigar_filter)

    if args.chromosome and args.sample:
        analyzer.plot_chrom_sample(args.saving_folder, args.reference, args.chromosome, args.sample,
                                   args.cigar_filter, args.Ns_count)
        analyzer.plot_norm_data_chr_sample(args.saving_folder, args.reference, args.chromosome, args.sample,
                                           args.cigar_filter,
                                           args.Ns_count)

    elif args.chromosome:
        analyzer.plot_chromosome(args.saving_folder, args.reference, args.chromosome, args.cigar_filter,
                                 args.Ns_count)
        analyzer.plot_norm_data_chr(args.saving_folder, args.reference, args.chromosome, args.cigar_filter,
                                    args.Ns_count)

    elif args.sample:
        analyzer.plot_sample(args.saving_folder, args.reference, args.sample, args.cigar_filter, args.Ns_count)
        analyzer.plot_norm_data_sample(args.saving_folder, args.reference, args.sample, args.cigar_filter,
                                       args.Ns_count)

    else:
        analyzer.plot_all(args.saving_folder, args.reference, args.cigar_filter, args.Ns_count)
        analyzer.plot_norm_data_all(args.saving_folder, args.reference, args.cigar_filter, args.Ns_count)

    analyzer.plot_sig_data(args.saving_folder, args.fold_change, args.p_value, args.cigar_filter)