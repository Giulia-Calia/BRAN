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
import math
# import pickle
import re
import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from BinReadCounter import BinReadCounter
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from plotly.subplots import make_subplots

# import plotly
# import time
# plotly.io.orca.ensure_server()
# time.sleep(10)

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
    -   normalize read_counts calculation
    -   fold-change calculation
    -   plot read_counts for all chromosomes, for a specific one or for
        a specific sample
    -   plot normalized_read_counts for all chromosomes, for a specific
        one or for a specific sample
    -   plot fold-change for all chromosomes, for a specific one or for
        a specific sample, taking in account the pairwise (-pw) parameter

    Attributes:
            folder_path: path to the folder in which the .bam files are
                         stored
            bin_size: bins' length in which each chromosome is divided
                      (in bp)
            flags: a list of string, each corresponding to a particular
                   pairwise-flag in the .bam/.sam format
            ref_genome: a string being the name of the reference_genome file
            cigar_filter: a filter of reads with a specific cigar, specified
                          by the parameter -c gives access to additional graphs
            out_pickle: default is None because it takes the name of an
                        already existing file or the name of a new file
                        created during the computation
            parameters: all the parameters of the already existed/new pickle
                        file
            norm: data structure build up by the normalization function
                  of the class
            log_norm: log2 applied to normalized counts
            norm_clip: normalized clipped read counts
            log_norm_clip: log2 applied to normalized clipped counts
            norm_unmapped:
            fold_change:
            clipped_fold_change:
            sig_data:
            clip_sig_data:

    """

    def __init__(self, bam_folder_path, bin_size, ref_genome, flags, cigar_filter, out_pickle):
        self.folder = bam_folder_path
        self.bin_size = bin_size
        self.flags = flags
        self.ref = ref_genome
        self.cigar_filter = cigar_filter
        self.out = out_pickle
        self.parameters = None
        self.norm = None
        self.log_norm = None
        self.norm_clip = None
        self.log_norm_clip = None
        self.norm_unmapped = None
        self.fold_change = None
        self.clipped_fold_change = None
        self.sig_data = None
        self.clip_sig_data = None

    def set_parameters(self, param):
        """set parameters as actual parameters"""
        self.parameters = param

    def set_read_counts(self, sort_read_counts):
        self.parameters["read_counts"] = sort_read_counts

    def set_norm(self, norm):
        """set norm counts from normalization function"""
        self.norm = norm

    def set_log_norm(self, log_norm):
        """set log2 of norm counts"""
        self.log_norm = log_norm

    def set_norm_clip(self, norm_clip):
        """set normalized clipped counts from normalization"""
        self.norm_clip = norm_clip

    def set_log_norm_clip(self, log_norm_clip):
        """set log2 of norm_clip_counts"""
        self.log_norm_clip = log_norm_clip

    def set_norm_unmapped(self, norm_unmapped):
        """set normalized counts of unmapped reads"""
        self.norm_unmapped = norm_unmapped

    def set_fold_change(self, diff):
        """set fold-change for read counts"""
        self.fold_change = diff

    def set_clipped_fold_change(self, clip_fc):
        """set fold-change for clipped read counts"""
        self.clipped_fold_change = clip_fc

    def set_sig_data(self, sig_data):
        """set the data structure of significant read fold-change counts"""
        self.sig_data = sig_data

    def set_clip_sig_data(self, clip_sig_data):
        """set the data structure of significant clipped fold-change counts"""
        self.clip_sig_data = clip_sig_data

    def load_data(self, cigar, cigar_filter=None, reference=None, read_info=None, unmapped=None,
                  verbose=False):
        """If the pickle file is already present in the output pickle folder and the
        parameters are not changed, the data structures is kept, otherwise the
        methods _export_pickle() and load_data() of the imported class BinReadCounter are
        used in order to create a new pickle file and to import the corresponding
        data-structures

        Args:
            cigar(bool): true if the parameter -c is given
            cigar_filter (str): other cigar filter other than default ones [S, H]
            reference (bool): true if the reference is necessary
            read_info (bool): true if read info are necessary
            unmapped (bool): true if the unmapped dictionary is given
            verbose (bool): default is False, but if set to True, it allows to print
                            and see if a pickle file is already present or if the
                            Analyzer is running to create a new one

        Returns:
            a dictionary which includes: bin_size, flags, list of .bam file names,
            length of each chromosome, cigar:True/False, cigar_filters, read_counts
            data structure, unmapped: True/False, unmapped read counts, reference name
            for N count, N_counts, read_info: True/False
        """
        counter = BinReadCounter(self.folder, self.bin_size, self.flags, self.ref, self.cigar_filter, self.out)
        found = None
        bam_files = []
        i = 0
        list_pick_dir = os.listdir(self.out)
        for f in self.folder:
            list_bam_dir = os.listdir(f)
            for file in list_bam_dir:
                if file.endswith(".bam"):
                    bam_files.append(file[:file.find(".bam")])
        # verify if a pickle file is present in the desired folder
        for file in list_pick_dir:
            if file.endswith(".p"):
                # if so, checks if the parameters for the building-up of
                # the previous pickle file, are the same as the actual ones
                parameters = counter._load_pickle(file)
                if parameters["bin_size"] == self.bin_size and \
                        parameters["flags"] == self.flags and \
                        all(el in parameters["bam"] for el in bam_files) and \
                        all(el in bam_files for el in parameters["bam"]) and \
                        parameters["cigar"] == cigar and \
                        parameters["cigar_filter"] == cigar_filter and \
                        parameters["unmapped"] == unmapped and \
                        parameters["ref"] == counter.get_ref_name() and \
                        parameters["info"] == read_info:
                    print(parameters)
                    found = True
                    if verbose:
                        print("Same parameters; import from: ", file)

                    if parameters["info"]:
                        counter._load_read_ID(cigar)

                    self.set_parameters(parameters)
                    return self.parameters

                else:
                    continue

        if not found:
            # if not found, none of the pickle files in the current directory
            # have the same parameters of the actual running or a file pickle is not p
            # resent at all in the directory thus the algorithm uses the modules of the
            # BinReadCounter to calculate the new data structures and it saves a new pickle file

            if verbose:
                print("No other pickle files exist in this directory or none has the same actual parameters",
                      "\nBinReadCounter is running with actual parameters",
                      "\nIT COULD TAKE A WHILE to create and import the pickle file")

            counter._export_pickle(cigar, reference, read_info, unmapped)
            name_pickle = counter.pickle_file_name(cigar, reference, read_info, unmapped)
            parameters = counter._load_pickle(name_pickle)

            if parameters["info"]:
                counter._load_read_ID(cigar)

            self.set_parameters(parameters)
            print(self.parameters)
            return self.parameters

    def sorted_chromosomes(self, column):
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(column, key=alphanum_key)

    def no_repeats_df_transformation(self):
        read_counts = self.parameters["read_counts"]
        # print(read_counts)
        chr_list = list(read_counts["chr"].value_counts().index)
        sorted_chr_list = self.sorted_chromosomes(chr_list)

        sorted_df = pd.DataFrame()
        for chrom in sorted_chr_list:
            single_df = read_counts[read_counts["chr"] == chrom]
            sorted_df = pd.concat([sorted_df, single_df])
        sorted_df = sorted_df.reset_index(drop=True)
        self.set_read_counts(sorted_df)
        # print("sorted_read_counts: \n", self.parameters["read_counts"])
        sorted_df.to_csv("sorted_df.txt", sep="\t")
        return self.parameters["read_counts"]

    def normalize_bins(self, control_name):
        """This method handles the normalization of the raw read_counts
        using an R package, edgeR, imported thanks to rpy2 that provides
        an already implemented function, cpm, for the normalization of a
        table of counts, as well as a series of other function specifically
        implemented for RNA-Seq; the normalization for clipped counts is done
        manually"""
        read_counts = self.parameters["read_counts"]
        read_counts.to_csv(str(self.bin_size) + "read_counts_check.txt", sep="\t")
        # read_counts = self.no_repeats_df_transformation()
        # the edgeR package is imported using rpy2 syntax to access to all its built-in functions
        read_counts_edger = {}  # a dictionary of sample: vector_of_counts to work with edger
        clipped_count = {}

        for col in read_counts:
            if col != "chr" and col != 'bin' and "cig_filt" not in col:
                # creates an element of the dictionary [sample_name]: r_vector_of_counts
                read_counts_edger[col] = robjects.IntVector(read_counts[col])

            elif "cig_filt" in col:
                clipped_count[col] = read_counts[col]

        read_counts_edger_df = robjects.DataFrame(read_counts_edger)  # R data frame
        norm_counts = edger.cpm(read_counts_edger_df, normalized_lib_sizes=True)
        log_norm_counts = edger.cpm(read_counts_edger_df, log=True)

        norm_counts_dict = {}
        log_norm_counts_dict = {}

        for i in range(1, norm_counts.ncol + 1):
            norm_counts_dict[norm_counts.colnames[i - 1]] = list(norm_counts.rx(True, i))
            log_norm_counts_dict[log_norm_counts.colnames[i - 1]] = list(log_norm_counts.rx(True, i))

        norm_counts_df = pd.DataFrame(norm_counts_dict)
        norm_counts_df = pd.concat([read_counts["chr"], read_counts['bin'], norm_counts_df], axis=1)
        # print(read_counts[["chr", "bin"]])
        # print(norm_counts_df)
        log_norm_counts_df = pd.DataFrame(log_norm_counts_dict)
        log_norm_counts_df = pd.concat([read_counts["chr"], read_counts['bin'], log_norm_counts_df], axis=1)

        self.set_norm(norm_counts_df)
        self.set_log_norm(log_norm_counts_df)
        # print(self.norm)

        clipped_count_df = pd.DataFrame(clipped_count)
        norm_clip = {}
        log_norm_clip = {}

        for col in clipped_count_df:
            read_counts_col = col[:col.find("_cig_filt")]
            norm_clip[col] = clipped_count_df[col] / (sum(read_counts[read_counts_col]) / 1000000)
            log_norm_clip[col] = []
            for row in clipped_count_df[col]:
                if row == 0:
                    log_norm_clip[col].append(0)
                else:
                    log_norm_clip[col].append(math.log2(row))

        norm_clip_df = pd.DataFrame(norm_clip)
        log_norm_clip_df = pd.DataFrame(log_norm_clip)
        log_norm_clip_df = pd.concat([read_counts["chr"], read_counts['bin'], log_norm_clip_df], axis=1)

        self.set_norm_clip(norm_clip_df)
        self.set_log_norm_clip(log_norm_clip_df)

        unmapped = self.parameters["unmapped_reads"]
        norm_unmapped = {}
        for el in unmapped:
            norm_unmapped[el] = unmapped[el] / (sum(read_counts[el]) / 1000000)

        self.set_norm_unmapped(norm_unmapped)

        return self.norm, self.log_norm, self.norm_clip, self.log_norm_clip, self.norm_unmapped

    def calc_fold_change(self, control_name, pairwise=False):
        # fc = self.log_norm["test6_alignSort_REDONE"] - self.log_norm["reference30x_alignSort"]

        if pairwise:
            fc_read_counts = {"chr": self.log_norm["chr"], "bin": self.log_norm["bin"]}
            for col in self.log_norm:
                if col != control_name and col != "bin" and col != "chr":
                    pw_fc = self.log_norm[col] - self.log_norm[control_name]
                    fc_read_counts[col + "-" + control_name] = list(pw_fc)

            fc_clipped_counts = {"chr": self.log_norm_clip["chr"], "bin": self.log_norm_clip["bin"]}
            for col in self.log_norm_clip:
                if col != control_name + "_cig_filt" and col != "bin" and col != "chr":
                    pw_clipped_fc = self.log_norm_clip[col] - self.log_norm_clip[control_name + "_cig_filt"]
                    fc_clipped_counts[col + "-" + control_name] = list(pw_clipped_fc)

            fc_df = pd.DataFrame(fc_read_counts)
            fc_clip_df = pd.DataFrame(fc_clipped_counts)
            print(fc_clip_df)
            # print(fc_df)
            self.set_fold_change(fc_df)
            self.set_clipped_fold_change(fc_clip_df)

            return self.fold_change, self.clipped_fold_change

        else:
            fc_all = self.log_norm[["chr", "bin"]]
            clones_norm_df = self.log_norm.drop(columns=["chr", "bin", control_name]).mean(axis=1)
            tmp_fc = pd.DataFrame({"fc": clones_norm_df - self.log_norm[control_name]})
            fc_all = pd.concat([fc_all, tmp_fc], axis=1)

            fc_clip_all = self.log_norm_clip[["chr", "bin"]]
            clones_clip_norm_df = self.log_norm_clip.drop(columns=["chr", "bin", control_name + "_cig_filt"]).mean(
                axis=1)
            tmp_clip_fc = pd.DataFrame(
                {"clip_fc": clones_clip_norm_df - self.log_norm_clip[control_name + "_cig_filt"]})
            fc_clip_all = pd.concat([fc_clip_all, tmp_clip_fc], axis=1)
            # print(fc_all)
            # print(fc_clip_all)
            self.set_fold_change(fc_all)
            self.set_clipped_fold_change(fc_clip_all)
            return self.fold_change, self.clipped_fold_change

    def plot_background(self, fig):  # , df_counts
        read_counts = self.parameters["read_counts"]
        # # print(read_counts)
        # sorted_df = read_counts.sort_values(["chr", "bin"])
        # # print(sorted_df)
        # # self.sorted_chromosomes(read_counts["chr"])
        # print(self.sorted_chromosomes(read_counts["chr"]))
        # sorted_df.to_csv("read_counts.txt", sep="\t")
        # read_counts = df_counts
        coordinates_x = []
        coordinates_y = []
        chromosomes = []
        start = 0

        for ch in read_counts["chr"]:
            # if ch[ch.find("c"):] not in chromosomes:
            if ch not in chromosomes:
                chromosomes.append(ch)
        # if the chromosomes are not sorted in ascending order, using the sort_
        # chromosomes method, we sort alphanumeric strings
        chromosomes = self.sorted_chromosomes(chromosomes)
        # print(chromosomes)
        for chrom in chromosomes:
            single_df = read_counts[read_counts["chr"] == chrom]
            # here start and end are updated every time, to allow concatenation of chromosomes on the plot
            # the average position within the interval of each chromosome take in account only the length of the
            # interval itself
            length = (single_df["bin"].iloc[-1] + 1) * self.parameters["bin_size"]
            tmp_end = start + length
            # print("start", start)
            # print("end", tmp_end)
            avg = start + length / 2
            coordinates_x.append(avg)
            coordinates_y.append(-5)
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
        # fig.show()
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

    def plot_counts_distributions(self, saving_folder):  # , cigar
        read_counts = self.parameters["read_counts"]
        counts_to_plot = []
        col_labels = []

        # --plot_raw_counts_distribution--
        fig_all = go.Figure()
        fig_all.update_yaxes(title_text="Raw_Counts")

        for col in read_counts.columns:
            if col != "chr" and col != 'bin' and "cig_filt" not in col:
                counts_to_plot.append(read_counts[col])
                col_labels.append(col)

                fig_all.add_trace(go.Violin(y=read_counts[col],
                                            box_visible=True,
                                            meanline_visible=True,
                                            name=col))

        # --plot_norm_counts_distribution--

        fig_all_norm = go.Figure()
        fig_all_norm.update_traces(opacity=0.75)
        fig_all_norm.update_yaxes(title_text="Normalized_Counts")

        for col in self.norm:
            if col != "chr" and col != "bin":
                fig_all_norm.add_trace(go.Violin(y=self.norm[col],
                                                 box_visible=True,
                                                 meanline_visible=True,
                                                 name=col))

        # --plot_log_scaled_norm_counts_distribution--

        fig_log_norm = go.Figure()
        fig_log_norm.update_traces(opacity=0.75)
        fig_log_norm.update_yaxes(title_text="Log_Normalized_Counts")

        for col in self.log_norm:
            # fig_log_norm.add_trace(go.Box(x=list(self.log_norm.rx(True, i + 1)),
            #                               name=str(self.log_norm.colnames[i])))
            if col != "chr" and col != "bin":
                fig_log_norm.add_trace(go.Violin(y=self.log_norm[col],
                                                 box_visible=True,
                                                 meanline_visible=True,
                                                 name=col))

        # --Titles--

        fig_all.update_layout(title_text="Violin Plot RAW data - BRAN" + str(self.bin_size) + " - all samples",
                              legend_orientation="h")

        fig_all_norm.update_layout(title_text="Violin Plot NORMALIZED data - BRAN" +
                                              str(self.bin_size) + " - all samples",
                                   legend_orientation="h")

        fig_log_norm.update_layout(title_text="Violin Plot LOG-SCALE NORMALIZED data - BRAN" +
                                              str(self.bin_size) + " - all samples",
                                   legend_orientation="h")

        # fig_all.show()
        # fig_all_norm.show()
        # fig_log_norm.show()
        # save_fig_all = fig_all.write_image(saving_folder + "all_sample_dist_" + str(self.bin_size) + ".jpeg",
        #                                    width=1280,
        #                                    height=1024)
        # save_fig_all_norm = fig_all_norm.write_image(saving_folder + "all_sample_norm_dist_" +
        #                                              str(self.bin_size) + ".jpeg",
        #                                              width=1280,
        #                                              height=1024)
        # save_fig_log_norm = fig_log_norm.write_image(saving_folder + "all_sample_norm_log_dist_" +
        #                                              str(self.bin_size) + ".jpeg",
        #                                              width=1280,
        #                                              height=1024)

    def plot_chrom_sample(self, saving_folder, reference, chrom, sample, template, ns=False,
                          fig=go.Figure()):  # , cigar
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

        fig.update_xaxes(title_text="Chromosomes_Position")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        c_name = None
        for c in list(read_counts["chr"]):
            if c.endswith(str(chrom)):
                c_name = c

        single_chrom = read_counts[read_counts["chr"] == c_name]
        col_list = list(single_chrom.columns)
        hover_pos = single_chrom["bin"] * self.bin_size
        for i in range(len(col_list[:col_list.index(sample)]) + 1):
            if col_list[i] == sample:
                fig.add_trace(go.Scatter(x=hover_pos,  # list(single_chrom.index * self.bin_size),
                                         y=single_chrom[sample],
                                         hovertext=hover_pos,
                                         hovertemplate=
                                         "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                         mode="markers",
                                         name=str(col_list[i])))

        if ns:
            self.add_ns_trace(fig, reference=reference, chrom=chrom)

        fig.update_layout(title="Read Counts - Clone: " +
                                sample +
                                " - Chr: " + str(chrom) +
                                " - Bin Size: " + str(self.bin_size),
                          template=template,
                          legend_orientation="h")

        # fig.show()
        # save_fig = fig.write_image(saving_folder +
        #                            "scatter_counts_chr" +
        #                            str(chrom) +
        #                            sample +
        #                            str(self.bin_size) +
        #                            ".jpeg",
        #                            width=1280,
        #                            height=1024)

    def plot_chromosome(self, saving_folder, reference, chrom, template, ns=False, fig=go.Figure()):  # , cigar
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

        fig.update_xaxes(title_text="Chromosomes_Position")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        c_name = None
        for c in list(read_counts["chr"]):
            if c.endswith(str(chrom)):
                c_name = c

        single_chrom = read_counts[read_counts["chr"] == c_name]
        col_list = list(single_chrom.columns)
        hover_pos = single_chrom["bin"] * self.bin_size
        for i in range(len(col_list)):
            if col_list[i] != "chr" and col_list[i] != "bin" and "cig_filt" not in col_list[i]:
                fig.add_trace(go.Scatter(x=hover_pos,  # list(single_chrom.index * self.bin_size),
                                         y=single_chrom[col_list[i]],
                                         hovertext=hover_pos,
                                         hovertemplate=
                                         "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                         mode="markers",
                                         name=str(col_list[i])))
        if ns:
            self.add_ns_trace(fig, reference=reference, chrom=chrom)

        fig.update_layout(title="Read Counts - All Clones - Chr: " + str(chrom) +
                                " - Bin Size: " + str(self.bin_size),
                          template=template,
                          legend_orientation="h")

        # fig.show()
        # save_fig = fig.write_image(saving_folder +
        #                            "scatter_counts_chr_" +
        #                            str(chrom) + "_" +
        #                            str(self.bin_size) +
        #                            ".jpeg",
        #                            width=1280,
        #                            height=1024)

    def plot_sample(self, saving_folder, reference, sample, template, ns=False, fig=go.Figure()):  # , cigar
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

        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        col_read_counts = list(read_counts.columns)
        hover_pos = read_counts["bin"] * self.bin_size
        for i in range(len(col_read_counts[:col_read_counts.index(sample)]) + 1):
            if col_read_counts[i] == sample:
                fig.add_trace(go.Scatter(x=list(read_counts.index * self.bin_size),
                                         y=read_counts[sample],
                                         hovertext=hover_pos,
                                         hovertemplate=
                                         "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                         mode="markers",
                                         name=str(col_read_counts[i])))

        if ns:
            self.add_ns_trace(fig, reference=reference)

        self.plot_background(fig)

        fig.update_layout(title="Read Counts - Clone: " + sample +
                                " - Bin Size: " + str(self.bin_size),
                          template=template,
                          legend_orientation="h")

        # fig.show()
        # save_fig = fig.write_image(saving_folder +
        #                            "scatter_counts_" +
        #                            sample + "_" +
        #                            str(self.bin_size) + ".jpeg",
        #                            width=1280,
        #                            height=1024)

    def plot_all(self, saving_folder, reference, template, ns=False, fig=go.Figure()):  # df_counts, bin_size / , cigar
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
        col_list = list(read_counts.columns)
        hover_pos = read_counts["bin"] * self.bin_size

        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Raw_Read_Count_Per_Bin")

        for i in range(len(col_list)):
            if col_list[i] != "chr" and col_list[i] != "bin" and "cig_filt" not in col_list[i]:
                fig.add_trace(go.Scatter(x=list(read_counts.index * self.bin_size),
                                         y=read_counts[col_list[i]],
                                         hovertext=hover_pos,
                                         hovertemplate=
                                         "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                         mode="markers",
                                         # name=str(col_list[i][:col_list[i].find("_Illumina")])))
                                         name=str(col_list[i][:col_list[i].find("_")])))
                # general_application------------------------------------------------------------------
                fig.update_layout(title="Read Counts - All Clones - All Chromosomes - Bin Size: " +
                                        str(self.bin_size),
                                  template=template,
                                  legend_orientation="h")

                # test_application----------------------------------------------------------------------
                # if "ref" not in col_list[i]:
                #     fig.update_layout(title="{}<br>Read Counts - All Chromosomes - Bin Size: {}".format(col_list[i][:col_list[i].find("_")],
                #                                                                                 str(self.bin_size
                #                                                                                     )),
                #                       template=template,
                #                       legend_orientation="h")
        if ns:
            self.add_ns_trace(fig, reference=reference)

        self.plot_background(fig)

        fig.show()

        save_fig = fig.write_image(saving_folder +
                                   "scatter_all_counts_"
                                   + str(self.bin_size) +
                                   ".jpeg",
                                   width=1920,
                                   height=1080)

    def plot_norm_data_chr_sample(self, saving_folder, reference, chrom, sample, template, ns=False,
                                  fig=go.Figure()):  # , cigar
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
        fig.update_xaxes(title_text="Chromosome_Position")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

        c_name = None
        for c in self.norm["chr"].value_counts().index:
            if c.endswith(str(chrom)):
                c_name = c

        single_chrom = self.norm[self.norm["chr"] == c_name]
        hover_pos = single_chrom["bin"] * self.bin_size
        for col in single_chrom:
            if col == sample:
                fig.add_trace(go.Scatter(x=hover_pos,  # list(single_chrom.index * self.bin_size),
                                         y=single_chrom[col],
                                         hovertext=hover_pos,
                                         hovertemplate=
                                         "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y:.0f}",
                                         mode="markers",
                                         name=col))
        if ns:
            self.add_ns_trace(fig, reference=reference, chrom=chrom)

        fig.update_layout(title="Normalized Read Counts - Clone: " +
                                sample +
                                " - Chr: " +
                                str(chrom) +
                                " - Bin Size: " +
                                str(self.bin_size),
                          template=template,
                          legend_orientation="h")

        # fig.show()

        # save_fig = fig.write_image(saving_folder + "scatter_norm_counts_chr_" +
        #                            str(chrom) + "_" +
        #                            sample + "_" +
        #                            str(self.bin_size) + ".jpeg",
        #                            width=1280,
        #                            height=1024)

    def plot_norm_data_chr(self, saving_folder, reference, chrom, template, ns=False, fig=go.Figure()):  # , cigar
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
        fig.update_xaxes(title_text="Chromosome_position")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

        c_name = None
        for c in self.norm["chr"].value_counts().index:
            if c.endswith(str(chrom)):
                c_name = c

        single_chrom = self.norm[self.norm["chr"] == c_name]
        hover_pos = single_chrom["bin"] * self.bin_size
        for col in single_chrom:
            if col != "chr" and col != "bin":
                fig.add_trace(go.Scatter(x=hover_pos,  # list(single_chrom.index * self.bin_size),
                                         y=single_chrom[col],
                                         hovertext=hover_pos,
                                         hovertemplate=
                                         "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y:.0f}",
                                         mode="markers",
                                         name=col))
        if ns:
            self.add_ns_trace(fig, reference=reference, chrom=chrom)

        fig.update_layout(title="Normalized Read Counts - Clone: all - Chr: " +
                                str(chrom) +
                                " - Bin Size: " +
                                str(self.bin_size),
                          template=template,
                          legend_orientation="h")

        # fig.show()
        #
        # save_fig = fig.write_image(saving_folder +
        #                            "scatter_norm_counts_chr_" +
        #                            str(chrom) + "_" +
        #                            str(self.bin_size) +
        #                            ".jpeg",
        #                            width=1280,
        #                            height=1024)

    def plot_norm_data_sample(self, saving_folder, reference, sample, template, control_name, ns=False,
                              fig=go.Figure()):  # cigar,
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
        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

        hover_pos = self.norm["bin"] * self.bin_size
        for col in self.norm:
            if col == sample:
                fig.add_trace(go.Scatter(x=list(self.norm.index * self.bin_size),
                                         y=self.norm[col],
                                         hovertext=hover_pos,
                                         hovertemplate=
                                         "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y:.0f}",
                                         mode="markers",
                                         name=col))

        if ns:
            self.add_ns_trace(fig, reference=reference)

        self.plot_background(fig)

        fig.update_layout(title="Normalized Read Counts - Clone: " +
                                sample +
                                " - Chr: all - Bin Size: " +
                                str(self.bin_size),
                          template=template,
                          legend_orientation="h")

        # fig.show()
        # save_fig = fig.write_image(saving_folder +
        #                            "scatter_norm_counts_" +
        #                            sample + "_" +
        #                            str(self.bin_size) +
        #                            ".jpeg",
        #                            width=1280,
        #                            height=1024)

    def plot_norm_data_all(self, saving_folder, reference, template, ns=False, fig=go.Figure()):  # df_counts, bin_size / cigar
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
        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

        hover_pos = self.norm["bin"] * self.bin_size
        for col in self.norm:
            if col != "chr" and col != "bin":
                fig.add_trace(go.Scatter(x=list(self.norm.index * self.bin_size),
                                         y=self.norm[col],
                                         hovertext=hover_pos,
                                         hovertemplate=
                                         "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y:.0f}",
                                         mode="markers",
                                         # name=col[:col.find("_Illumina")]))
                                         name=col[:col.find("_")]))
                # general_application------------------------------------------------------------------
                fig.update_layout(title="Normalized Read Counts - Clone: all - Chr: all - Bin Size: " +
                                        str(self.bin_size),
                                  template=template,
                                  legend_orientation="h")

                # test_application-----------------------------------------------------------------------
                # if "ref" not in col:
                #     fig.update_layout(title="{}<br>Normalized Read Counts - All Chromosomes - Bin Size: {}".format(col[:col.find("_")],
                #                                                                                              str(self.bin_size
                #                                                                                                  )),
                #                   template=template,
                #                   legend_orientation="h")

        if ns:
            self.add_ns_trace(fig, reference=reference)

        self.plot_background(fig)

        fig.show()
        save_fig = fig.write_image(saving_folder + "scatter_norm_all_counts_" +
                                   str(self.bin_size) +
                                   ".jpeg",
                                   width=1920,
                                   height=1080)

    def plot_fold_change_chr_sample(self, pairwise, fc, chrom, sample, control_name, saving_folder):
        """"""
        if pairwise:
            fig = go.Figure()

            fig.update_xaxes(title_text="Chromosome_Position")
            fig.update_yaxes(title_text="Fold-Change")

            c_name = None
            for c in self.fold_change["chr"].value_counts().index:
                if c.endswith(str(chrom)):
                    c_name = c

            single_chrom = self.fold_change[self.fold_change["chr"] == c_name]
            for col in single_chrom:
                if col == sample + "-" + control_name:
                    sig_data_pos = single_chrom[["chr", "bin", col]][single_chrom[col] > fc]
                    sig_data_neg = single_chrom[["chr", "bin", col]][single_chrom[col] < -fc]
                    sig_data = pd.concat([sig_data_pos, sig_data_neg])
                    not_sig_data = single_chrom[["chr", "bin", col]].drop(
                        list(sig_data_pos.index) + list(sig_data_neg.index))

                    hover_pos_sig = sig_data["bin"] * self.bin_size
                    hover_pos_no_sig = not_sig_data["bin"] * self.bin_size

                    fig.add_trace(go.Scatter(x=hover_pos_sig,  # list(sig_data.index * self.bin_size),
                                             y=sig_data[col],
                                             mode="markers",
                                             hovertext=hover_pos_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             name=col))
                    fig.add_trace(go.Scatter(x=hover_pos_no_sig,  # list(not_sig_data.index * self.bin_size),
                                             y=not_sig_data[col],
                                             mode="markers",
                                             marker=dict(size=5, color="rgb(176, 196, 222)"),
                                             hovertext=hover_pos_no_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             legendgroup="group",
                                             name=col))

                    fig.update_layout(title="Pairwise Fold Change - Chromosome: " + c_name +
                                            " - Clone: " + str(sample) + " vs " + str(control_name) +
                                            " - Bin Size: " + str(self.bin_size) +
                                            " - Threshold_FC: " + str(fc),
                                      legend_orientation="h")

                    # save_fig = fig.write_image(saving_folder +
                    #                            "pairwise_fold_change_" +
                    #                            c_name + "_" +
                    #                            sample + "_" +
                    #                            str(self.bin_size) +
                    #                            ".jpeg",
                    #                            width=1280,
                    #                            height=1024)

            # fig.show()

        else:
            print("""ATTENTION: if parameter '-pw' not give, its impossible to retrieve graphical information 
                  on single sample fold-change. \nPlease TRY AGAIN specifying '-pw' or '--pairwise' in command line""")

    def plot_fold_change_chr(self, pairwise, fc, chrom, saving_folder):
        """"""
        fig = go.Figure()

        fig.update_xaxes(title_text="Chromosome_Position")
        fig.update_yaxes(title_text="Fold-Change")

        c_name = None
        for c in self.fold_change["chr"].value_counts().index:
            if c.endswith(str(chrom)):
                c_name = c
        single_chrom = self.fold_change[self.fold_change["chr"] == c_name]
        if pairwise:
            for col in single_chrom:
                if col != "chr" and col != "bin":
                    sig_data_pos = single_chrom[["chr", "bin", col]][single_chrom[col] > fc]
                    sig_data_neg = single_chrom[["chr", "bin", col]][single_chrom[col] < -fc]
                    sig_data = pd.concat([sig_data_pos, sig_data_neg])
                    not_sig_data = single_chrom[["chr", "bin", col]].drop(
                        list(sig_data_pos.index) + list(sig_data_neg.index))

                    hover_pos_sig = sig_data["bin"] * self.bin_size
                    hover_pos_no_sig = not_sig_data["bin"] * self.bin_size

                    fig.add_trace(go.Scatter(x=hover_pos_sig,  # list(sig_data.index * self.bin_size),
                                             y=sig_data[col],
                                             mode="markers",
                                             hovertext=hover_pos_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             name=col))
                    fig.add_trace(go.Scatter(x=hover_pos_no_sig,  # list(not_sig_data.index * self.bin_size),
                                             y=not_sig_data[col],
                                             mode="markers",
                                             marker=dict(size=5, color="rgb(176, 196, 222)"),
                                             hovertext=hover_pos_no_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             legendgroup="group",
                                             name=col))

                    fig.update_layout(title="Pairwise Fold Change - Chromosome: " + c_name +
                                            " - Each clone vs reference - Bin Size: " + str(self.bin_size) +
                                            " - Threshold_FC: " + str(fc),
                                      legend_orientation="h")

                    # save_fig = fig.write_image(saving_folder +
                    #                            "pairwise_fold_change_" +
                    #                            c_name + "_" +
                    #                            str(self.bin_size) +
                    #                            ".jpeg",
                    #                            width=1280,
                    #                            height=1024)

            # fig.show()

        else:
            for col in list(single_chrom.columns):
                if col != "bin" and col != "chr":
                    sig_data_pos = single_chrom[single_chrom[col] > fc]
                    sig_data_neg = single_chrom[single_chrom[col] < -fc]
                    sig_data = pd.concat([sig_data_pos, sig_data_neg])
                    not_sig_data = single_chrom.drop(list(sig_data_pos.index) + list(sig_data_neg.index))
                    hover_pos_sig = sig_data["bin"] * self.bin_size
                    hover_pos_no_sig = not_sig_data["bin"] * self.bin_size

                    fig.add_trace(go.Scatter(x=hover_pos_sig,  # list(sig_data.index * self.bin_size),
                                             y=sig_data[col],
                                             mode="markers",
                                             hovertext=hover_pos_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             name="Significant Differences"))
                    fig.add_trace(go.Scatter(x=hover_pos_no_sig,  # list(not_sig_data.index * self.bin_size),
                                             y=not_sig_data[col],
                                             mode="markers",
                                             marker=dict(size=5, color="rgb(176, 196, 222)"),
                                             hovertext=hover_pos_no_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             legendgroup="group",
                                             name="Not Significant Differences"))

                    fig.update_layout(title="Fold Change - Chromosome: " + c_name + "Clone: all vs Reference - "
                                                                                    "Bin Size: " + str(
                        self.bin_size) + " - Threshold_FC: " + str(fc),
                                      legend_orientation="h")

                # save_fig = fig.write_image(saving_folder +
                #                            "fold_change_" +
                #                            c_name + "_" +
                #                            str(self.bin_size) +
                #                            ".jpeg",
                #                            width=1280,
                #                            height=1024)

            # fig.show()

    def plot_fold_change_sample(self, pairwise, fc, sample, control_name, saving_folder):
        """"""
        if pairwise:
            fig = go.Figure()

            fig.update_xaxes(title_text="Genome_Position")
            fig.update_yaxes(title_text="Fold_change")

            for col in self.fold_change:
                if col == sample + "-" + control_name:
                    sig_data_pos = self.fold_change[["chr", "bin", col]][self.fold_change[col] > fc]
                    sig_data_neg = self.fold_change[["chr", "bin", col]][self.fold_change[col] < -fc]
                    sig_data = pd.concat([sig_data_pos, sig_data_neg])
                    not_sig_data = self.fold_change[["chr", "bin", col]].drop(list(sig_data_pos.index) +
                                                                              list(sig_data_neg.index))

                    hover_pos_sig = sig_data["bin"] * self.bin_size
                    hover_pos_no_sig = not_sig_data["bin"] * self.bin_size

                    fig.add_trace(go.Scatter(x=list(sig_data.index * self.bin_size),
                                             y=sig_data[col],
                                             mode="markers",
                                             hovertext=hover_pos_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             name=col))
                    fig.add_trace(go.Scatter(x=list(not_sig_data.index * self.bin_size),
                                             y=not_sig_data[col],
                                             mode="markers",
                                             marker=dict(size=5, color="rgb(176, 196, 222)"),
                                             hovertext=hover_pos_no_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             legendgroup="group",
                                             name=col))

                    fig.update_layout(title="Pairwise Fold Change - All chromosomes - Clone: " +
                                            str(sample) + " vs " + str(control_name) +
                                            " - Bin Size: " + str(self.bin_size) +
                                            " - Threshold_FC: " + str(fc),
                                      legend_orientation="h")

                    # save_fig = fig.write_image(saving_folder +
                    #                            "pairwise_fold_change_sample_"
                    #                            + str(self.bin_size) +
                    #                            ".jpeg",
                    #                            width=1280,
                    #                            height=1024)

            self.add_threshold_fc(fig, fc)
            self.plot_background(fig)
            # fig.show()

        else:
            print("""ATTENTION: if parameter '-pw' not give, its impossible to retrieve graphical information 
                  on single sample fold-change. \nPlease TRY AGAIN specifying '-pw' or '--pairwise' in command line""")

    def plot_fold_change(self, fc, saving_folder, pairwise, control_name):
        """"""
        fig = go.Figure()

        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="log2_Fold-Change")
        summary_sig_data = {"chr": [], "start_pos": [], "end_pos": [], "clone_name": [], "type": [], "fc": []}

        if pairwise:
            for col in self.fold_change:
                if col != "bin" and col != "chr":
                    sig_data_pos = self.fold_change[["chr", "bin", col]][self.fold_change[col] > fc]
                    sig_data_neg = self.fold_change[["chr", "bin", col]][self.fold_change[col] < -fc]
                    sig_data = pd.concat([sig_data_pos, sig_data_neg])

                    # print(sig_data)
                    summary_sig_data["chr"] += list(sig_data["chr"])
                    # print(len(list(sig_data["chr"])))
                    summary_sig_data["start_pos"] += (list(sig_data["bin"] * self.bin_size))
                    # print(len(sig_data["bin"] * self.bin_size))
                    summary_sig_data["end_pos"] += (list((sig_data["bin"] * self.bin_size) + self.bin_size))
                    # print(len((sig_data["bin"] * self.bin_size) + self.bin_size))
                    summary_sig_data["clone_name"] += [col] * len(sig_data)
                    # print(len([col] * len(sig_data)))
                    summary_sig_data["type"] += ["read_count"] * len(sig_data)
                    # print(len(["clipped_count"] * len(sig_data)))
                    summary_sig_data["fc"] += ["+"] * len(sig_data_pos) + ["-"] * len(sig_data_neg)
                    # print(["+"] * len(sig_data_pos))
                    not_sig_data = self.fold_change[["chr", "bin", col]].drop(list(sig_data_pos.index) +
                                                                              list(sig_data_neg.index))

                    hover_pos_sig = sig_data["bin"] * self.bin_size
                    hover_pos_no_sig = not_sig_data["bin"] * self.bin_size

                    fig.add_trace(go.Scatter(x=list(sig_data.index * self.bin_size),
                                             y=sig_data[col],
                                             mode="markers",
                                             hovertext=hover_pos_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             # name=col[:col.find("_Illumina")],
                                             legendgroup= "group",
                                             name=col[:col.find("_")]))
                    fig.add_trace(go.Scatter(x=list(not_sig_data.index * self.bin_size),
                                             y=not_sig_data[col],
                                             mode="markers",
                                             marker=dict(size=5, color="rgb(176, 196, 222)"),
                                             hovertext=hover_pos_no_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             legendgroup="group",
                                             # name=col[:col.find("_Illumina")] + "_no_significant"))
                                             name=col[:col.find("_")] + "_no_significant"))

                    # general application ----------------------------------------------------------------------
                    fig.update_layout(
                            title="Each <i>vs</i> {}<br>Pairwise log2 Fold Change - Bin Size: {} - Threshold_FC: ".format(control_name[:control_name.find("_Illumina")],
                                                                                                               str(self.bin_size),
                                                                                                               str(fc)),
                            legend_orientation="h")

                    # test_application ---------------------------------------------------------------------------------
                    # fig.update_layout(
                    #     title="{} <i>vs</i> {}<br>Pairwise log2 Fold Change - Bin Size: {} - Threshold_FC: ".format(col[:col.find("_")],
                    #                                                                                        control_name[:control_name.find("_")],
                    #                                                                                        str(self.bin_size),
                    #                                                                                        str(fc)),
                    #     legend_orientation="h")

            self.add_threshold_fc(fig, fc)
            self.plot_background(fig)

            fig.show()
            save_fig = fig.write_image(saving_folder +
                                       "pairwise_fold_change_" +
                                       str(self.bin_size) +
                                       ".jpeg",
                                       width=1920,
                                       height=1080)

            summary_sig_data = pd.DataFrame(summary_sig_data)
            # print(summary_sig_data)
            self.set_sig_data(summary_sig_data)
            return self.sig_data

        else:
            for col in list(self.fold_change.columns):
                if col != "bin" and col != "chr":
                    sig_data_pos = self.fold_change[self.fold_change[col] > fc]
                    sig_data_neg = self.fold_change[self.fold_change[col] < -fc]
                    sig_data = pd.concat([sig_data_pos, sig_data_neg])

                    # print(sig_data)
                    summary_sig_data["chr"] += list(sig_data["chr"])
                    # print(len(list(sig_data["chr"])))
                    summary_sig_data["start_pos"] += (list(sig_data["bin"] * self.bin_size))
                    # print(len(sig_data["bin"] * self.bin_size))
                    summary_sig_data["end_pos"] += (list((sig_data["bin"] * self.bin_size) + self.bin_size))
                    # print(len((sig_data["bin"] * self.bin_size) + self.bin_size))
                    summary_sig_data["clone_name"] += [col] * len(sig_data)
                    # print(len([col] * len(sig_data)))
                    summary_sig_data["event"] += ["read_count"] * len(sig_data)
                    # print(len(["clipped_count"] * len(sig_data)))

                    not_sig_data = self.fold_change.drop(list(sig_data_pos.index) + list(sig_data_neg.index))
                    hover_pos_sig = sig_data["bin"] * self.bin_size
                    hover_pos_no_sig = not_sig_data["bin"] * self.bin_size

                    fig.add_trace(go.Scatter(x=list(sig_data.index * self.bin_size),
                                             y=sig_data[col],
                                             mode="markers",
                                             hovertext=hover_pos_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             name="Significant Differences"))
                    fig.add_trace(go.Scatter(x=list(not_sig_data.index * self.bin_size),
                                             y=not_sig_data[col],
                                             mode="markers",
                                             marker=dict(size=5, color="rgb(176, 196, 222)"),
                                             hovertext=hover_pos_no_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             legendgroup="group",
                                             name="Not Significant Differences"))

                    fig.update_layout(title="log2 Fold Change - Clone: all vs Reference - "
                                            "Bin Size: " + str(self.bin_size) + " - Threshold_FC: " + str(fc),
                                      legend_orientation="h")

            self.add_threshold_fc(fig, fc)
            self.plot_background(fig)

            fig.show()

            save_fig = fig.write_image(saving_folder +
                                       "fold_change_" +
                                       str(self.bin_size) +
                                       ".jpeg",
                                       width=1920,
                                       height=1080)

            summary_sig_data = pd.DataFrame(summary_sig_data)
            # print(summary_sig_data)
            self.set_sig_data(summary_sig_data)
            return self.sig_data

    def plot_filtered_reads(self, saving_folder):
        read_counts = self.parameters["read_counts"]
        fig = go.Figure()

        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Norm_Clipped_Read_counts")

        for col in read_counts:
            if "cig_filt" in col:  # and col[:col.find("cig_filt")] != control_ref:
                hover_pos = read_counts["bin"] * self.bin_size
                fig.add_trace(go.Scatter(x=list(self.norm_clip[col].index * self.bin_size),
                                         y=self.norm_clip[col],
                                         hovertext=hover_pos,
                                         hovertemplate=
                                         "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                         hoverinfo="text",
                                         mode="markers",
                                         # name=col[:col.find("_Illumina")]))
                                        name=col[:col.find("_")]))

                # general application--------------------------------------------------------------------------------
                fig.update_layout(title="Norm_Clipped_Read_Counts - Bin Size: {}".format(str(self.bin_size)),
                                  legend_orientation="h")
                # test_application-------------------------------------------------------------------------------------
                # if "ref" not in col:
                #     fig.update_layout(title="{}<br>Norm_Clipped_Read_Counts - Bin Size: {}".format(col[:col.find("_")],
                #                                                                              str(self.bin_size)),
                #                   legend_orientation="h")

        self.plot_background(fig)

        fig.show()
        save_fig = fig.write_image(saving_folder + "clipped_reads_counts" + str(self.bin_size) + ".jpeg",
                                   width=1920,
                                   height=1080)

    def plot_filt_reads_fold_change(self, fc, pairwise, control_name, saving_folder):
        fig = go.Figure()

        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Clipped_log2_Fold-Change")

        # print(self.fold_change[self.fold_change > 1.5])

        summary_sig_data = {"chr": [], "start_pos": [], "end_pos": [], "clone_name": [], "type": [], "fc": []}

        if pairwise:
            for col in self.clipped_fold_change:
                if col != "bin" and col != "chr":
                    sig_data_pos = self.clipped_fold_change[["chr", "bin", col]][self.clipped_fold_change[col] > fc]
                    sig_data_neg = self.clipped_fold_change[["chr", "bin", col]][self.clipped_fold_change[col] < -fc]
                    sig_data = pd.concat([sig_data_pos, sig_data_neg])

                    # print(sig_data)
                    summary_sig_data["chr"] += list(sig_data["chr"])
                    # print(len(list(sig_data["chr"])))
                    summary_sig_data["start_pos"] += (list(sig_data["bin"] * self.bin_size))
                    # print(len(sig_data["bin"] * self.bin_size))
                    summary_sig_data["end_pos"] += (list((sig_data["bin"] * self.bin_size) + self.bin_size))
                    # print(len((sig_data["bin"] * self.bin_size) + self.bin_size))
                    summary_sig_data["clone_name"] += [col.replace("_cig_filt", "")] * len(sig_data)
                    # print(len([col] * len(sig_data)))
                    summary_sig_data["type"] += ["clipped_count"] * len(sig_data)
                    # print(len(["clipped_count"] * len(sig_data)))
                    summary_sig_data["fc"] += ["+"] * len(sig_data_pos) + ["-"] * len(sig_data_neg)

                    not_sig_data = self.clipped_fold_change[["chr", "bin", col]].drop(list(sig_data_pos.index) +
                                                                                      list(sig_data_neg.index))

                    hover_pos_sig = sig_data["bin"] * self.bin_size
                    hover_pos_no_sig = not_sig_data["bin"] * self.bin_size

                    fig.add_trace(go.Scatter(x=list(sig_data.index * self.bin_size),
                                             y=sig_data[col],
                                             mode="markers",
                                             hovertext=hover_pos_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             legendgroup="group",
                                             # name=col[:col.find("_Illumina")]))
                                             name=col[:col.find("_")]))
                    fig.add_trace(go.Scatter(x=list(not_sig_data.index * self.bin_size),
                                             y=not_sig_data[col],
                                             mode="markers",
                                             marker=dict(size=5, color="rgb(176, 196, 222)"),
                                             hovertext=hover_pos_no_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             legendgroup="group",
                                             # name=col[:col.find("_Illumina")] + "_no_significant"))
                                             name=col[:col.find("_")] + "_no_significant"))
                    # all_sig_data_name = pd.concat([all_sig_data, name], axis=1)
                    # print(list(name.index))
                    # general application -----------------------------------------------------------------------
                    fig.update_layout(
                        title="Each <i>vs</i> {}<br>Clipped Reads Pairwise log2 Fold Change - Bin Size: {} - "
                              "Threshold_FC: {}".format(control_name[:control_name.find("_Illumina")],
                                                        str(self.bin_size),
                                                        str(fc)),
                        legend_orientation="h")
                    # test_application -----------------------------------------------------------------------------
                    # fig.update_layout(title="{} <i>vs</i> {}<br>Clipped Reads Pairwise log2 Fold Change - Bin Size: {} - "
                    #                         "Threshold_FC: {}".format(col[:col.find("_")],
                    #                                                   control_name[:control_name.find("_")],
                    #                                                   str(self.bin_size),
                    #                                                   str(fc)),
                    #                   legend_orientation="h")

            self.add_threshold_fc(fig, fc)
            self.plot_background(fig)
            fig.show()
            save_fig = fig.write_image(saving_folder +
                                       "pairwise_clipped_fold_change_" +
                                       str(self.bin_size) +
                                       ".jpeg",
                                       width=1920,
                                       height=1080)

            save_fig = fig.write_image(saving_folder +
                                       "pairwise_clipped_fold_change_" +
                                       str(self.bin_size) +
                                       ".svg",
                                       width=1920,
                                       height=1080)


            # print(len(summary_sig_data["chr"]))
            # print(len(summary_sig_data["start_pos"]))
            # print(len(summary_sig_data["end_pos"]))
            # print(len(summary_sig_data["clone_name"]))
            # print(len(summary_sig_data["event"]))
            print(summary_sig_data["fc"])

            summary_sig_data = pd.DataFrame(summary_sig_data)
            print(summary_sig_data)
            self.set_clip_sig_data(summary_sig_data)
            return self.clip_sig_data

        else:
            for col in list(self.clipped_fold_change.columns):
                print(col)
                if col != "bin" and col != "chr":
                    print(col)
                    sig_data_pos = self.clipped_fold_change[self.clipped_fold_change[col] > fc]
                    sig_data_neg = self.clipped_fold_change[self.clipped_fold_change[col] < -fc]
                    sig_data = pd.concat([sig_data_pos, sig_data_neg])

                    # print(sig_data)
                    summary_sig_data["chr"] += list(sig_data["chr"])
                    # print(len(list(sig_data["chr"])))
                    summary_sig_data["start_pos"] += (list(sig_data["bin"] * self.bin_size))
                    # print(len(sig_data["bin"] * self.bin_size))
                    summary_sig_data["end_pos"] += (list((sig_data["bin"] * self.bin_size) + self.bin_size))
                    # print(len((sig_data["bin"] * self.bin_size) + self.bin_size))
                    summary_sig_data["clone_name"] += [col.replace("_cig_filt", "")] * len(sig_data)
                    # print(len([col] * len(sig_data)))
                    summary_sig_data["event"] += ["clipped_count"] * len(sig_data)
                    # print(len(["clipped_count"] * len(sig_data)))

                    not_sig_data = self.clipped_fold_change.drop(list(sig_data_pos.index) + list(sig_data_neg.index))
                    # print(sig_data)
                    # print(self.clipped_fold_change)
                    hover_pos_sig = sig_data["bin"] * self.bin_size
                    hover_pos_no_sig = not_sig_data["bin"] * self.bin_size

                    fig.add_trace(go.Scatter(x=list(sig_data.index * self.bin_size),
                                             y=sig_data[col],
                                             mode="markers",
                                             hovertext=hover_pos_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             name="Significant Differences"))
                    fig.add_trace(go.Scatter(x=list(not_sig_data.index * self.bin_size),
                                             y=not_sig_data[col],
                                             mode="markers",
                                             marker=dict(size=5, color="rgb(176, 196, 222)"),
                                             hovertext=hover_pos_no_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             legendgroup="group",
                                             name="Not Significant Differences"))

                    fig.update_layout(title="Clipped Reads log2 Fold Change - Clone: all vs Reference - "
                                            "Bin Size: " + str(self.bin_size) + " - Threshold_FC: " + str(fc),
                                      legend_orientation="h")

            self.add_threshold_fc(fig, fc)
            self.plot_background(fig)
            fig.show()
            save_fig = fig.write_image(saving_folder +
                                       "clipped_fold_change_" +
                                       str(self.bin_size) +
                                       ".jpeg",
                                       width=1280,
                                       height=1024)

            summary_sig_data = pd.DataFrame(summary_sig_data)
            # print(summary_sig_data)
            self.set_clip_sig_data(summary_sig_data)
            return self.clip_sig_data

    def sig_positions_output_file(self, fc, control_name, file_output_path):

        sig_df = pd.concat([self.sig_data, self.clip_sig_data])
        sig_df = sig_df.sort_values(by=["clone_name", "chr", "start_pos"])

        header = "# Chromosome positions in which fold-change is > " + str(fc) + \
                 " or < -" + str(fc) + " with respect to : " + control_name + "\n"

        with open(file_output_path + "BRAN" + str(self.bin_size) + "_significant_changes_regions.tsv", "w") as file:
            file.write(header)
            sig_df.to_csv(path_or_buf=file, sep="\t", index=False)

        print(sig_df)

    def plot_violin_dist_counts(self, saving_folder):
        """"""
        # only for the scope of analysis on chardonnay
        # [ corn_flower_blue, dark_cyan, dark_orange ,crimson]
        possible_color_violin = ["rgb(100,149,237)", "rgb(0,139,139)", "rgb(255,140,0)", "rgb(220,20,60)"]
        read_counts = self.parameters["read_counts"]
        fig = make_subplots(rows=1, cols=2, subplot_titles=("Normalized Counts",
                                                            "Soft_Hard Clipped Read Counts"))
        i = 0
        j = 0
        hover_pos = self.norm["bin"] * self.bin_size
        for col in read_counts:
            if "cig_filt" in col:
                fig.add_trace(go.Violin(y=self.norm_clip[col],
                                        box_visible=True,
                                        meanline_visible=True,
                                        hovertext=hover_pos,
                                        text=self.norm["chr"],
                                        hovertemplate=
                                        "<b>Chrom</b>: %{text}" +
                                        "<br><b>Position</b>: %{hovertext:,}" +
                                        "<br>Count: %{y:.0f}",
                                        # fillcolor="rgb(0,139,139)",
                                        line_color=possible_color_violin[i],
                                        # name=col[:col.find("_Illumina")]),
                                        name=col[:col.find("_")]),
                              row=1,
                              col=2)
                fig.update_yaxes(title_text="Norm_Clipped_counts", row=1, col=2)
                fig.update_xaxes(tickangle=45)
                i += 1

            elif "cig_filt" not in col and col != "chr" and col != "bin":
                fig.add_trace(go.Violin(y=self.norm[col],
                                        box_visible=True,
                                        meanline_visible=True,
                                        hovertext=hover_pos,
                                        text=self.norm["chr"],
                                        hovertemplate=
                                        "<b>Chrom</b>: %{text}" +
                                        "<br><b>Position</b>: %{hovertext:,}" +
                                        "<br>Count: %{y:.0f}",
                                        # fillcolor="rgb(255,160,122)",
                                        line_color=possible_color_violin[j],
                                        # name=col[:col.find("_Illumina")]),
                                        name=col[:col.find("_")]),
                              row=1,
                              col=1)
                fig.update_yaxes(title_text="Norm_Read_counts", row=1, col=1)

                j += 1

        fig.update_layout(showlegend=False,
                          title="Comparison Between Read Counts and Only Clipped Read Counts per Sample" +
                                "- Bin Size: " + str(self.bin_size))

        fig.update_traces(opacity=0.75)
        fig.show()

        save_fig = fig.write_image(saving_folder + "comparison_norm_clipped_reads_counts" + str(self.bin_size) + ".jpeg",
                                   width=1920,
                                   height=1080)

    def plot_bar_chart(self, saving_folder, cigar, unmapped):
        read_counts = self.parameters["read_counts"]
        summary_read_counts = []
        summary_clipped_counts = []
        summary_unmapped_reads = list(self.parameters["unmapped_reads"].values())
        x_labels = list(self.parameters["unmapped_reads"].keys())
        perc_read_counts = []
        total_reads = []
        fig = go.Figure()

        for col in read_counts:
            if col != "chr" and col != "bin":
                if "cig_filt" in col:
                    summary_clipped_counts.append(sum(read_counts[col]))
                else:
                    summary_read_counts.append(sum(read_counts[col]))

        for i in range(len(summary_read_counts)):
            total_reads.append(summary_read_counts[i] + summary_clipped_counts[i] + summary_unmapped_reads[i])
            perc_read_counts.append((summary_read_counts[i] / total_reads[i]) * 100)

        fig.add_trace(go.Bar(x=x_labels,
                             y=perc_read_counts,
                             text=[str(int(el)) + "%" for el in perc_read_counts],
                             textposition="auto",
                             # marker_color="rgb(0, 139, 139)",  # dark cyan
                             name="%read_counts"))
        if unmapped:
            for i in range(len(summary_unmapped_reads)):
                summary_unmapped_reads[i] = (summary_unmapped_reads[i] / total_reads[i]) * 100

            fig.add_trace(go.Bar(x=x_labels,
                                 y=summary_unmapped_reads,
                                 text=[str(int(el)) + "%" for el in summary_unmapped_reads],
                                 textposition="auto",
                                 # marker_color="rgb(220,20,60)",  # crimson
                                 name="%unmapped_reads"))

        if cigar:
            for i in range(len(summary_clipped_counts)):
                summary_clipped_counts[i] = (summary_clipped_counts[i] / total_reads[i]) * 100

            fig.add_trace(go.Bar(x=x_labels,
                                 y=summary_clipped_counts,
                                 text=[str(int(el)) + "%" for el in summary_clipped_counts],
                                 textposition="auto",
                                 marker_color="rgb(220,20,60)",  # crimson
                                 name="%clipped_reads"))

        fig.update_layout(barmode="stack",
                          title_text="Proportion of Clipped Reads, Unmapped Reads and other Reads Counts - Bin Size: " +
                                     str(self.bin_size))
        fig.update_traces(opacity=0.6)
        fig.show()
        save_fig = fig.write_image(saving_folder + "bar_chart_proportion" + str(self.bin_size) + ".jpeg",
                                   width=1920,
                                   height=1080)


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
                        nargs="+",
                        default=["./"],
                        help="The path to the folder in which are located the files to be analyzed (.bam)")

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
                        # required=True,
                        help="""The name of the control group for the fold_change analysis, it should be the same 
                        name as the column name in the read_counts data structure, thus the name of the alignment 
                        file used as baseline without the ".bam" string""")

    parser.add_argument("-c", "--cigar",
                        action="store_true",
                        help="If specified, it allows the application of all the filters on cigar_string, per read")

    parser.add_argument("-cf", "--cigar_filter",
                        nargs="+",
                        default=["S", "H"],
                        help="""If specified, the reads mapped with soft and hard clipping (S and H) by default, are taken out 
                            form the read counts; it returns a data frame with same structure of the default one.
                            \n(Specify other filters like e.g. "I" "D")""")

    parser.add_argument("-fc", "--fold_change",
                        type=float,
                        default=1.5,
                        help="An integer to set the fold-change cut-off")

    parser.add_argument("-u", "--unmapped",
                        action="store_true",
                        help="""If specified, also a .txt file is  created, with all the unmapped reads and, as 
                        last raw, the counts for each sample""")

    parser.add_argument("-sf", "--saving_folder",
                        type=str,
                        default="./plots/",
                        help="""Path to the directory in which create the plots folder and save all images; 
                        if not specified, a directory 'plots' will be created in the current one""")

    parser.add_argument("-pw", "--pairwise",
                        action="store_true",
                        help="""If specified, the fold change is calculated pairwise between 
                        each clone and the reference""")

    parser.add_argument("-gi", "--general_info",
                        action="store_true",
                        help="""If specified, distribution_plot, norm_plot_all, (filtered_plot_all) and 
                        fold_change_plot are displayed""")

    pio.templates.default = "seaborn+none"
    template = "seaborn+none"


    args = parser.parse_args()
    dict_args = vars(parser.parse_args([]))

    if args.flag_list != dict_args["flag_list"]:
        flags = dict_args["flag_list"] + args.flag_list
    else:
        flags = args.flag_list

    if args.folder != dict_args["folder"]:
        args.folder = dict_args["folder"] + args.folder

    analyzer = TestingBinReadAnalyzer(args.folder,
                                      args.bin_size,
                                      args.reference,
                                      flags,
                                      args.cigar_filter,
                                      args.output_pickle)

    analyzer.load_data(cigar=args.cigar,
                       cigar_filter=args.cigar_filter,
                       reference=args.reference,
                       read_info=args.read_info,
                       unmapped=args.unmapped,
                       verbose=True)

    if not os.path.exists(args.saving_folder):
        os.mkdir(args.saving_folder)

    analyzer.no_repeats_df_transformation()
    analyzer.normalize_bins(args.control_name)
    # exit(1)
    analyzer.calc_fold_change(args.control_name, args.pairwise)

    if args.general_info and args.cigar:
        analyzer.plot_violin_dist_counts(args.saving_folder)
        analyzer.plot_bar_chart(args.saving_folder, args.cigar, args.unmapped)
        analyzer.plot_all(args.saving_folder, args.reference, template, args.Ns_count)
        analyzer.plot_norm_data_all(args.saving_folder, args.reference, template)
        analyzer.plot_filtered_reads(args.saving_folder)
        analyzer.plot_fold_change(args.fold_change, args.saving_folder, args.pairwise, args.control_name)
        analyzer.plot_filt_reads_fold_change(args.fold_change, args.pairwise, args.control_name, args.saving_folder)
        analyzer.sig_positions_output_file(args.fold_change, args.control_name, args.output_pickle)

    elif args.general_info:
        analyzer.plot_counts_distributions(args.saving_folder)
        analyzer.plot_bar_chart(args.saving_folder, args.cigar, args.unmapped)
        analyzer.plot_all(args.saving_folder, args.reference, template, args.Ns_count)
        analyzer.plot_norm_data_all(args.saving_folder, args.reference, template)
        analyzer.plot_fold_change(args.fold_change, args.saving_folder, args.pairwise)

    else:
        if args.chromosome and args.sample:
            analyzer.plot_chrom_sample(args.saving_folder, args.reference, args.chromosome,
                                       args.sample, template, args.Ns_count)
            analyzer.plot_norm_data_chr_sample(args.saving_folder, args.reference, args.chromosome,
                                               args.sample, template, args.Ns_count)
            analyzer.plot_fold_change_chr_sample(args.pairwise, args.fold_change, args.chromosome,
                                                 args.sample, args.control_name, args.saving_folder)

        elif args.chromosome:
            analyzer.plot_chromosome(args.saving_folder, args.reference, args.chromosome, template, args.Ns_count)
            analyzer.plot_norm_data_chr(args.saving_folder, args.reference, args.chromosome, template, args.Ns_count)
            analyzer.plot_fold_change_chr(args.pairwise, args.fold_change, args.chromosome, args.saving_folder)

        else:
            analyzer.plot_sample(args.saving_folder, args.reference, args.sample, template, args.Ns_count)
            analyzer.plot_norm_data_sample(args.saving_folder, args.reference, args.sample, template, args.Ns_count)
            analyzer.plot_fold_change_sample(args.pairwise, args.fold_change, args.sample,
                                             args.control_name, args.saving_folder)
