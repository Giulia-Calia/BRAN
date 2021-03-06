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
import progressbar
import re
import plotly.graph_objects as go
# Disable the orca response timeout.
import plotly.io._orca
import retrying

import plotly.express as px
import plotly.io as pio
import pandas as pd
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from testing_code_counter import TestingBinReadCounter
from testing_code_visualizer import TestingBinReadVisualizer
from testing_code_identifier import TestingBinReadIdentifier
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from plotly.subplots import make_subplots

# import plotly
# import time
# plotly.io.orca.ensure_server()
# time.sleep(10)
unwrapped = plotly.io._orca.request_image_with_retrying.__wrapped__
wrapped = retrying.retry(wait_random_min=1000)(unwrapped)
plotly.io._orca.request_image_with_retrying = wrapped

pandas2ri.activate()
base = importr('base')
utils = importr('utils')
stats = importr('stats')
edger = rpackages.importr('edgeR')


# --------------------------
# install tzlocal for pandas2ri
# --------------------------


def sorted_chromosomes(column):
    """Sort a given list of string having numbers on it"""
    # convert a number of type string into a number of type int in order to sort it
    convert = lambda text: int(text) if text.isdigit() else text
    # regular expression used to retrieve chromosome number of type string and pass them to the converter
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]  # any n. having digits from 0 to 9
    # as result a sorted list of chromosomes is returned applying the lambda functions above
    # (chr1, chr2, chr3, ..., chr_n)
    return sorted(column, key=alphanum_key)


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
    """

    def __init__(self, bam_folder_path, bin_size, ref_genome, flags, cigar_filter, out_pickle):
        self.folder = bam_folder_path
        self.bin_size = bin_size
        self.flags = flags
        self.ref = ref_genome
        self.cigar_filter = cigar_filter
        self.out = out_pickle
        self.parameters = None
        self.read_counts = None
        self.norm = None
        self.log_norm = None
        self.norm_clip = None
        self.log_norm_clip = None
        self.norm_unmapped = None
        self.fold_change = None
        self.clipped_fold_change = None
        self.sig_bins = None
        self.nosig_bins = None
        # self.sig_clip_bins = None
        # self.nosig_clip_bins = None

    def set_parameters(self, param):
        """set parameters as actual parameters"""
        self.parameters = param

    def set_read_counts(self, sort_read_counts):
        self.read_counts = sort_read_counts

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

    def set_sig_bins(self, sig_data):
        """set the data structure of significant read fold-change counts"""
        self.sig_bins = sig_data

    def set_no_sig_bins(self, no_sig_data):
        self.nosig_bins = pd.concat([no_sig_data])

    # def set_sig_clip_bins(self, sig_clip_data):
    #     """set the data structure of significant clipped fold-change counts"""
    #     self.sig_clip_bins = sig_clip_data
    #
    # def set_no_sig_clip_bins(self, no_sig_clip_data):
    #     self.nosig_clip_bins = pd.concat([no_sig_clip_data])

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
        counter = TestingBinReadCounter(self.folder, self.bin_size, self.flags, self.ref, self.cigar_filter, self.out)
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
                    found = True
                    if verbose:
                        print("Same parameters; import from: ", file)
                        print("\n", parameters)

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

    def sorted_df(self, saving_folder):
        """Return a sorted data-frame on a chromosome order base"""
        read_counts = self.parameters["read_counts"]
        chr_list = list(read_counts["chr"].value_counts().index)
        # use the lambda functions below to retrieve a sorted list of chromosomes from 1 to whatever
        sorted_chr_list = sorted_chromosomes(chr_list)
        sorted_df = pd.DataFrame()
        # to rebuilt the data-frame in a sorted order of chromosomes have to iterate and concatenate a single chromosome
        # df at a time
        for chrom in sorted_chr_list:
            single_df = read_counts[read_counts["chr"] == chrom]
            sorted_df = pd.concat([sorted_df, single_df])
        sorted_df = sorted_df.reset_index(drop=True)
        self.set_read_counts(sorted_df)
        # print("sorted_read_counts: \n", self.parameters["read_counts"])
        # sorted_df.to_csv(saving_folder + "sorted_df.txt", sep="\t")
        return self.read_counts

    def norm_structures(self, counts_edger):

        counts_edger_df = robjects.DataFrame(counts_edger)  # R data frame of raw counts
        norm_counts = edger.cpm(counts_edger_df, normalized_lib_sizes=True)  # R object of norm counts
        log_norm_counts = edger.cpm(counts_edger_df, log=True)  # R object of log norm counts

        # in order to pass from an R dataframe to a python dataframe in the most secure way, the R structure has to
        # be decomposed and recomposed into the python one
        norm_dict = {}
        log_norm_dict = {}
        for i in range(1, norm_counts.ncol + 1):
            norm_dict[norm_counts.colnames[i - 1]] = list(norm_counts.rx(True, i))
            log_norm_dict[log_norm_counts.colnames[i - 1]] = list(log_norm_counts.rx(True, i))

        return self.norm_dfs(norm_dict), self.norm_dfs(log_norm_dict)

    def norm_dfs(self, counts_dict):
        # pandas data frame of normalized counts
        counts_df = pd.DataFrame(counts_dict)
        counts_df = pd.concat([self.read_counts["chr"], self.read_counts["bin"], counts_df], axis=1)

        return counts_df

    def normalize_bins(self, control_name, saving_folder):
        """This method handles the normalization of the raw read_counts
        using an R package, edgeR, imported thanks to rpy2 that provides
        an already implemented function, cpm, for the normalization of a
        table of counts, as well as a series of other function specifically
        implemented for RNA-Seq; the normalization for clipped counts is done
        manually"""
        # self.read_counts.to_csv(str(self.bin_size) + "_norm_counts_check.txt", sep="\t")
        # read_counts = self.no_repeats_df_transformation()
        # the edgeR package is imported using rpy2 syntax to access to all its built-in functions
        read_counts_edger = {}  # a dictionary of sample: vector_of_counts to work with edger
        # clipped_count = {}
        clipped_counts_edger = {}
        print("\n")
        print("Normalization process:\n")

        norm_bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
        update_bar = 0

        for col in self.read_counts:
            if col != "chr" and col != "bin" and "cig_filt" not in col:
                # creates an element of the dictionary [sample_name]: r_vector_of_counts
                read_counts_edger[col] = robjects.IntVector(self.read_counts[col])
                # print(read_counts_edger[col])

            elif "cig_filt" in col:
                # print(list(read_counts[col]))
                # clipped_count[col] = self.read_counts[col]
                clipped_counts_edger[col] = robjects.IntVector(self.read_counts[col])
            else:
                continue

            update_bar += 1
            norm_bar.update(update_bar)

        norm_counts_df = self.norm_structures(read_counts_edger)[0]
        log_norm_counts_df = self.norm_structures(read_counts_edger)[1]
        self.set_norm(norm_counts_df)
        self.set_log_norm(log_norm_counts_df)

        norm_clip_df = self.norm_structures(clipped_counts_edger)[0]
        log_norm_clip_df = self.norm_structures(clipped_counts_edger)[1]
        self.set_norm_clip(norm_clip_df)
        self.set_log_norm_clip(log_norm_clip_df)

        unmapped = self.parameters["unmapped_reads"]
        norm_unmapped = {}
        for u in unmapped:
            norm_unmapped[u] = unmapped[u] / (sum(self.read_counts[u]) / 1000000)
            update_bar += 1
            norm_bar.update(update_bar)

        self.set_norm_unmapped(norm_unmapped)

        # -----------uncomment to have files for a check on the normalization results--------------------
        # norm_counts_df.to_csv(saving_folder + "norm_mod_float_counts.tsv", sep="\t")
        # norm_clip_df.to_csv(saving_folder + "norm_mod_float_norm_clipped_counts.tsv", sep="\t")
        # log_norm_counts_df.to_csv(saving_folder + "norm_mod_log_norm_counts.tsv", sep="\t")
        # log_norm_clip_df.to_csv(saving_folder + "norm_mod_log_norm_clipped_counts.tsv", sep="\t")

        # print(clipped_counts_edger)
        # read_counts_edger_df = robjects.DataFrame(read_counts_edger)  # R data frame of raw counts
        # norm_counts = edger.cpm(read_counts_edger_df, normalized_lib_sizes=True)  # R object of norm counts
        # log_norm_counts = edger.cpm(read_counts_edger_df, log=True)  # R object of log norm counts
        #
        # norm_counts_dict = {}
        # log_norm_counts_dict = {}
        #
        # clipped_counts_edger_df = robjects.DataFrame(clipped_counts_edger)
        # norm_clipped_counts = edger.cpm(clipped_counts_edger_df, normalized_lib_sizes=True)
        # log_norm_clipped_counts = edger.cpm(clipped_counts_edger_df, log=True)
        #
        # norm_clipped_dict = {}
        # log_norm_clipped_dict = {}
        #
        # for i in range(1, norm_counts.ncol + 1):
        #     norm_counts_dict[norm_counts.colnames[i - 1]] = list(norm_counts.rx(True, i))
        #     log_norm_counts_dict[log_norm_counts.colnames[i - 1]] = list(log_norm_counts.rx(True, i))
        #
        # for j in range(1, norm_clipped_counts.ncol + 1):
        #     norm_clipped_dict[norm_clipped_counts.colnames[j - 1]] = list(norm_clipped_counts.rx(True, j))
        #     log_norm_clipped_dict[norm_clipped_counts.colnames[j - 1]] = list(log_norm_clipped_counts.rx(True, j))
        #
        # norm_counts_df = pd.DataFrame(norm_counts_dict)
        # # pandas data frame of normalized counts
        # norm_counts_df = pd.concat([self.read_counts["chr"], self.read_counts["bin"], norm_counts_df], axis=1)
        #
        # norm_clipped_counts_df = pd.DataFrame(norm_clipped_dict)
        # norm_clipped_counts_df = pd.concat([self.read_counts["chr"], self.read_counts["bin"], norm_clipped_counts_df],
        #                                    axis=1)
        # # print(norm_clipped_counts_df)
        # # to transform float counts, that are not truthful for read counts, into integers
        # for col in norm_counts_df:
        #     if col != "chr" and col != "bin":
        #         norm_counts_df[[col]] = norm_counts_df[[col]].astype(int)
        #
        # log_norm_counts_df = pd.DataFrame(log_norm_counts_dict)
        # # pandas data frame of log normalized counts (to be used in fold change calc)
        # log_norm_counts_df = pd.concat([self.read_counts["chr"], self.read_counts['bin'], log_norm_counts_df], axis=1)

        # print(self.norm)

        # log_norm_clip_df = pd.DataFrame(log_norm_clipped_dict)
        # log_norm_clip_df = pd.concat([self.read_counts["chr"], self.read_counts['bin'], log_norm_clip_df], axis=1)

        # clipped_count_df = pd.DataFrame(clipped_count)  # clipped_count filled before
        # norm_clip = {}
        # log_norm_clip = {}
        #
        # for col in clipped_count_df:
        #     # for normalization scope (dividing by the sum(counts) of sample of interest)
        #     read_counts_col = col[:col.find("_cig_filt")]
        #     norm_clip[col] = clipped_count_df[col] / (sum(self.read_counts[read_counts_col]) / 1000000)
        #     log_norm_clip[col] = []
        #     for row in norm_clip[col]:  # this line caused the wrong calculation of fold change in clipped counts,
        #         # norm_clip is the correct structure on which iterate and not clipped_count_df
        #         if row == 0:  # int is useful to avoid that bin with norm count at 0.07... are considered as
        #             # significant in fold change calculation with respect to a 0 count in the reference
        #             log_norm_clip[col].append(0)
        #         else:
        #             log_norm_clip[col].append(math.log2(row))
        #     update_bar += 1
        #     norm_bar.update(update_bar)
        #
        # norm_clip_df = pd.DataFrame(norm_clip)
        # norm_clip_df = pd.concat([self.read_counts["chr"], self.read_counts["bin"], norm_clip_df], axis=1)
        #
        # # for col in norm_clip_df:
        # #     if col != "chr" and col != "bin":
        # #         norm_clip_df[[col]] = norm_clip_df[[col]].astype(int)
        # #         print(norm_clip_df[col])
        #
        # log_norm_clip_df = pd.DataFrame(log_norm_clip)
        # log_norm_clip_df = pd.concat([self.read_counts["chr"], self.read_counts['bin'], log_norm_clip_df], axis=1)
        #

        return self.norm, self.log_norm, self.norm_clip, self.log_norm_clip, self.norm_unmapped

    def calc_fold_change(self, control_name, pairwise=False):
        # fc = self.log_norm["test6_alignSort_REDONE"] - self.log_norm["reference30x_alignSort"]
        if pairwise:
            print("\nPairwise Fold Change calculation:\n")
            pw_fc_bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
            update_bar = 0
            fc_read_counts = {"chr": self.log_norm["chr"], "bin": self.log_norm["bin"],
                              control_name: self.norm[control_name]}
            for col in self.log_norm:
                if col != control_name and col != "bin" and col != "chr":
                    # a column with the normalized value
                    fc_read_counts[col] = self.norm[col]
                    # self.log_norm[[col]] = self.log_norm[[col]].astype(int)
                    # # print(self.log_norm[[col]])
                    # self.log_norm[[control_name]] = self.log_norm[[control_name]].astype(int)
                    pw_fc = self.log_norm[col] - self.log_norm[control_name]
                    # a column with the fold change value
                    fc_read_counts[col + "-" + control_name] = list(pw_fc)
                    # print(list(pw_fc))
                    update_bar += 1
                    pw_fc_bar.update(update_bar)

            fc_clipped_counts = {"chr": self.log_norm_clip["chr"], "bin": self.log_norm_clip["bin"],
                                 control_name + "_cig_filt": self.norm_clip[control_name + "_cig_filt"]}
            for col in self.log_norm_clip:
                if col != control_name + "_cig_filt" and col != "bin" and col != "chr":
                    fc_clipped_counts[col] = self.norm_clip[col]

                    # self.log_norm_clip[[col]] = self.log_norm_clip[[col]].astype(int)
                    # self.log_norm_clip[[control_name + "_cig_filt"]] = \
                    #     self.log_norm_clip[[control_name + "_cig_filt"]].astype(int)

                    pw_clipped_fc = self.log_norm_clip[col] - self.log_norm_clip[control_name + "_cig_filt"]
                    fc_clipped_counts[col + "-" + control_name] = list(pw_clipped_fc)

                    update_bar += 1
                    pw_fc_bar.update(update_bar)

            fc_df = pd.DataFrame(fc_read_counts)
            # print("\nfc_df\n", fc_df)
            fc_clip_df = pd.DataFrame(fc_clipped_counts)
            # print("\nclipped_fc\n", fc_clip_df)
            # print(fc_df.columns)
            # print(fc_clip_df.columns)
            self.set_fold_change(fc_df)
            fc_df.to_csv("./fc_values.tsv", sep="\t")
            self.set_clipped_fold_change(fc_clip_df)
            fc_clip_df.to_csv("./clipped_fc_values.tsv", sep="\t")

            return self.fold_change, self.clipped_fold_change

        else:
            print("\nOne vs All Fold Change calculation:\n")
            fc_bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
            fc_update_bar = 0
            fc_all = self.log_norm[["chr", "bin"]]
            clones_norm_df = self.log_norm.drop(columns=["chr", "bin", control_name]).mean(axis=1)
            for el in range(len(clones_norm_df)):
                fc_update_bar += 1
                fc_bar.update(fc_update_bar)
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

    def build_summary_dict(self, control_name, fc_value, fc_structure, sig=True, clip=False):
        summary_dict = {"chr": [], "start_pos": [], "end_pos": [], "count": [], "control_count": [],
                        "clone_name": [], "type": [], "fc": [], "fc_values": []}

        for col in fc_structure:
            if col != "chr" and col != "bin" and col != control_name and "-" + control_name in col:
                sig_data_pos = fc_structure[fc_structure[col] > fc_value]
                sig_data_neg = fc_structure[fc_structure[col] < -fc_value]
                sig_data = pd.concat([sig_data_pos, sig_data_neg])
                nosig_data = fc_structure.drop(sig_data.index)

                if not sig:
                    if clip:
                        self.fill_summary_dict(summary_dict, nosig_data, col, control_name, clip=clip)
                    else:
                        self.fill_summary_dict(summary_dict, nosig_data, col, control_name)
                else:
                    if clip:
                        self.fill_summary_dict(summary_dict, sig_data, col, control_name,
                                               clip=clip, data_pos=sig_data_pos, data_neg=sig_data_neg)
                    else:
                        self.fill_summary_dict(summary_dict, sig_data, col, control_name,
                                               data_pos=sig_data_pos, data_neg=sig_data_neg)

        return summary_dict

    def fill_summary_dict(self, summary_dict, data, col, control_name, clip=False, data_pos=None, data_neg=None):

        summary_dict["chr"] += list(data["chr"])
        # print(len(list(data["chr"])))
        summary_dict["start_pos"] += list(data["bin"] * self.bin_size)
        # print(len(data["bin"] * self.bin_size))
        summary_dict["end_pos"] += list((data["bin"] * self.bin_size) + self.bin_size)
        # print(len((data["bin"] * self.bin_size) + self.bin_size))
        summary_dict["count"] += list(data[col[:col.find("-" + control_name)]])
        # print(len(list(data[col[:col.find("-")]])))
        if clip:
            summary_dict["control_count"] += list(data[control_name + "_cig_filt"])
            summary_dict["type"] += list(["clipped_count"] * len(data))
        else:
            summary_dict["control_count"] += list(data[control_name])
            summary_dict["type"] += list(["read_count"] * len(data))
            # print(len(["read_count"] * len(data)))
        # print(len(list(data[control_name])))
        summary_dict["clone_name"] += list([col] * len(data))
        # print(len([col] * len(data)))
        summary_dict["fc_values"] += list(data[col])
        if data_pos is not None and data_neg is not None:
            summary_dict["fc"] += list(["+"] * len(data_pos) + ["-"] * len(data_neg))
            # print(len(["+"] * len(data_pos) + ["-"] * len(data_neg)))
        else:
            summary_dict["fc"] += list("n" * len(data))

    def int_summary(self, summary_df):
        summary_df[["count", "control_count"]] = summary_df[["count", "control_count"]].astype(int)

    def concat_summary(self, sum_prop, sum_clip):
        data_df = pd.concat([sum_prop, sum_clip])
        data_df = data_df.sort_values(by=["chr", "clone_name"])
        return data_df

    def summary_sig_nosig_bins(self, fc, control_name, saving_folder):
        # dataframe for information on significant BINS

        summary_sig_data = self.build_summary_dict(control_name, fc, self.fold_change)
        summary_nosig_data = self.build_summary_dict(control_name, fc, self.fold_change, sig=False)

        summary_sig_clip_data = self.build_summary_dict(control_name, fc, self.clipped_fold_change,
                                                        clip=True)
        summary_nosig_clip_data = self.build_summary_dict(control_name, fc, self.clipped_fold_change,
                                                          sig=False, clip=True)

        sum_sig_bins = pd.DataFrame(summary_sig_data)
        sum_nosig_bins = pd.DataFrame(summary_nosig_data)
        sum_sig_clip_bins = pd.DataFrame(summary_sig_clip_data)
        sum_nosig_clip_bins = pd.DataFrame(summary_nosig_clip_data)

        sig_df = self.concat_summary(sum_sig_bins, sum_sig_clip_bins)
        nosig_df = self.concat_summary(sum_nosig_bins, sum_nosig_clip_bins)

        # ---- uncomment to transform float counts coming from normalization process into integer counts ----
        # self.int_summary(sum_sig_bins)
        # self.int_summary(sum_nosig_bins)
        # ---- uncomment to check the subdataframe of significant bins ----
        # sum_sig_bins.to_csv(saving_folder + "method_mod_sum_sig_bins.tsv", sep="\t")
        # sum_nosig_bins.to_csv(saving_folder + "method_mod_sum_NOsig_bins.tsv", sep="\t")

        # ---- uncomment to transform float counts coming from normalization process into integer counts ----
        # self.int_summary(sum_sig_clip_bins)
        # self.int_summary(sum_nosig_clip_bins)
        # ---- uncomment to check the subdataframe of significant bins ----
        # sum_sig_clip_bins.to_csv(saving_folder + "method_mod_sum_sig_clip_bins.tsv", sep="\t")
        # sum_nosig_clip_bins.to_csv(saving_folder + "/method_mod_sum_NOsig_clip_bins.tsv", sep="\t")

        self.set_sig_bins(sig_df)
        # self.set_sig_clip_bins(sum_sig_clip_bins)
        self.set_no_sig_bins(nosig_df)
        # self.set_no_sig_clip_bins(sum_nosig_clip_bins)

        return self.sig_bins, self.nosig_bins  #, self.sig_clip_bins, self.nosig_clip_bins

    def output_sig_positions(self, fc, control_name, file_output_path):

        header = "# Chromosome positions in which fold-change is > " + str(fc) + \
                 " or < -" + str(fc) + " with respect to: " + control_name + "\n"
        with open(file_output_path + "BRAN" + str(self.bin_size) + "_significant_changes_regions_" + str(fc) +
                  ".tsv", "w") as file:
            file.write(header)
            self.sig_bins.to_csv(path_or_buf=file, sep="\t", index=False)

        nosig_header = "# Chromosome positions in which fold-change is between the two thresholds of " + str(fc) + \
                       " / " + str(fc) + " respect to: " + control_name + "\n"
        with open(file_output_path + "BRAN" + str(self.bin_size) + "_NOT_significant_changes_regions_" + str(fc) +
                  ".tsv", "w") as file:
            file.write(nosig_header)
            self.nosig_bins.to_csv(path_or_buf=file, sep="\t", index=False)

        print("\nSignificant changes regions\n", self.sig_bins)
        print("\nNOT significant changes regions\n", self.nosig_bins)

    def plot(self, saving_folder, saving_format, cigar, unmapped, ref_genome, fc, pairwise, control_name, violin_bar,
             scatter, fold_change_pl, chr_name=None, sample=None):
        visualizer = TestingBinReadVisualizer(self.bin_size, self.read_counts, self.norm, self.log_norm,
                                              self.norm_clip, self.log_norm_clip, self.parameters["unmapped_reads"],
                                              self.norm_unmapped, self.fold_change, self.clipped_fold_change,
                                              saving_folder, saving_format)

        if chr_name and sample:
            if violin_bar:
                print("Sorry but for violin distribution and bar plot the entire set of sample is necessary, "
                      "\nPlease retry without the '-vb' parameter")
            elif scatter and not fold_change_pl:
                visualizer.plot_chr_sample(chr_name, sample, cigar)
                # print("\nok9")
                visualizer.plot_norm_chr_sample(chr_name, sample, cigar)
                # print("\nok10")
            elif fold_change_pl and not scatter:
                visualizer.plot_fold_change_chr_sample(pairwise, fc, chr_name, sample, control_name, cigar)
                # print("\nok15")
            else:
                visualizer.plot_chr_sample(chr_name, sample, cigar)
                # print("\nok9")
                visualizer.plot_norm_chr_sample(chr_name, sample, cigar)
                # print("\nok10")
                visualizer.plot_fold_change_chr_sample(pairwise, fc, chr_name, sample, control_name, cigar)
                # print("\nok15")

        elif chr_name and not sample:
            if violin_bar:
                print("Sorry but for violin distribution and bar plot the entire set of sample is necessary, "
                      "\nPlease retry without the '-vb' parameter")
            elif scatter and not fold_change_pl:
                visualizer.plot_chr(chr_name, cigar)
                # print("\nok11")
                visualizer.plot_norm_chr(chr_name, cigar)
                # print("\nok12")
            elif fold_change_pl and not scatter:
                visualizer.plot_fold_change_chr(pairwise, fc, chr_name, control_name, cigar)
                # print("\nok16")
            else:
                visualizer.plot_chr(chr_name, cigar)
                # print("\nok11")
                visualizer.plot_norm_chr(chr_name, cigar)
                # print("\nok12")
                visualizer.plot_fold_change_chr(pairwise, fc, chr_name, control_name, cigar)
                # print("\nok16")

        elif sample and not chr_name:
            if violin_bar:
                print("Sorry but for violin distribution and bar plot the entire set of sample is necessary, "
                      "\nPlease retry without the '-vb' parameter")
            elif scatter and not fold_change_pl:
                visualizer.plot_sample(sample, cigar)
                # print("\nok13")
                visualizer.plot_norm_sample(sample, cigar)
                # print("\nok14")
            elif fold_change_pl and not scatter:
                visualizer.plot_fold_change_sample(pairwise, fc, sample, control_name, cigar)
                # print("\nok17")
            else:
                visualizer.plot_sample(sample, cigar)
                # print("\nok13")
                visualizer.plot_norm_sample(sample, cigar)
                # print("\nok14")
                visualizer.plot_fold_change_sample(pairwise, fc, sample, control_name, cigar)
                # print("\nok17")

        else:
            if violin_bar and not scatter and not fold_change_pl:
                visualizer.plot_violin()
                # print("\nok1")
                visualizer.plot_bar(cigar, unmapped)
                # print("\nok2")
            elif violin_bar and scatter and not fold_change_pl:
                visualizer.plot_violin()
                # print("\nok1")
                visualizer.plot_bar(cigar, unmapped)
                # print("\nok2")
                visualizer.plot_scatter()
                # print("\nok3")
                visualizer.plot_norm_scatter()
                # print("\nok4")
                visualizer.plot_clipped_scatter()
                # print("\nok5")
                visualizer.plot_norm_clipped_scatter()
                # print("\nok6")
            elif violin_bar and fold_change_pl and not scatter:
                visualizer.plot_violin()
                # print("\nok1")
                visualizer.plot_bar(cigar, unmapped)
                # print("\nok2")
                # visualizer.fold_change_colors()
                visualizer.plot_fold_change(fc, pairwise, control_name)
                # print("\nok7")
                visualizer.plot_clip_fold_change(fc, pairwise, control_name)
                # print("\nok8")

            elif fold_change_pl and not violin_bar and not scatter:
                # visualizer.fold_change_colors()
                visualizer.plot_fold_change(fc, pairwise, control_name)
                # print("\nok7")
                visualizer.plot_clip_fold_change(fc, pairwise, control_name)
                # print("\nok8")
            elif fold_change_pl and scatter and not violin_bar:
                visualizer.plot_scatter()
                # print("\nok3")
                visualizer.plot_norm_scatter()
                # print("\nok4")
                visualizer.plot_clipped_scatter()
                # print("\nok5")
                visualizer.plot_norm_clipped_scatter()
                # print("\nok6")
                # visualizer.fold_change_colors()
                visualizer.plot_fold_change(fc, pairwise, control_name)
                # print("\nok7")
                visualizer.plot_clip_fold_change(fc, pairwise, control_name)
                # print("\nok8")
            elif fold_change_pl and violin_bar and not scatter:
                visualizer.plot_violin()
                # print("\nok1")
                visualizer.plot_bar(cigar, unmapped)
                # print("\nok2")
                # visualizer.fold_change_colors()
                visualizer.plot_fold_change(fc, pairwise, control_name)
                # print("\nok7")
                visualizer.plot_clip_fold_change(fc, pairwise, control_name)
                # print("\nok8")
            elif scatter and not fold_change_pl and not violin_bar:
                visualizer.plot_scatter()
                # print("\nok3")
                visualizer.plot_norm_scatter()
                # print("\nok4")
                visualizer.plot_clipped_scatter()
                # print("\nok5")
                visualizer.plot_norm_clipped_scatter()
                # print("\nok6")
            else:
                visualizer.plot_violin()
                # print("\nok1")
                visualizer.plot_bar(cigar, unmapped)
                # print("\nok2")
                visualizer.plot_scatter()
                # print("\nok3")
                visualizer.plot_norm_scatter()
                # print("\nok4")
                visualizer.plot_clipped_scatter()
                # print("\nok5")
                visualizer.plot_norm_clipped_scatter()
                # print("\nok6")
                # visualizer.fold_change_colors()
                visualizer.plot_fold_change(fc, pairwise, control_name)
                # print("\nok7")
                visualizer.plot_clip_fold_change(fc, pairwise, control_name)
                # print("\nok8")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(usage="%(prog)s [options]",
                                     # formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="",  # before the argument help
                                     epilog="")  # after the argument help
    # parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-bs", "--bin_size",
                        type=int,
                        default=250000,
                        help="The length in bp, of the segments that split the chromosomes equally")

    parser.add_argument("-fl", "--flag_list",
                        nargs="+",
                        default=["99", "147", "163", "83"],
                        help="""A list of the bitwise-flags in BAM/SAM format that identify the reads to be counted 
                        during analyses; if different flags wants to be added, add them as strings 
                        (e.g. "177" "129")""")

    parser.add_argument("-fc", "--fold_change",
                        type=float,
                        default=1.5,
                        help="An number to set as fold-change cut-off for bins significance")

    parser.add_argument("-co", "--control_name",
                        type=str,
                        # required=True,
                        help="""The name of the control sample(s) for the fold_change analysis, it should be the same 
                        name as the column name in the read_counts data structure, so the name of the alignment 
                        file used as a baseline, without the '.bam'""")

    parser.add_argument("-pw", "--pairwise",
                        action="store_true",
                        help="""If specified, the fold change is calculated pairwise between 
                        each sample and the control_sample""")

    parser.add_argument("-c", "--cigar",
                        action="store_true",
                        help="If specified, it allows the application of filters on cigar_string, per read")

    parser.add_argument("-cf", "--cigar_filter",
                        nargs="+",
                        default=["S", "H"],
                        help="""If specified, the reads mapped with soft and hard clipping (S and H) by default, 
                            are taken out form the read counts; it returns a data frame with same structure of the 
                            default one but with 2 different columns for sample, one without counting the reads having
                            these filters and the other counting only the reads having these filters.
                            \n(Specify other filters using e.g. "I" "D")""")

    parser.add_argument("-u", "--unmapped",
                        action="store_true",
                        help="""If specified, also a .txt file is  created, with all the unmapped reads and, as 
                        last raw, the counts for each sample""")

    parser.add_argument("-ch", "--chromosome",
                        default=None,
                        type=str,
                        help="""The number of the chromosome of interest to obtain a restricted view of the data""")

    parser.add_argument("-s", "--sample",
                        default=None,
                        type=str,
                        help="""The name of the clone of interest to obtain a restricted view of the data""")

    parser.add_argument("-bc", "--bin_chromosomes",
                        nargs="+",
                        default=[],
                        type=str,
                        help="""The name of the chromosome for each interesting bin (no repetitions) from which retrieve
                             IDs information""")

    parser.add_argument("-bp", "--bin_positions",
                        action="append",
                        default=[],
                        # type=int,
                        # dest="list",
                        help="""The bin position on the corresponding chromosome, be careful that for each 
                            position/list of positions there is one and only one chromosome""")

    parser.add_argument("-id", "--identifier",
                        action="store_true",
                        help="If the retrieval of read IDs is needed")

    parser.add_argument("-f", "--folder",
                        nargs="+",
                        default=["./"],
                        help="The path to the folder in which are located the files to be analyzed (.bam)")

    parser.add_argument("-op", "--output_pickle",
                        default='./',
                        help="The path to the folder in which search the pickle file already created")

    parser.add_argument("-sf", "--saving_folder",
                        type=str,
                        default="./",
                        help="""Path to the directory in which save all plots""")

    parser.add_argument("-fo", "--saving_format",
                        type=str,
                        default="svg",
                        help="""file format for saved plot images, the choice is between:\n
                             ["svg", "jpeg", "pdf", "png"]\n
                             default format is: svg""")

    parser.add_argument("-vb", "--violin_bar",
                        action="store_true",
                        help="If specified BRAN return and save only the violin and the bar plots")

    parser.add_argument("-sc", "--scatter",
                        action="store_true",
                        help="If specified BRAN return and save only the scatter distribution plots")

    parser.add_argument("-cp", "--fold_change_pl",
                        action="store_true",
                        help="If specified BRAN return and save only the fold change plots")

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
    # TODO add this parameter to main and/or to visualizer
    parser.add_argument("--show",
                        action="store_true",
                        default=False,
                        help="""If specified allows BRAN to also show the plots in an HTML interactive page other then
                        save them in a static mode. By default BRAN saves plots only in a static mode""")
    # ---------------------------------------
    # parser.add_argument("-N", "--Ns_count",
    #                     action="store_true",
    #                     help="""Specify if the Ns counts has to be included in the plot of the read counts""")

    # parser.add_argument("-gi", "--general_info",
    #                     action="store_true",
    #                     help="""If specified, distribution_plot, norm_plot_all, (filtered_plot_all) and
    #                     fold_change_plot are displayed""")

    args = parser.parse_args()
    dict_args = vars(parser.parse_args([]))

    if args.flag_list != dict_args["flag_list"]:
        flags = dict_args["flag_list"] + args.flag_list
    else:
        flags = args.flag_list

    if args.folder != dict_args["folder"]:
        bam_folder = dict_args["folder"] + args.folder
    else:
        bam_folder = args.folder

    if args.identifier:
        if not os.path.exists(args.saving_folder):
            os.mkdir(args.saving_folder)
        else:
            pass
        # print(args.bin_chromosomes)
        # print(args.bin_positions)
        bin_pos = []
        for el in args.bin_positions:
            if len(el) > 1:
                bin_pos.append(el.split(","))
            else:
                bin_pos.append(el)
        # print(bin_pos)
        bin_dictionary = dict(zip(args.bin_chromosomes, bin_pos))
        # print(bin_dictionary)
        ide = TestingBinReadIdentifier(args.bin_size,
                                       flags,
                                       bam_folder,
                                       args.saving_folder,
                                       bin_dictionary,
                                       args.cigar,
                                       args.cigar_filter)
        # ide.load_bam()
        ide.get_read_ids()

    else:
        if args.control_name:
            analyzer = TestingBinReadAnalyzer(bam_folder,
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
                os.makedirs(args.saving_folder)
            else:
                pass

            plots_folder = "{}plots/{}/".format(args.saving_folder, str(args.bin_size))
            os.makedirs(plots_folder, exist_ok=True)
            if not os.path.exists(plots_folder):
                os.mkdir(plots_folder)
            else:
                pass
            analyzer.sorted_df(args.saving_folder)

            analyzer.normalize_bins(args.control_name, args.saving_folder)
            analyzer.calc_fold_change(args.control_name, args.pairwise)
            analyzer.summary_sig_nosig_bins(args.fold_change, args.control_name, args.saving_folder)
            analyzer.output_sig_positions(args.fold_change, args.control_name, args.saving_folder)
            # exit(1)
            analyzer.plot(plots_folder, args.saving_format, args.cigar,
                          args.unmapped, args.reference, args.fold_change, args.pairwise, args.control_name,
                          args.violin_bar, args.scatter, args.fold_change_pl, chr_name=args.chromosome,
                          sample=args.sample)
        else:
            print("Argument '-co/--control_name' not passed or wrong, "
                  "it has to be passed and exist in order to calculate the fold change")

    # ---------------------------complete_file_pickle-------------------------------------------------------------------
    # with open("../all_samples_pickles/BRAN250000_df.p", "rb") as input_param:
    #     param = pickle.load(input_param)
    #     chr6 = param["read_counts"][param["read_counts"]["chr"] == "CH.chr6"]
    #     chr1 = param["read_counts"][param["read_counts"]["chr"] == "CH.chr1"]
    #     chr19 = param["read_counts"][param["read_counts"]["chr"] == "CH.chr19"]
    #     chr8 = param["read_counts"][param["read_counts"]["chr"] == "CH.chr8"]
    #     chr13 = param["read_counts"][param["read_counts"]["chr"] == "CH.chr13"]
    #
    #     # print(chr6[["Ch15_3_1_ref_IlluminaPE_aligns_Primary", "bin"]])
    #     print("all\n", param["read_counts"][["Ch15_3_1_ref_IlluminaPE_aligns_Primary",
    #     "Ch15_2_1_ref_IlluminaPE_aligns_Primary"]])
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

    # ----------------------------------tries---------------------------------------------------------------------------
    # analyzer.plot_sh_clipping(args.control_name, args.saving_folder)
    #
    #
    # analyzer.normalize_bins()
    # # analyzer.get_fold_change(args.control_name)
    # fc = analyzer.get_fold_change(args.control_name)
    # analyzer.get_sig_pos(args.fold_change, args.p_value)
    # analyzer.get_no_sig_pos(args.fold_change, args.p_value)
    #
    # counts = params["read_counts"]
    # start = 15309000 // args.bin_size
    # # print(start)
    # end = 15309706 // args.bin_size
    # # print(end)
    #
    # data = counts[counts["chr"] == "CH.chr1"]
    # data = data[data['bin'] >= start]
    # data = data[data['bin'] <= end]
    # print(data)
    # # fc_start_filt = data.index.iloc[0]
    # # fc_end_filt = data.index.iloc[-1]
    #
    # # for clone in fc:
    # #     print(fc[clone].iloc[fc_start_filt:fc_end_filt + 1])
    #
    # # exit(1)
