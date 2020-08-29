# BinReadAnalyzer class
#
# This class aims to connect BinReadCounter, BinReadVisualizer and BinReadIdentifier all in one, to analyze data
# structure, visualize data plots and retrieve IDs information of certain specified reads.
# One can consider this class as the very effector of the analysis on the genome of the of interest.
# For analysis is intended, first of all, the normalization of the raw_read_counts and than the detection of
# significant differences in terms of number of reads in bin, that could lead to the identification of part of
# chromosome or even an entire chromosome that have a significantly different count of reads from the origin plant;
# proving that something in that region or chromosome is happened.

import os
import argparse
import math
import progressbar
import re
import plotly.graph_objects as go
# Disable the orca response timeout.
import plotly.io._orca
import retrying
import pandas as pd
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from BinReadCounter import BinReadCounter
from BinReadVisualizer import BinReadVisualizer
from BinReadIdentifier import BinReadIdentifier
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

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
    return sorted(column, key=alphanum_key)


class BinReadAnalyzer:
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
        self.sig_clip_bins = None
        self.nosig_clip_bins = None

    def set_parameters(self, param):
        """set parameters as actual parameters"""
        self.parameters = param

    def set_read_counts(self, sort_read_counts):
        self.read_counts = sort_read_counts

    def set_norm_dfs(self, norm):
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

    def set_sig_clip_bins(self, sig_clip_data):
        """set the data structure of significant clipped fold-change counts"""
        self.sig_clip_bins = sig_clip_data

    def set_no_sig_clip_bins(self, no_sig_clip_data):
        self.nosig_clip_bins = pd.concat([no_sig_clip_data])

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

    def normalize_bins(self, control_name):
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
        clipped_counts_edger = {}
        print("\n")
        print("Normalization process:\n")

        norm_bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
        update_bar = 0

        for col in self.read_counts:
            if col != "chr" and col != "bin" and "cig_filt" not in col:
                # creates an element of the dictionary [sample_name]: r_vector_of_counts
                read_counts_edger[col] = robjects.IntVector(self.read_counts[col])

            elif "cig_filt" in col:
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
        # -----------------------------------------------------------------------------------------------
        return self.norm, self.log_norm, self.norm_clip, self.log_norm_clip, self.norm_unmapped

    def calc_fold_change(self, control_name, pairwise=False):

        if pairwise:
            print("\nPairwise Fold Change calculation:\n")
            pw_fc_bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
            update_bar = 0
            fc_read_counts = {"chr": self.log_norm["chr"], "bin": self.log_norm["bin"],
                              control_name: self.norm[control_name]}
            fc_clipped_counts = {"chr": self.log_norm_clip["chr"], "bin": self.log_norm_clip["bin"],
                                 control_name + "_cig_filt": self.norm_clip[control_name + "_cig_filt"]}

            for col in self.log_norm:
                if col != control_name and col != "bin" and col != "chr":
                    fc_read_counts[col] = self.norm[col]
                    pw_fc = self.log_norm[col] - self.log_norm[control_name]
                    fc_read_counts[col + "-" + control_name] = list(pw_fc)
                    update_bar += 1
                    pw_fc_bar.update(update_bar)

            for col in self.log_norm_clip:
                if col != control_name + "_cig_filt" and col != "bin" and col != "chr":
                    fc_clipped_counts[col] = self.norm_clip[col]
                    pw_clipped_fc = self.log_norm_clip[col] - self.log_norm_clip[control_name + "_cig_filt"]
                    fc_clipped_counts[col + "-" + control_name] = list(pw_clipped_fc)
                    update_bar += 1
                    pw_fc_bar.update(update_bar)

            fc_df = pd.DataFrame(fc_read_counts)
            fc_clip_df = pd.DataFrame(fc_clipped_counts)

            self.set_fold_change(fc_df)
            print(fc_df)
            self.set_clipped_fold_change(fc_clip_df)

            return self.fold_change, self.clipped_fold_change

        else:
            print("\nOne vs All Fold Change calculation:\n")
            fc_bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
            fc_update_bar = 0
            fc_all = self.log_norm[["chr", "bin"]]
            clones_norm_df = self.log_norm.drop(columns=["chr", "bin", control_name]).mean(axis=1)

            for r in range(len(clones_norm_df)):
                fc_update_bar += 1
                fc_bar.update(fc_update_bar)
            tmp_fc = pd.DataFrame({"fc": clones_norm_df - self.log_norm[control_name]})
            fc_all = pd.concat([fc_all, tmp_fc], axis=1)

            fc_clip_all = self.log_norm_clip[["chr", "bin"]]

            clones_clip_norm_df = self.log_norm_clip.drop(
                columns=["chr", "bin", control_name + "_cig_filt"]).mean(axis=1)

            tmp_clip_fc = pd.DataFrame(
                {"clip_fc": clones_clip_norm_df - self.log_norm_clip[control_name + "_cig_filt"]})
            fc_clip_all = pd.concat([fc_clip_all, tmp_clip_fc], axis=1)

            self.set_fold_change(fc_all)
            self.set_clipped_fold_change(fc_clip_all)
            return self.fold_change, self.clipped_fold_change

    def summary_sig_bins(self, fc, control_name):
        # dataframe for information on significant BINS

        summary_sig_clip_data = {"chr": [], "start_pos": [], "end_pos": [], "count": [], "control_count": [],
                                 "clone_name": [], "type": [], "fc": []}
        summary_sig_data = {"chr": [], "start_pos": [], "end_pos": [], "count": [], "control_count": [],
                            "clone_name": [], "type": [], "fc": []}

        summary_nosig_data = {"chr": [], "start_pos": [], "end_pos": [], "count": [], "control_count": [],
                              "clone_name": [], "type": [], "fc": []}
        summary_nosig_clip_data = {"chr": [], "start_pos": [], "end_pos": [], "count": [], "control_count": [],
                                   "clone_name": [], "type": [], "fc": []}

        for col in self.fold_change.columns:
            if col != "chr" and col != "bin" and col != control_name and "-" in col:
                sig_data_pos = self.fold_change[self.fold_change[col] > fc]
                sig_data_neg = self.fold_change[self.fold_change[col] < -fc]
                sig_data = pd.concat([sig_data_pos, sig_data_neg])
                nosig = self.fold_change.drop(sig_data.index)
                # sig_data.to_csv("./" + col + ".tsv", sep="\t")
                # nosig.to_csv("./" + col + "nosig.tsv", sep="\t")
                print(col)
                print("SIG")
                summary_sig_data["chr"] += list(sig_data["chr"])
                print(len(list(sig_data["chr"])))
                summary_sig_data["start_pos"] += list(sig_data["bin"] * self.bin_size)
                print(len(sig_data["bin"] * self.bin_size))
                summary_sig_data["end_pos"] += list((sig_data["bin"] * self.bin_size) + self.bin_size)
                print(len((sig_data["bin"] * self.bin_size) + self.bin_size))
                summary_sig_data["count"] += list(sig_data[col[:col.find("-")]])
                print(len(list(sig_data[col[:col.find("-")]])))
                summary_sig_data["control_count"] += list(sig_data[control_name])
                print(len(list(sig_data[control_name])))
                summary_sig_data["clone_name"] += list([col] * len(sig_data))
                print(len([col] * len(sig_data)))
                summary_sig_data["type"] += list(["read_count"] * len(sig_data))
                print(len(["read_count"] * len(sig_data)))
                summary_sig_data["fc"] += list(["+"] * len(sig_data_pos) + ["-"] * len(sig_data_neg))
                print(len(["+"] * len(sig_data_pos) + ["-"] * len(sig_data_neg)))

                print("NO SIG")
                summary_nosig_data["chr"] += list(nosig["chr"])
                print(len(list(nosig["chr"])))
                summary_nosig_data["start_pos"] += list(nosig["bin"] * self.bin_size)
                print(len(nosig["bin"] * self.bin_size))
                summary_nosig_data["end_pos"] += list((nosig["bin"] * self.bin_size) + self.bin_size)
                print(len((nosig["bin"] * self.bin_size) + self.bin_size))
                summary_nosig_data["count"] += list(nosig[col[:col.find("-")]])
                print(len(list(nosig[col[:col.find("-")]])))
                summary_nosig_data["control_count"] += list(nosig[control_name])
                print(len(list(nosig[control_name])))
                summary_nosig_data["clone_name"] += list([col] * len(nosig))
                print(len([col] * len(nosig)))
                summary_nosig_data["type"] += list(["read_count"] * len(nosig))
                print(len(["read_count"] * len(nosig)))
                summary_nosig_data["fc"] += list("n" * len(nosig))
                print(len(list("n" * len(nosig))))

        sum_sig_bins = pd.DataFrame(summary_sig_data)
        sum_nosig_bins = pd.DataFrame(summary_nosig_data)

        # ---- to transform float counts coming from normalization process into integer counts ----
        sum_sig_bins[["count", "control_count"]] = sum_sig_bins[["count", "control_count"]].astype(int)
        # ---- to check the subdataframe of significant bins ----
        # sum_sig_bins.to_csv("./sum_sig_bins.tsv", sep="\t")
        # sum_nosig_bins.to_csv("./sum_NOsig_bins.tsv", sep="\t")

        for col in self.clipped_fold_change:
            if col != "chr" and col != "bin" and col != control_name and "-" in col:
                sig_clip_data_pos = self.clipped_fold_change[self.clipped_fold_change[col] > fc]
                sig_clip_data_neg = self.clipped_fold_change[self.clipped_fold_change[col] < -fc]

                sig_clip_data = pd.concat([sig_clip_data_pos, sig_clip_data_neg])
                nosig_clip_data = self.clipped_fold_change.drop(sig_clip_data.index)
                # print(sig_bins)
                summary_sig_clip_data["chr"] += list(sig_clip_data["chr"])
                # print(len(list(sig_clip_data["chr"])))
                summary_sig_clip_data["start_pos"] += (list(sig_clip_data["bin"] * self.bin_size))
                # print(len(sig_clip_data["bin"] * self.bin_size))
                summary_sig_clip_data["end_pos"] += (list((sig_clip_data["bin"] * self.bin_size) + self.bin_size))
                # print(len((sig_clip_data["bin"] * self.bin_size) + self.bin_size))
                summary_sig_clip_data["count"] += (list(sig_clip_data[col[:col.find("-")]]))
                # print(len([col] * len(sig_clip_data)))
                summary_sig_clip_data["control_count"] += (list(sig_clip_data[control_name + "_cig_filt"]))
                summary_sig_clip_data["clone_name"] += [col.replace("_cig_filt", "")] * len(sig_clip_data)
                summary_sig_clip_data["type"] += ["clipped_count"] * len(sig_clip_data)
                # print(len(["clipped_count"] * len(sig_clip_data)))
                summary_sig_clip_data["fc"] += ["+"] * len(sig_clip_data_pos) + ["-"] * len(sig_clip_data_neg)

                summary_nosig_clip_data["chr"] += list(nosig_clip_data["chr"])
                summary_nosig_clip_data["start_pos"] += (list(nosig_clip_data["bin"] * self.bin_size))
                summary_nosig_clip_data["end_pos"] += (list((nosig_clip_data["bin"] * self.bin_size) + self.bin_size))
                summary_nosig_clip_data["count"] += (list(nosig_clip_data[col[:col.find("-")]]))
                summary_nosig_clip_data["control_count"] += (list(nosig_clip_data[control_name + "_cig_filt"]))
                summary_nosig_clip_data["clone_name"] += [col.replace("_cig_filt", "")] * len(nosig_clip_data)
                summary_nosig_clip_data["type"] += ["clipped_count"] * len(nosig_clip_data)
                summary_nosig_clip_data["fc"] += list("n" * len(nosig_clip_data))

        sum_sig_clip_bins = pd.DataFrame(summary_sig_clip_data)
        sum_nosig_clip_bins = pd.DataFrame(summary_nosig_clip_data)
        # ---- to transform float counts coming from normalization process into integer counts ----
        # sum_sig_clip_bins[["count", "control_count"]] = sum_sig_clip_bins[["count", "control_count"]].astype(int)

        # ---- to check the subdataframe of significant bins ----
        # sum_sig_clip_bins.to_csv("./sum_sig_clip_bins.tsv", sep="\t")
        # sum_nosig_clip_bins.to_csv("./sum_NOsig_clip_bins.tsv", sep="\t")

        self.set_sig_bins(sum_sig_bins)
        self.set_sig_clip_bins(sum_sig_clip_bins)
        self.set_no_sig_bins(sum_nosig_bins)
        self.set_no_sig_clip_bins(sum_nosig_clip_bins)
        return self.sig_bins, self.sig_clip_bins, self.nosig_bins, self.nosig_clip_bins

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

    def output_sig_positions(self, fc, control_name, file_output_path):

        sig_df = pd.concat([self.sig_bins, self.sig_clip_bins])
        sig_df = sig_df.sort_values(by=["clone_name", "chr", "start_pos"])

        nosig_df = pd.concat([self.nosig_bins, self.nosig_clip_bins])
        nosig_df = nosig_df.sort_values(by=["clone_name", "chr", "start_pos"])

        header = "# Chromosome positions in which fold-change is > " + str(fc) + \
                 " or < -" + str(fc) + " with respect to: " + control_name + "\n"
        with open(file_output_path + "BRAN" + str(self.bin_size) + "_significant_changes_regions.tsv", "w") as file:
            file.write(header)
            sig_df.to_csv(path_or_buf=file, sep="\t", index=False)

        nosig_header = "# Chromosome positions in which fold-change is between the two thresholds of " + str(fc) + \
                       " / " + str(fc) + " respect to: " + control_name + "\n"
        with open(file_output_path + "BRAN" + str(self.bin_size) + "_NOT_significant_changes_regions.tsv", "w") as file:
            file.write(nosig_header)
            nosig_df.to_csv(path_or_buf=file, sep="\t", index=False)

        print("\nSignificant changes regions\n", sig_df)
        print("\nNOT significant changes regions\n", nosig_df)

    def plot(self, saving_folder, saving_format, cigar, unmapped, ref_genome, fc, pairwise, control_name, violin_bar,
             scatter, fold_change_pl, chr_name=None, sample=None):
        visualizer = BinReadVisualizer(self.bin_size, self.read_counts, self.norm, self.log_norm,
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
                        help="""Path to the directory in which save all outputs.
                         ATTENTION: for plots, new directories is automatically created with default name: 
                         /plots/[bin_size]/""")

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
        ide = BinReadIdentifier(args.bin_size,
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
            if not os.path.exists(args.output_pickle):
                os.mkdir(args.output_pickle)
            else:
                pass

            analyzer = BinReadAnalyzer(bam_folder,
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
                print("\nargs.saving_folder path actually not present. The new folder path has been created\n")
                os.mkdir(args.saving_folder)
            else:
                pass

            plots_folder = "{}/plots/{}/".format(args.saving_folder, str(args.bin_size))
            if not os.path.exists(plots_folder):
                os.mkdir(plots_folder)
            else:
                pass
            analyzer.sorted_df(args.saving_folder)

            # exit(1)
            analyzer.normalize_bins(args.control_name)
            analyzer.calc_fold_change(args.control_name, args.pairwise)
            analyzer.summary_sig_bins(args.fold_change, args.control_name)
            analyzer.output_sig_positions(args.fold_change, args.control_name, args.saving_folder)
            analyzer.plot(plots_folder, args.saving_format, args.cigar,
                          args.unmapped, args.reference, args.fold_change, args.pairwise, args.control_name,
                          args.violin_bar, args.scatter, args.fold_change_pl, chr_name=args.chromosome,
                          sample=args.sample)
        else:
            print("Argument '-co/--control_name' not passed, "
                  "it has to be passed in order for fold_change to be calculated")
# # version 27/08/2020
