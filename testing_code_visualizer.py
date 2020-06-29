# testing BinReadVisualizer

# This class furnish different ways to visualize data coming from BinReadAnalyzer

import os
import argparse
import math
import re
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots


class TestingBinReadVisualizer:
    def __init__(self, bin_size, counts, norm, log_norm, norm_clip, log_norm_clip, unmapped,
                 norm_unmapped, fold_change, clip_fold_change, saving_folder, saving_format):
        self.bin_size = bin_size
        self.read_counts = counts
        self.norm_counts = norm
        self.log_norm_counts = log_norm
        self.norm_clip_counts = norm_clip
        self.log_norm_clip = log_norm_clip
        self.unmapped = unmapped
        self.norm_unmapped = norm_unmapped
        self.fold_change = fold_change
        self.clip_fold_change = clip_fold_change
        self.saving_folder = saving_folder
        self.saving_format = saving_format

    def sorted_chromosomes(self, column):
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(column, key=alphanum_key)

    def color_palette(self, fc_colors=False):
        template = pio.templates["seaborn"]
        template.layout["font"]["color"] = "#2a3f5f"  # to verify: "#DC143C"
        color_palette = ["rgb(183, 10, 48)", "rgb(255, 186, 8)", "rgb(63, 163, 197)", "rgb(3, 84, 99)",
                         "rgb(110, 29, 93)", "rgb(235, 181, 155)", "rgb(188, 178, 215)", "rgb(196, 235, 112)",
                         "rgb(196, 90, 140)", "rgb(32, 100, 186)", "rgb(255, 0, 75)", "rgb(109, 89, 122)",
                         "rgb(165, 255, 186)", "rgb(196, 51, 122)", "rgb(255, 132, 11)", "rgb(121, 110, 192)"]

        if fc_colors:
            fc_palette = []
            for color in color_palette:
                fc_palette += [color, color]
            template.layout["colorway"] = fc_palette
        else:
            template.layout["colorway"] = color_palette

        return template

    def saving_plot(self, fig, description):
        acceptable_formats = ["svg", "jpeg", "pdf", "png"]
        if self.saving_format not in acceptable_formats:
            print("Sorry the specified image format is not acceptable, please try again with one of these:",
                  acceptable_formats)
        else:
            fig.write_image(self.saving_folder + description + "." + self.saving_format, width=1920, height=1080)

    def plot_background(self, fig):  # , df_counts

        coordinates_x = []
        coordinates_y = []
        chromosomes = []
        start = 0

        for ch in self.read_counts["chr"]:
            # if ch[ch.find("c"):] not in chromosomes:
            if ch not in chromosomes:
                chromosomes.append(ch)
        # if the chromosomes are not sorted in ascending order, using the sort_
        # chromosomes method, which sort alphanumeric strings
        chromosomes = self.sorted_chromosomes(chromosomes)
        # print(chromosomes)
        for chrom in chromosomes:
            single_df = self.read_counts[self.read_counts["chr"] == chrom]
            # print(single_df)
            # here start and end are updated every time, to allow concatenation of chromosomes on the plot
            # the average position within the interval of each chromosome take in account only the length of the
            # interval itself
            length = (single_df["bin"].iloc[-1] + 1) * self.bin_size
            tmp_end = start + length

            avg = start + length / 2
            coordinates_x.append(avg)
            coordinates_y.append(-5)

            if chrom[-1].isdigit() is True:
                if int(chrom[-1]) % 2 == 0:
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
                                                  fillcolor="rgb(245,245,245)",
                                                  # "rgb(240, 248, 255)",  # ~light_mint_green
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
                                              fillcolor="grey",
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

    def add_threshold_fc(self, fig, fc, len_x_axis):

        fig.add_shape(go.layout.Shape(type="line",
                                      x0=0,
                                      y0=fc,
                                      x1=len_x_axis * self.bin_size,
                                      y1=fc))

        fig.add_shape(go.layout.Shape(type="line",
                                      x0=0,
                                      y0=-fc,
                                      x1=len_x_axis * self.bin_size,
                                      y1=-fc))

        fig.update_shapes(dict(xref="x",
                               yref="y",
                               line=dict(color="crimson",  # dark red
                                         width=1)))

    def plot_violin(self):
        """"""

        fig = make_subplots(rows=1, cols=2, subplot_titles=("Normalized Counts",
                                                            "Soft_Hard Clipped Read Counts"))
        i = 0
        j = 0
        hover_pos = self.norm_counts["bin"] * self.bin_size
        for col in self.read_counts:
            if "cig_filt" in col:
                fig.add_trace(go.Violin(y=self.norm_clip_counts[col],
                                        box_visible=True,
                                        meanline_visible=True,
                                        hovertext=hover_pos,
                                        text=self.norm_counts["chr"],
                                        hovertemplate=
                                        "<b>Chrom</b>: %{text}" +
                                        "<br><b>Position</b>: %{hovertext:,}" +
                                        "<br>Count: %{y:.0f}",
                                        # fillcolor="rgb(0,139,139)",
                                        # line_color=possible_color_violin[i],
                                        # name=col[:col.find("_Illumina")]),
                                        line_color=self.color_palette().layout["colorway"][i],
                                        name=col),
                              row=1,
                              col=2)
                fig.update_yaxes(title_text="Norm_Clipped_counts", row=1, col=2)
                fig.update_xaxes(tickangle=45)
                i += 1

            elif "cig_filt" not in col and col != "chr" and col != "bin":
                fig.add_trace(go.Violin(y=self.norm_counts[col],
                                        box_visible=True,
                                        meanline_visible=True,
                                        hovertext=hover_pos,
                                        text=self.norm_counts["chr"],
                                        hovertemplate=
                                        "<b>Chrom</b>: %{text}" +
                                        "<br><b>Position</b>: %{hovertext:,}" +
                                        "<br>Count: %{y:.0f}",
                                        # fillcolor="rgb(255,160,122)",
                                        # line_color=possible_color_violin[j],
                                        # name=col[:col.find("_Illumina")]),
                                        line_color= self.color_palette().layout["colorway"][j],
                                        name=col),
                              row=1,
                              col=1)
                fig.update_yaxes(title_text="Norm_Read_counts", row=1, col=1)

                j += 1

        fig.update_layout(showlegend=False,
                          title="Comparison Between Read Counts and Only Clipped Read Counts per Sample" +
                                "- Bin Size: " + str(self.bin_size), title_x=0.5)

        fig.update_traces(opacity=0.75)

        fig.show()

        self.saving_plot(fig, description="violin_plot_" + str(self.bin_size))

    def plot_bar(self, cigar, unmapped):

        summary_read_counts = []
        summary_clipped_counts = []
        summary_unmapped_reads = list(self.unmapped.values())
        x_labels = list(self.unmapped.keys())
        perc_read_counts = []
        total_reads = []
        fig = go.Figure()

        for col in self.read_counts:
            if col != "chr" and col != "bin":
                if "cig_filt" in col:
                    summary_clipped_counts.append(sum(self.read_counts[col]))
                else:
                    summary_read_counts.append(sum(self.read_counts[col]))

        for i in range(len(summary_read_counts)):
            total_reads.append(summary_read_counts[i] + summary_clipped_counts[i] + summary_unmapped_reads[i])
            perc_read_counts.append((summary_read_counts[i] / total_reads[i]) * 100)

        fig.add_trace(go.Bar(x=x_labels,
                             y=perc_read_counts,
                             text=["{:.3f}%".format(el) for el in perc_read_counts],
                             textposition="auto",
                             marker_color="#4C78A8",  # dark cyan
                             name="%properly_mapped_reads"))
        if unmapped:
            for i in range(len(summary_unmapped_reads)):
                summary_unmapped_reads[i] = (summary_unmapped_reads[i] / total_reads[i]) * 100

            fig.add_trace(go.Bar(x=x_labels,
                                 y=summary_unmapped_reads,
                                 text=["{:.3f}%".format(el) for el in summary_unmapped_reads],
                                 textposition="auto",
                                 marker_color="#FF9900",  # dark yellow/orange
                                 name="%unmapped_reads"))

        if cigar:
            for i in range(len(summary_clipped_counts)):
                summary_clipped_counts[i] = (summary_clipped_counts[i] / total_reads[i]) * 100

            fig.add_trace(go.Bar(x=x_labels,
                                 y=summary_clipped_counts,
                                 text=["{:.3f}%".format(el) for el in summary_clipped_counts],
                                 textposition="auto",
                                 marker_color="rgb(220,20,60)",  # crimson
                                 name="%clipped_reads"))

        fig.update_layout(barmode="stack",
                          title_text="% Clipped Reads - Unmapped Reads - Properly Mapped Reads Per Sample- Bin Size: " +
                                     str(self.bin_size))

        fig.update_traces(opacity=0.75)

        fig.show()

        self.saving_plot(fig, description="bar_chart_" + str(self.bin_size))

    def scatter_traces(self, fig, x_coord, y_coord, hover_text, trace_name):
        fig.add_trace(go.Scatter(x=x_coord,  # itw it begins from 0 any t.
                                 y=y_coord,
                                 hovertext=hover_text,
                                 text=self.read_counts["chr"],
                                 hovertemplate=
                                 "<b>%{text}</b>" +
                                 "<br>Chrom_position</b>: %{hovertext:,}" +
                                 "<br>Count: %{y}",
                                 mode="markers",
                                 name=trace_name))

    def scatter_layout(self, fig, title):
        fig.update_layout(title=title,
                          template=self.color_palette(),
                          legend_orientation="h",
                          legend=dict(x=-0.01, y=1.05))

    def plot_scatter(self, ref_genome=False, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of raw_read_counts
        of all chromosomes and all samples

        Args:
            ref_genome (bool): true if the reference is declared
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of counts
        """

        col_list = list(self.read_counts.columns)
        hover_pos = self.read_counts["bin"] * self.bin_size

        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Raw_Read_Count_Per_Bin")

        for i in range(len(col_list)):
            if col_list[i] != "chr" and col_list[i] != "bin" and "cig_filt" not in col_list[i]:
                self.scatter_traces(fig,
                                    list(self.read_counts.index * self.bin_size),  # otw it begins from 0 any t.
                                    self.read_counts[col_list[i]],
                                    hover_pos,
                                    str(col_list[i]))

                self.scatter_layout(fig, title="Raw_Read Counts - All Clones - All Chromosomes - Bin Size: " +
                                               str(self.bin_size))

        # make a subplot if N trace is required
        # if ns:
        #     self.add_ns_trace(fig, reference=reference)

        self.plot_background(fig)

        fig.show()

        self.saving_plot(fig, description="scatter_all_counts_" + str(self.bin_size))

    def plot_norm_scatter(self, ref_genome=False, ns=False, fig=go.Figure()):
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

        hover_pos = self.norm_counts["bin"] * self.bin_size
        for col in self.norm_counts:
            if col != "chr" and col != "bin":
                self.scatter_traces(fig,
                                    list(self.norm_counts.index * self.bin_size),
                                    self.norm_counts[col],
                                    hover_pos,
                                    col)

                self.scatter_layout(fig, title="Normalized Read Counts - Clone: all - Chr: all - Bin Size: " +
                                               str(self.bin_size))

        # if ns:
        #     self.add_ns_trace(fig, reference=reference)

        self.plot_background(fig)

        fig.show()

        self.saving_plot(fig, description="scatter_norm_all_counts_" + str(self.bin_size))

    def plot_clipped_scatter(self, ref_genome=False, ns=False, fig=go.Figure()):
        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Raw_Clipped_Read_counts")

        for col in self.read_counts:
            if "cig_filt" in col:
                hover_pos = self.read_counts["bin"] * self.bin_size
                self.scatter_traces(fig,
                                    list(self.read_counts[col].index * self.bin_size),
                                    self.read_counts[col],
                                    hover_pos,
                                    col)

                self.scatter_layout(fig, title="Raw_Clipped Read Counts - All Clones - All Chromosomes - Bin Size: " +
                                               str(self.bin_size))

        self.plot_background(fig)

        fig.show()

        self.saving_plot(fig, description="scatter_clip_counts" + str(self.bin_size))

    def plot_norm_clipped_scatter(self, ref_genome=False, ns=False, fig=go.Figure()):
        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Norm_Clipped_Read_counts")

        for col in self.read_counts:
            if "cig_filt" in col:
                hover_pos = self.read_counts["bin"] * self.bin_size
                self.scatter_traces(fig,
                                    list(self.norm_clip_counts[col].index * self.bin_size),
                                    self.norm_clip_counts[col],
                                    hover_pos,
                                    col)

                self.scatter_layout(fig, title="Normalized Clipped Read Counts - All Clones - "
                                               "All Chromosomes - Bin Size: " +
                                               str(self.bin_size))

        self.plot_background(fig)

        fig.show()

        self.saving_plot(fig, description="scatter_norm_clip_counts" + str(self.bin_size))

    def plot_chr_sample(self, chr_name, sample, cigar, ref_genome=False, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of raw read_counts
        for a specific chromosome and a specific sample

        Args:
            ref_genome (bool): true if the reference is declared
            chr_name (int): a number representing the chromosome of interest
            sample (str): the name of the sample of interest (it corresponds
                          to the name of the column in the data structure)
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of counts
        """
        fig.update_xaxes(title_text="Chromosomes_Position")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        counts_chr = self.read_counts[self.read_counts["chr"].str.contains(chr_name) == True]
        hover_pos = counts_chr["bin"] * self.bin_size
        self.scatter_traces(fig,
                            hover_pos,
                            counts_chr[sample],
                            hover_pos,
                            sample)

        self.scatter_layout(fig, title="Read Counts - Bin Size: {} <br> Clone: {} <br> Chr: {}".format(
            str(self.bin_size),
            sample,
            str(chr_name)))

        fig.show()

        self.saving_plot(fig, description="scatter_counts_{}_{}_{}".format(str(chr_name),
                                                                           sample,
                                                                           str(self.bin_size)))

        if cigar:
            fig2 = go.Figure()
            fig2.update_xaxes(title_text="Chromosomes_Position")
            fig2.update_yaxes(title_text="Read_Clipped_Count_Per_Bin")

            self.scatter_traces(fig2,
                                hover_pos,
                                counts_chr[sample + "_cig_filt"],
                                hover_pos,
                                sample + "clipped")

            self.scatter_layout(fig2, title="Clipped Read Counts - Bin Size: {} <br> Clone: {} <br> Chr: {}".format(
                str(self.bin_size),
                sample,
                str(chr_name)))

            fig2.show()

            self.saving_plot(fig2, description="scatter_clip_counts_{}_{}_{}".format(str(chr_name),
                                                                                     sample,
                                                                                     str(self.bin_size)))

    def plot_norm_chr_sample(self, chr_name, sample, cigar, ref_genome=False, ns=False,
                             fig=go.Figure()):
        """This method allows to obtain a scatter-plot of normalized_read_counts
        of a specific chromosome of a specific sample

        Args:
            ref_genome (bool): true if the reference is declared
            chr_name (int): a number representing the chromosome of interest
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

        norm_counts_chr = self.norm_counts[self.norm_counts["chr"].str.contains(chr_name)]
        hover_pos = norm_counts_chr["bin"] * self.bin_size

        self.scatter_traces(fig,
                            hover_pos,  # list(single_chrom.index * self.bin_size),
                            norm_counts_chr[sample],
                            hover_pos,
                            sample)

        # if ns:
        #     self.add_ns_trace(fig, reference=reference, chrom=chrom)
        self.scatter_layout(fig, title="Normalized Read Counts - Bin Size: {} "
                                       "<br> Clone: {} "
                                       "<br> Chr: {}".format(str(self.bin_size),
                                                             sample,
                                                             str(chr_name)))

        fig.show()
        self.saving_plot(fig, description="scatter_norm_counts_{}_{}_{}".format(str(chr_name),
                                                                                sample,
                                                                                str(self.bin_size)))

        if cigar:
            fig2 = go.Figure()
            fig2.update_xaxes(title_text="Chromosomes_Position")
            fig2.update_yaxes(title_text="Read_Clipped_Count_Per_Bin")

            norm_clipped_counts_chr = self.norm_clip_counts[self.norm_clip_counts["chr"].str.contains(chr_name)]
            hover_pos = norm_clipped_counts_chr["bin"] * self.bin_size

            self.scatter_traces(fig2,
                                hover_pos,
                                norm_clipped_counts_chr[sample + "_cig_filt"],
                                hover_pos,
                                sample + "clipped")

            self.scatter_layout(fig2, title="Normalized Clipped Read Counts - Bin Size: {} "
                                            "<br> Clone: {} "
                                            "<br> Chr: {}".format(str(self.bin_size),
                                                                  sample,
                                                                  str(chr_name)))

            fig2.show()

            self.saving_plot(fig2, description="scatter_norm_clip_counts_{}_{}_{}".format(str(chr_name),
                                                                                          sample,
                                                                                          str(self.bin_size)))

    def plot_chr(self, chr_name, cigar, ref_genome=False, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of raw_read_counts
        of all samples, but for a specific chromosome of interest

        Args:
            ref_genome (bool): true if the reference is declared
            chr_name (int): a number representing the chromosome of interest
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of counts
        """
        fig.update_xaxes(title_text="Chromosomes_Position")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        counts_chr = self.read_counts[self.read_counts["chr"].str.contains(chr_name) == True]
        hover_pos = counts_chr["bin"] * self.bin_size
        for col in counts_chr:
            if col != "chr" and col != "bin" and "cig_filt" not in col:
                self.scatter_traces(fig,
                                    hover_pos,  # list(single_chrom.index * self.bin_size),
                                    counts_chr[col],
                                    hover_pos,
                                    col)
        # if ns:
        #     self.add_ns_trace(fig, reference=reference, chrom=chrom)

        self.scatter_layout(fig, title="Read Counts - Bin Size: {} "
                                       "<br> Chr: {}".format(str(self.bin_size),
                                                             str(chr_name)))

        fig.show()

        self.saving_plot(fig, description="scatter_read_counts_{}_all_{}".format(str(chr_name),
                                                                                 str(self.bin_size)))

        if cigar:
            fig2 = go.Figure()
            fig2.update_xaxes(title_text="Chromosomes_Position")
            fig2.update_yaxes(title_text="Read_Clipped_Count_Per_Bin")

            for col in counts_chr:
                if col != "chr" and col != "bin" and "cig_filt" in col:
                    self.scatter_traces(fig2,
                                        hover_pos,
                                        counts_chr[col],
                                        hover_pos,
                                        col)

            self.scatter_layout(fig2, title="Clipped Read Counts - Bin Size: {} "
                                            "<br> Chr: {}".format(str(self.bin_size),
                                                                  str(chr_name)))

            fig2.show()

            self.saving_plot(fig2, description="scatter_clip_read_counts_{}_all_{}".format(str(chr_name),
                                                                                           str(self.bin_size)))

    def plot_norm_chr(self, chr_name, cigar, ref_genome=False, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of normalized_read_counts
        of a specific chromosome for all samples

        Args:
            ref_genome (bool): true if the reference is declared
            chr_name (int): a number representing the chromosome of interest
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of normalized counts
        """
        fig.update_xaxes(title_text="Chromosome_position")
        fig.update_yaxes(title_text="Norm_Read_Count_Per_Bin")

        counts_chr = self.norm_counts[self.norm_counts["chr"].str.contains(chr_name) == True]
        hover_pos = counts_chr["bin"] * self.bin_size
        for col in counts_chr:
            if col != "chr" and col != "bin":
                self.scatter_traces(fig,
                                    hover_pos,  # list(single_chrom.index * self.bin_size),
                                    counts_chr[col],
                                    hover_pos,
                                    col)

        # if ns:
        #     self.add_ns_trace(fig, reference=reference, chrom=chrom)

        self.scatter_layout(fig, title="Normalized Read Counts - Bin Size: {}"
                                       "<br> Chr: {}".format(str(self.bin_size), chr_name))

        fig.show()

        self.saving_plot(fig, description="scatter_counts_chr_{}_all_{}".format(chr_name, str(self.bin_size)))

        if cigar:
            fig2 = go.Figure()
            fig2.update_xaxes(title_text="Chromosome_position")
            fig2.update_yaxes(title_text="Norm_Clipped_Read_Count_Per_Bin")

            clip_counts_chr = self.norm_clip_counts[self.norm_clip_counts["chr"].str.contains(chr_name) == True]
            for col in clip_counts_chr:
                if col != "chr" and col != "bin":
                    self.scatter_traces(fig2,
                                        hover_pos,  # list(single_chrom.index * self.bin_size),
                                        clip_counts_chr[col],
                                        hover_pos,
                                        col)

            self.scatter_layout(fig2, title="Normalized Clipped Read Counts - Bin Size: {}"
                                            "<br> Chr: {}".format(str(self.bin_size), chr_name))

            fig2.show()

            self.saving_plot(fig2, description="scatter_clip_counts_chr_{}_all_{}".format(chr_name, str(self.bin_size)))

    def plot_sample(self, sample, cigar, ref_genome=False, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of raw_read_counts
        of all chromosomes, but for a specific sample of interest

        Args:
            ref_genome (bool): true if the reference is declared
            sample (str): the name of the sample of interest (it corresponds
                          to the name of the column in the data structure)
            ns (bool): by default sets to False, but is one want to include
                       the Ns count trace, it has to set to True
            fig (obj): is a go.Figure() object for the building of the plot

        Returns:
            A scatter-plot of counts
        """
        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Read_Count_Per_Bin")

        hover_pos = self.read_counts["bin"] * self.bin_size

        self.scatter_traces(fig,
                            list(self.read_counts.index * self.bin_size),  # otw it begins from 0 for any chr
                            self.read_counts[sample],
                            hover_pos,  # useful to know the bin pos within the chromosome
                            sample)

        # if ns:
        #     self.add_ns_trace(fig, reference=reference)

        self.scatter_layout(fig, title="Read Counts - Bin Size: {}"
                                       "<br> Clone: {}".format(str(self.bin_size),
                                                               sample))

        self.plot_background(fig)

        fig.show()

        self.saving_plot(fig, description="scatter_counts_{}_{}".format(sample, str(self.bin_size)))

        if cigar:
            fig2 = go.Figure()
            fig2.update_xaxes(title_text="Genome_Position")
            fig2.update_yaxes(title_text="Read_Count_Per_Bin")

            hover_pos = self.read_counts["bin"] * self.bin_size

            self.scatter_traces(fig2,
                                list(self.read_counts.index * self.bin_size),
                                self.read_counts[sample + "_cig_filt"],
                                hover_pos,  # useful to know the bin pos within the chromosome
                                sample + "_cig_filt")

            # if ns:
            #     self.add_ns_trace(fig2, reference=reference)

            self.scatter_layout(fig2, title="Clipped Read Counts - Bin Size: {}"
                                            "<br> Clone: {}".format(str(self.bin_size),
                                                                    sample))

            self.plot_background(fig2)

            fig2.show()

            self.saving_plot(fig2, description="scatter_clip_counts_{}_{}".format(sample, str(self.bin_size)))

    def plot_norm_sample(self, sample, cigar, ref_genome=False, ns=False, fig=go.Figure()):
        """This method allows to obtain a scatter-plot of normalized_read_counts
        in all chromosomes of a specific sample

        Args:
            ref_genome (bool): true if the reference is declared
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

        hover_pos = self.norm_counts["bin"] * self.bin_size

        self.scatter_traces(fig,
                            list(self.norm_counts.index * self.bin_size),
                            self.norm_counts[sample],
                            hover_pos,
                            sample)

        # if ns:
        #     self.add_ns_trace(fig, reference=reference)

        self.plot_background(fig)

        self.scatter_layout(fig, title="Normalized Read Counts - Bin Size: {}"
                                       "<br> Clone: {}".format(str(self.bin_size), sample))

        fig.show()

        self.saving_plot(fig, description="scatter_norm_counts_{}_{}".format(sample, str(self.bin_size)))

        if cigar:
            fig2 = go.Figure()
            fig2.update_xaxes(title_text="Genome_Position")
            fig2.update_yaxes(title_text="Norm_Clipped_Read_Count_Per_Bin")

            hover_pos = self.norm_clip_counts["bin"] * self.bin_size

            self.scatter_traces(fig2,
                                list(self.norm_clip_counts.index * self.bin_size),
                                self.norm_clip_counts[sample + "_cig_filt"],
                                hover_pos,
                                sample + "_cig_filt")

            # if ns:
            #     self.add_ns_trace(fig2, reference=reference)

            self.plot_background(fig2)

            self.scatter_layout(fig2, title="Normalized Clipped Read Counts - Bin Size: {}"
                                            "<br> Clone: {}".format(str(self.bin_size), sample))

            fig2.show()

            self.saving_plot(fig2, description="scatter_norm_clip_counts_{}_{}".format(sample, str(self.bin_size)))

    # def fold_change_colors(self):
    #     palette = self.template.layout["colorway"]
    #     fc_palette = []
    #     for color in palette:
    #         fc_palette += [color, color]
    #     fc_palette_array = np.array(fc_palette)
    #     return fc_palette_array

    def fold_change_traces(self, pairwise, fig, x_sig, y_sig,
                           x_no_sig, y_no_sig, hover_sig, hover_no_sig, text_sig, text_no_sig, trace_name):

        if pairwise:
            fig.add_trace(go.Scatter(x=x_sig,
                                     y=y_sig,
                                     mode="markers",
                                     hovertext=hover_sig,
                                     text=text_sig,
                                     hovertemplate="<b>%{text}</b>" +
                                                   "<br>Chrom_position</b>: %{hovertext:,}" +
                                                   "<br>Count: %{y}",
                                     name=trace_name))

            fig.add_trace(go.Scatter(x=x_no_sig,
                                     y=y_no_sig,
                                     mode="markers",
                                     opacity=0.1,
                                     hovertext=hover_no_sig,
                                     text=text_no_sig,
                                     hovertemplate="<b>%{text}</b>" +
                                                   "<br>Chrom_position</b>: %{hovertext:,}" +
                                                   "<br>Count: %{y}",
                                     name=trace_name + "_<b>not_significant</b>"))

        else:
            fig.add_trace(go.Scatter(x=x_sig,
                                     y=y_sig,
                                     mode="markers",
                                     hovertext=hover_sig,
                                     text=text_sig,
                                     hovertemplate="<b>%{text}</b>" +
                                                   "<br>Chrom_position</b>: %{hovertext:,}" +
                                                   "<br>Count: %{y}",
                                     name="Significant Differences"))

            fig.add_trace(go.Scatter(x=x_no_sig,
                                     y=y_no_sig,
                                     mode="markers",
                                     opacity=0.1,
                                     hovertext=hover_no_sig,
                                     text=text_no_sig,
                                     hovertemplate="<b>%{text}</b>" +
                                                   "<br>Chrom_position</b>: %{hovertext:,}" +
                                                   "<br>Count: %{y}",
                                     name="Not Significant Differences"))

    def fold_change_layout(self, fig, title):
        fig.update_layout(title=title,
                          template=self.color_palette(fc_colors=True),
                          legend_orientation="h")

    def plot_fold_change(self, fc, pairwise, control_name):
        """"""
        fig = go.Figure()
        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="log2_Fold-Change")
        description = ""
        # colors = self.fold_change_colors()
        for col in self.fold_change:
            if col != "chr" and col != "bin":
                sig_data_pos = self.fold_change[self.fold_change[col] > fc]
                sig_data_neg = self.fold_change[self.fold_change[col] < -fc]

                sig_bins = pd.concat([sig_data_pos, sig_data_neg])
                no_sig_bins = self.fold_change.drop(list(sig_data_pos.index) + list(sig_data_neg.index))

                hover_pos_sig = sig_bins["bin"] * self.bin_size
                hover_pos_no_sig = no_sig_bins["bin"] * self.bin_size

                x_axis_sig = list(sig_bins.index * self.bin_size)
                x_axis_no_sig = list(no_sig_bins.index * self.bin_size)

                self.fold_change_traces(pairwise, fig, x_axis_sig, sig_bins[col], x_axis_no_sig, no_sig_bins[col],
                                        hover_pos_sig, hover_pos_no_sig, sig_bins["chr"], no_sig_bins["chr"],
                                        col[:col.find("-")])
                if pairwise:
                    self.fold_change_layout(fig,
                                            title="Each <i>vs</i> {}<br>Pairwise log2 Fold Change - Bin Size: {} - "
                                                  "Threshold_FC: {}".format(control_name,
                                                                            str(self.bin_size),
                                                                            str(fc)))

                    description = "pw_fold_change_" + str(self.bin_size)

                else:

                    self.fold_change_layout(fig, title="All Mean <i>vs</i> {}<br>log2 Fold Change - Bin Size: {} - "
                                                       "Threshold_FC: {}".format(control_name,
                                                                                 str(self.bin_size),
                                                                                 str(fc)))
                    description = "fold_change_" + str(self.bin_size)

        self.add_threshold_fc(fig, fc, len(self.fold_change))
        self.plot_background(fig)

        fig.show()

        self.saving_plot(fig, description=description)

    def plot_clip_fold_change(self, fc, pairwise, control_name):
        fig = go.Figure()
        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Clipped_log2_Fold-Change")
        description = ""
        # colors = self.fold_change_colors()
        for col in self.clip_fold_change:
            if col != "bin" and col != "chr":
                sig_clip_data_pos = self.clip_fold_change[self.clip_fold_change[col] > fc]
                sig_clip_data_neg = self.clip_fold_change[self.clip_fold_change[col] < -fc]

                sig_clip_bins = pd.concat([sig_clip_data_pos, sig_clip_data_neg])

                no_sig_clip_bins = self.clip_fold_change.drop(list(sig_clip_data_pos.index) +
                                                              list(sig_clip_data_neg.index))

                hover_pos_sig = sig_clip_bins["bin"] * self.bin_size
                hover_pos_no_sig = no_sig_clip_bins["bin"] * self.bin_size
                x_axis_sig = list(sig_clip_bins.index * self.bin_size)
                x_axis_no_sig = list(no_sig_clip_bins.index * self.bin_size)

                self.fold_change_traces(pairwise, fig, x_axis_sig, sig_clip_bins[col], x_axis_no_sig,
                                        no_sig_clip_bins[col], hover_pos_sig, hover_pos_no_sig, sig_clip_bins["chr"],
                                        no_sig_clip_bins["chr"], col[:col.find("-")])

                if pairwise:
                    self.fold_change_layout(fig,
                                            title="Each <i>vs</i> {}<br>".format(control_name) +
                                                  "Clipped Reads Pairwise log2 Fold Change - " +
                                                  "Bin Size: {} - Threshold_FC: {}".format(str(self.bin_size),
                                                                                           str(fc)))
                    description = "pw_clip_fold_change_" + str(self.bin_size)

                else:

                    self.fold_change_layout(fig,
                                            title="All <i>vs</i> {}<br>".format(control_name) +
                                                  "Clipped Reads log2 Fold Change - " +
                                                  "Bin Size: " + str(self.bin_size) + " - Threshold_FC: " + str(fc))

                    description = "clip_fold_change_" + str(self.bin_size)

        self.add_threshold_fc(fig, fc, len_x_axis=len(self.clip_fold_change))
        self.plot_background(fig)

        fig.show()

        self.saving_plot(fig, description=description)

    def plot_fold_change_chr_sample(self, pairwise, fc, chr_name, sample, control_name, cigar):
        """"""
        if pairwise:
            sample_pw = sample + "-" + control_name
            fig = go.Figure()

            fig.update_xaxes(title_text="Chromosome_Position")
            fig.update_yaxes(title_text="Fold-Change")

            fc_chr = self.fold_change[self.fold_change["chr"].str.contains(chr_name) == True]

            sig_bins_pos = fc_chr[fc_chr[sample_pw] > fc]
            sig_bins_neg = fc_chr[fc_chr[sample_pw] < -fc]
            sig_bins = pd.concat([sig_bins_pos, sig_bins_neg])
            not_sig_bins = fc_chr.drop(list(sig_bins_pos.index) + list(sig_bins_neg.index))

            hover_pos_sig = sig_bins["bin"] * self.bin_size  # no needs for the list because only one chr
            hover_pos_no_sig = not_sig_bins["bin"] * self.bin_size

            self.fold_change_traces(pairwise, fig, hover_pos_sig, sig_bins[sample_pw], hover_pos_no_sig,
                                    not_sig_bins[sample_pw], hover_pos_sig, hover_pos_no_sig, chr_name, chr_name, sample)

            self.fold_change_layout(fig,
                                    title="Pairwise Fold Change {} <i>vs</i> {}"
                                          "<br> Bin Size: {} - Threshold_FC: {} "
                                          "- Chr: {} ".format(sample,
                                                              control_name,
                                                              str(self.bin_size),
                                                              fc,
                                                              chr_name))

            x_axis = len(sig_bins) + len(not_sig_bins)
            self.add_threshold_fc(fig, fc, x_axis)

            fig.show()

            self.saving_plot(fig, description="pw_fc_{}_{}_{}".format(chr_name, sample_pw, str(self.bin_size)))

            if cigar:
                fig2 = go.Figure()
                sample_clip = sample + "_cig_filt-" + control_name

                fig2.update_xaxes(title_text="Chromosome_Position")
                fig2.update_yaxes(title_text="Clipped_log2_Fold-Change")

                clip_fc_chr = self.clip_fold_change[self.clip_fold_change["chr"].str.contains(chr_name) == True]

                sig_bins_pos = clip_fc_chr[clip_fc_chr[sample_clip] > fc]
                sig_bins_neg = clip_fc_chr[clip_fc_chr[sample_clip] < -fc]
                sig_bins = pd.concat([sig_bins_pos, sig_bins_neg])
                not_sig_bins = clip_fc_chr.drop(list(sig_bins_pos.index) + list(sig_bins_neg.index))

                hover_pos_sig = sig_bins["bin"] * self.bin_size
                hover_pos_no_sig = not_sig_bins["bin"] * self.bin_size

                self.fold_change_traces(pairwise, fig2, hover_pos_sig, sig_bins[sample_clip], hover_pos_no_sig,
                                        not_sig_bins[sample_clip], hover_pos_sig, hover_pos_no_sig, chr_name, chr_name,
                                        sample + "_cig_filt")

                self.fold_change_layout(fig2,
                                        title="{} <i>vs</i> {} - Chr: {} "
                                              "<br> Clipped Pairwise Fold Change - "
                                              "Bin Size: {} - Threshold_FC: {} ".format(sample,
                                                                                        control_name,
                                                                                        chr_name,
                                                                                        str(self.bin_size),
                                                                                        fc))

                x_axis = len(sig_bins) + len(not_sig_bins)
                self.add_threshold_fc(fig2, fc, x_axis)

                fig2.show()

                self.saving_plot(fig2, description="clipped_pw_fc_{}_{}_{}".format(chr_name,
                                                                                   sample_clip,
                                                                                   str(self.bin_size)))
        else:
            print("""ATTENTION: if parameter '-pw' not given, it's impossible to retrieve graphical information 
                  on single sample fold-change. \nPlease TRY AGAIN specifying '-pw' or '--pairwise' in command line""")

    def plot_fold_change_chr(self, pairwise, fc, chr_name, control_name, cigar):
        """"""
        fig = go.Figure()
        fig.update_xaxes(title_text="Chromosome_Position")
        fig.update_yaxes(title_text="log2_Fold-Change")

        fig2 = go.Figure()
        fig2.update_xaxes(title_text="Chromosome_Position")
        fig2.update_yaxes(title_text="Clipped_log2_Fold-Change")

        fc_chr = self.fold_change[self.fold_change["chr"].str.contains(chr_name) == True]
        clip_fc_chr = self.clip_fold_change[self.clip_fold_change["chr"].str.contains(chr_name) == True]

        for col in fc_chr:
            if col != "chr" and col != "bin":
                sig_bins_pos = fc_chr[fc_chr[col] > fc]
                sig_bins_neg = fc_chr[fc_chr[col] < -fc]
                sig_bins = pd.concat([sig_bins_pos, sig_bins_neg])
                not_sig_bins = fc_chr.drop(list(sig_bins_pos.index) + list(sig_bins_neg.index))

                hover_pos_sig = sig_bins["bin"] * self.bin_size  # no need for the list because is only one chr
                hover_pos_no_sig = not_sig_bins["bin"] * self.bin_size

                self.fold_change_traces(pairwise, fig, hover_pos_sig, sig_bins[col], hover_pos_no_sig, not_sig_bins[col],
                                        hover_pos_sig, hover_pos_no_sig, chr_name, chr_name, col[:col.find("-")])

        if cigar:
            for col in clip_fc_chr:
                if col != "chr" and col != "bin":
                    sig_bins_pos = clip_fc_chr[clip_fc_chr[col] > fc]
                    sig_bins_neg = clip_fc_chr[clip_fc_chr[col] < -fc]
                    sig_bins = pd.concat([sig_bins_pos, sig_bins_neg])
                    not_sig_bins = clip_fc_chr.drop(list(sig_bins_pos.index) + list(sig_bins_neg.index))

                    hover_pos_sig = sig_bins["bin"] * self.bin_size
                    hover_pos_no_sig = not_sig_bins["bin"] * self.bin_size

                    self.fold_change_traces(pairwise, fig2, hover_pos_sig, sig_bins[col], hover_pos_no_sig,
                                            not_sig_bins[col], hover_pos_sig, hover_pos_no_sig, chr_name, chr_name,
                                            col[:col.find("-")])

        self.add_threshold_fc(fig, fc, len_x_axis=len(fc_chr))
        self.add_threshold_fc(fig2, fc, len_x_axis=len(clip_fc_chr))

        if pairwise:

            self.fold_change_layout(fig,
                                    title="Pairwise Fold Change Each <i>vs</i> {}"
                                          "<br> Bin Size: {} - Threshold_FC:  {} "
                                          "- Chr: {}".format(control_name,
                                                             str(self.bin_size),
                                                             fc,
                                                             chr_name))

            self.fold_change_layout(fig2,
                                    title="Clipped Pairwise Fold Change Each <i>vs</i> {}"
                                          "<br> Bin Size: {} - Threshold_FC:  {} "
                                          "- Chr: {}".format(control_name,
                                                             str(self.bin_size),
                                                             fc,
                                                             chr_name))

            self.saving_plot(fig, description="pw_fc_{}_{}".format(chr_name, str(self.bin_size)))
            self.saving_plot(fig2, description="clip_pw_fc_{}_{}".format(chr_name, str(self.bin_size)))

        else:
            self.fold_change_layout(fig,
                                    title="Fold Change All Mean <i>vs</i> {}"
                                          "<br> Bin Size: {} - Threshold_FC:  {} "
                                          "- Chr: {}".format(control_name,
                                                             str(self.bin_size),
                                                             fc,
                                                             chr_name))

            self.fold_change_layout(fig2,
                                    title="Clipped Fold Change All Mean <i>vs</i> {}"
                                          "<br> Bin Size: {} - Threshold_FC:  {} "
                                          "- Chr: {}".format(control_name,
                                                             str(self.bin_size),
                                                             fc,
                                                             chr_name))

            self.saving_plot(fig, description="fc_{}_{}".format(chr_name, str(self.bin_size)))
            self.saving_plot(fig2, description="clip_fc_{}_{}".format(chr_name, str(self.bin_size)))

        fig.show()

        fig2.show()

    def plot_fold_change_sample(self, pairwise, fc, sample, control_name, cigar):
        """"""
        if pairwise:
            fig = go.Figure()

            fig.update_xaxes(title_text="Genome_Position")
            fig.update_yaxes(title_text="Fold_change")

            col = sample + "-" + control_name
            sig_bins_pos = self.fold_change[["chr", "bin", col]][self.fold_change[col] > fc]
            sig_bins_neg = self.fold_change[["chr", "bin", col]][self.fold_change[col] < -fc]
            sig_bins = pd.concat([sig_bins_pos, sig_bins_neg])
            not_sig_bins = self.fold_change[["chr", "bin", col]].drop(list(sig_bins_pos.index) +
                                                                      list(sig_bins_neg.index))

            hover_pos_sig = sig_bins["bin"] * self.bin_size
            hover_pos_no_sig = not_sig_bins["bin"] * self.bin_size
            x_axis_sig = list(sig_bins.index * self.bin_size)
            x_axis_no_sig = list(not_sig_bins.index * self.bin_size)

            self.fold_change_traces(pairwise, fig, x_axis_sig, sig_bins[col], x_axis_no_sig, not_sig_bins[col],
                                    hover_pos_sig, hover_pos_no_sig, sig_bins["chr"], not_sig_bins["chr"], sample)

            self.fold_change_layout(fig, title="Pairwise Fold Change - All Chromosomes - Bin Size: {} "
                                               "- Threshold_fc: {} <br> {} <i>vs</i> {}".format(str(self.bin_size),
                                                                                                fc,
                                                                                                sample,
                                                                                                control_name))

            self.add_threshold_fc(fig, fc, len_x_axis=len(self.fold_change[col]))
            self.plot_background(fig)

            fig.show()

            self.saving_plot(fig, description="pw_fc_{}_{}".format(sample, str(self.bin_size)))

            if cigar:
                sample_clip = sample + "_cig_filt-" + control_name
                fig2 = go.Figure()

                fig2.update_xaxes(title_text="Chromosome_Position")
                fig2.update_yaxes(title_text="Clipped_log2_Fold-Change")

                sig_bins_pos = self.clip_fold_change[["chr", "bin", sample_clip]][
                    self.clip_fold_change[sample_clip] > fc]
                sig_bins_neg = self.clip_fold_change[["chr", "bin", sample_clip]][
                    self.clip_fold_change[sample_clip] < -fc]
                sig_bins = pd.concat([sig_bins_pos, sig_bins_neg])
                not_sig_bins = self.clip_fold_change[["chr", "bin", sample_clip]].drop(list(sig_bins_pos.index) +
                                                                                       list(sig_bins_neg.index))

                hover_pos_sig = sig_bins["bin"] * self.bin_size
                hover_pos_no_sig = not_sig_bins["bin"] * self.bin_size
                x_axis_sig = list(sig_bins.index * self.bin_size)
                x_axis_no_sig = list(not_sig_bins.index * self.bin_size)

                self.fold_change_traces(pairwise, fig2, x_axis_sig, sig_bins[sample_clip], x_axis_no_sig,
                                        not_sig_bins[sample_clip], hover_pos_sig, hover_pos_no_sig, sig_bins["chr"],
                                        not_sig_bins["chr"], sample + "_cig_filt")

                self.fold_change_layout(fig2,
                                        title="{} <i>vs</i> {}"
                                              "<br> Clipped Pairwise Fold Change - "
                                              "Bin Size: {} - Threshold_FC: {} ".format(sample,
                                                                                        control_name,
                                                                                        str(self.bin_size),
                                                                                        fc))

                self.add_threshold_fc(fig2, fc, len_x_axis=len(self.clip_fold_change))
                self.plot_background(fig2)
                fig2.show()

                self.saving_plot(fig2, description="clipped_pw_fc_{}_{}".format(sample_clip,
                                                                                str(self.bin_size)))

        else:
            print("""ATTENTION: if parameter '-pw' not given, it's impossible to retrieve graphical information 
                  on single sample fold-change. \nPlease TRY AGAIN specifying '-pw' or '--pairwise' in command line""")
