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
from plotly.subplots import make_subplots


class TestingBinReadVisualizer:
    def __init__(self, bin_size, counts, norm, log_norm, norm_clip, log_norm_clip, unmapped,
                 norm_unmapped, fold_change, clip_fold_change, saving_folder, saving_format, template):
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
        self.template = template
        self.color_palette = px.colors.qualitative.T10 + px.colors.qualitative.Pastel + px.colors.qualitative.Safe

    def sorted_chromosomes(self, column):
        convert = lambda text: int(text) if text.isdigit() else text
        alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        return sorted(column, key=alphanum_key)

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
                                                  fillcolor="rgb(245,245,245)",  # "rgb(240, 248, 255)",  # ~light_mint_green
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

    def add_threshold_fc(self, fig, fc):

        fig.add_shape(go.layout.Shape(type="line",
                                      x0=0,
                                      y0=fc,
                                      x1=len(self.read_counts) * self.bin_size,
                                      y1=fc))

        fig.add_shape(go.layout.Shape(type="line",
                                      x0=0,
                                      y0=-fc,
                                      x1=len(self.read_counts) * self.bin_size,
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
                                        line_color=self.color_palette[i],
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
                                        line_color=self.color_palette[j],
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

    def scatter_layout(self, fig, title):
        fig.update_layout(title=title,
                          template=self.template,
                          legend_orientation="h",
                          legend=dict(x=-0.01, y=1.05))

    def fold_change_layout(self, fig, title, colors):
        fig.update_layout(title=title,
                          template=colors,
                          legend_orientation="h")

    def plot_scatter(self, ref_genome, ns=False, fig=go.Figure()):
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

        col_list = list(self.read_counts.columns)
        hover_pos = self.read_counts["bin"] * self.bin_size

        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Raw_Read_Count_Per_Bin")

        for i in range(len(col_list)):
            if col_list[i] != "chr" and col_list[i] != "bin" and "cig_filt" not in col_list[i]:
                fig.add_trace(go.Scatter(x=list(self.read_counts.index * self.bin_size),
                                         y=self.read_counts[col_list[i]],
                                         hovertext=hover_pos,
                                         hovertemplate=
                                         "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                         mode="markers",
                                         # name=str(col_list[i][:col_list[i].find("_Illumina")])))
                                         name=str(col_list[i])))

                self.scatter_layout(fig, title="Raw_Read Counts - All Clones - All Chromosomes - Bin Size: " +
                                               str(self.bin_size))

        # make a subplot if N trace is required
        # if ns:
        #     self.add_ns_trace(fig, reference=reference)

        self.plot_background(fig)

        fig.show()

        self.saving_plot(fig, description="scatter_all_counts_" + str(self.bin_size))

    def plot_norm_scatter(self, ref_genome, ns=False, fig=go.Figure()):
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
                fig.add_trace(go.Scatter(x=list(self.norm_counts.index * self.bin_size),
                                         y=self.norm_counts[col],
                                         hovertext=hover_pos,
                                         hovertemplate=
                                         "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y:.0f}",
                                         mode="markers",
                                         # name=col[:col.find("_Illumina")]))
                                         name=col))

                self.scatter_layout(fig, title="Normalized Read Counts - Clone: all - Chr: all - Bin Size: " +
                                               str(self.bin_size))

        # if ns:
        #     self.add_ns_trace(fig, reference=reference)

        self.plot_background(fig)

        fig.show()

        self.saving_plot(fig, description="scatter_norm_all_counts_" + str(self.bin_size))

    def plot_clipped_scatter(self, ref_genome, ns=False, fig=go.Figure()):
        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Raw_Clipped_Read_counts")

        for col in self.read_counts:
            if "cig_filt" in col:
                hover_pos = self.read_counts["bin"] * self.bin_size
                fig.add_trace(go.Scatter(x=list(self.read_counts[col].index * self.bin_size),
                                         y=self.read_counts[col],
                                         hovertext=hover_pos,
                                         hovertemplate=
                                         "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                         hoverinfo="text",
                                         mode="markers",
                                         name=col))

                self.scatter_layout(fig, title="Raw_Clipped Read Counts - All Clones - All Chromosomes - Bin Size: " +
                                               str(self.bin_size))

        self.plot_background(fig)

        fig.show()

        self.saving_plot(fig, description="scatter_clip_counts" + str(self.bin_size))

    def plot_norm_clipped_scatter(self, ref_genome, ns=False, fig=go.Figure()):
        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Norm_Clipped_Read_counts")

        for col in self.read_counts:
            if "cig_filt" in col:
                hover_pos = self.read_counts["bin"] * self.bin_size
                fig.add_trace(go.Scatter(x=list(self.norm_clip_counts[col].index * self.bin_size),
                                         y=self.norm_clip_counts[col],
                                         hovertext=hover_pos,
                                         hovertemplate=
                                         "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                         hoverinfo="text",
                                         mode="markers",
                                         name=col))

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
        fig.add_trace(go.Scatter(x=hover_pos,
                                 y=counts_chr[sample],
                                 hovertext=hover_pos,
                                 hovertemplate=
                                 "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                 mode="markers",
                                 name=sample))

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

            fig2.add_trace(go.Scatter(x=hover_pos,
                                      y=counts_chr[sample + "_cig_filt"],
                                      hovertext=hover_pos,
                                      hovertemplate=
                                      "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                      mode="markers",
                                      name=sample + "clipped"))

            self.scatter_layout(fig2, title="Clipped Read Counts - Bin Size: {} <br> Clone: {} <br> Chr: {}".format(
                str(self.bin_size),
                sample,
                str(chr_name)))

            fig2.show()

            self.saving_plot(fig2, description="scatter_clipped_counts_{}_{}_{}".format(str(chr_name),
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

        fig.add_trace(go.Scatter(x=hover_pos,  # list(single_chrom.index * self.bin_size),
                                 y=norm_counts_chr[sample],
                                 hovertext=hover_pos,
                                 hovertemplate=
                                 "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y:.0f}",
                                 mode="markers",
                                 name=sample))

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

            fig2.add_trace(go.Scatter(x=hover_pos,
                                      y=norm_clipped_counts_chr[sample + "_cig_filt"],
                                      hovertext=hover_pos,
                                      hovertemplate=
                                      "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                      mode="markers",
                                      name=sample + "clipped"))

            self.scatter_layout(fig2, title="Normalized Clipped Read Counts - Bin Size: {} "
                                            "<br> Clone: {} "
                                            "<br> Chr: {}".format(str(self.bin_size),
                                                                  sample,
                                                                  str(chr_name)))

            fig2.show()

            self.saving_plot(fig2, description="scatter_norm_clipped_counts_{}_{}_{}".format(str(chr_name),
                                                                                             sample,
                                                                                             str(self.bin_size)))

    def fold_change_colors(self):
        colors = px.colors.qualitative.T10
        colors[2] = "rgb(135,197,35)"  # light green
        colors[4] = "rgb(184, 0, 58)"  # dark magenta
        self.template.layout["colorway"] = colors
        return self.template

    def plot_fold_change(self, fc, pairwise, control_name):
        """"""
        fig = go.Figure()
        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="log2_Fold-Change")
        description = ""
        colors = self.fold_change_colors()
        for col in self.fold_change:
            if col != "chr" and col != "bin":
                sig_data_pos = self.fold_change[self.fold_change[col] > fc]
                sig_data_neg = self.fold_change[self.fold_change[col] < -fc]

                sig_bins = pd.concat([sig_data_pos, sig_data_neg])
                no_sig_bins = self.fold_change.drop(list(sig_data_pos.index) + list(sig_data_neg.index))

                hover_pos_sig = sig_bins["bin"] * self.bin_size
                hover_pos_no_sig = no_sig_bins["bin"] * self.bin_size

                if pairwise:
                    fig.add_trace(go.Scatter(x=list(sig_bins.index * self.bin_size),
                                             y=sig_bins[col],
                                             mode="markers",
                                             hovertext=hover_pos_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             name=col[:col.find("-")]))
                    fig.add_trace(go.Scatter(x=list(no_sig_bins.index * self.bin_size),
                                             y=no_sig_bins[col],
                                             mode="markers",
                                             marker=dict(size=5, color="rgb(176, 196, 222)"),
                                             hovertext=hover_pos_no_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             # showlegend=False,
                                             # legendgroup="group",
                                             # name=col[:col.find("_Illumina")] + "_no_significant"))
                                             name=col[:col.find("-")] + "_no_significant"))

                    self.fold_change_layout(fig,
                                            title="Each <i>vs</i> {}<br>Pairwise log2 Fold Change - Bin Size: {} - "
                                                  "Threshold_FC: ".format(control_name,
                                                                          str(self.bin_size),
                                                                          str(fc)),
                                            colors=colors)
                    description = "pw_fold_change_" + str(self.bin_size)
                else:
                    fig.add_trace(go.Scatter(x=list(sig_bins.index * self.bin_size),
                                             y=sig_bins[col],
                                             mode="markers",
                                             hovertext=hover_pos_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             name="Significant Differences"))
                    fig.add_trace(go.Scatter(x=list(no_sig_bins.index * self.bin_size),
                                             y=no_sig_bins[col],
                                             mode="markers",
                                             marker=dict(size=5, color="rgb(176, 196, 222)"),
                                             hovertext=hover_pos_no_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             # legendgroup="group",
                                             # showlegend=False,
                                             name="Not Significant Differences"))

                    self.fold_change_layout(fig, title="All <i>vs</i> {}<br> log2 Fold Change - Bin Size: {} - "
                                                       "Threshold_FC: ".format(control_name,
                                                                               str(self.bin_size),
                                                                               str(fc)),
                                            colors=colors)
                    description = "fold_change_" + str(self.bin_size)

        self.add_threshold_fc(fig, fc)
        self.plot_background(fig)

        fig.show()

        self.saving_plot(fig, description=description)

    def plot_clip_fold_change(self, fc, pairwise, control_name):
        fig = go.Figure()
        fig.update_xaxes(title_text="Genome_Position")
        fig.update_yaxes(title_text="Clipped_log2_Fold-Change")
        description = ""
        colors = self.fold_change_colors()
        for col in self.clip_fold_change:
            if col != "bin" and col != "chr":
                sig_clip_data_pos = self.clip_fold_change[self.clip_fold_change[col] > fc]
                sig_clip_data_neg = self.clip_fold_change[self.clip_fold_change[col] < -fc]

                sig_clip_bins = pd.concat([sig_clip_data_pos, sig_clip_data_neg])

                no_sig_clip_bins = self.clip_fold_change.drop(list(sig_clip_data_pos.index) +
                                                              list(sig_clip_data_neg.index))

                hover_pos_sig = sig_clip_bins["bin"] * self.bin_size
                hover_pos_no_sig = no_sig_clip_bins["bin"] * self.bin_size

                if pairwise:
                    fig.add_trace(go.Scatter(x=list(sig_clip_bins.index * self.bin_size),
                                             y=sig_clip_bins[col],
                                             mode="markers",
                                             hovertext=hover_pos_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             # legendgroup="group",
                                             name=col[:col.find("-")]))
                    fig.add_trace(go.Scatter(x=list(no_sig_clip_bins.index * self.bin_size),
                                             y=no_sig_clip_bins[col],
                                             mode="markers",
                                             marker=dict(size=5, color="rgb(176, 196, 222)"),
                                             hovertext=hover_pos_no_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             # legendgroup="group",
                                             name=col[:col.find("-")] + "_no_significant"))

                    self.fold_change_layout(fig,
                                            title="Each <i>vs</i> {}<br>".format(control_name) +
                                                  "Clipped Reads Pairwise log2 Fold Change - " +
                                                  "Bin Size: {} - Threshold_FC: {}".format(str(self.bin_size),
                                                                                           str(fc)),
                                            colors=colors)
                    description = "pw_clip_fold_change_" + str(self.bin_size)

                else:
                    fig.add_trace(go.Scatter(x=list(sig_clip_bins.index * self.bin_size),
                                             y=sig_clip_bins[col],
                                             mode="markers",
                                             hovertext=hover_pos_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             name="Significant Differences"))
                    fig.add_trace(go.Scatter(x=list(no_sig_clip_bins.index * self.bin_size),
                                             y=no_sig_clip_bins[col],
                                             mode="markers",
                                             marker=dict(size=5, color="rgb(176, 196, 222)"),
                                             hovertext=hover_pos_no_sig,
                                             hovertemplate=
                                             "<b>Chrom_position</b>: %{hovertext:,}" + "<br>Count: %{y}",
                                             # legendgroup="group",
                                             name="Not Significant Differences"))

                    self.fold_change_layout(fig,
                                            title="All <i>vs</i> {}<br>".format(control_name) +
                                                  "Clipped Reads log2 Fold Change - " +
                                                  "Bin Size: " + str(self.bin_size) + " - Threshold_FC: " + str(fc),
                                            colors=colors)

                    description = "clip_fold_change_" + str(self.bin_size)

        self.add_threshold_fc(fig, fc)
        self.plot_background(fig)

        fig.show()

        self.saving_plot(fig, description=description)
