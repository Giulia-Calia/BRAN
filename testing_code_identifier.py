# testing BinReadIdentifier

# This class is aimed to retrieve the ID of reads present in the sig_bins file
# given a certain chromosome and a certain position on the chromosome, together with the bin_size,
# it return a new file in which the ID or reads in that bin are stored
# another file is also created with the IDs of all unmapped reads of the sample
import os
import pysam
import progressbar
import pandas as pd


class TestingBinReadIdentifier:
    def __init__(self, bin_size, flags, folder, bins, cigar_filter):
        self.bin_size = bin_size
        self.flags = flags
        self.bam_folder = folder  # a list of folders
        self.bins = bins  # a dictionary of chromosome:bin_str_pos pairs
        self.cigar_filter = cigar_filter

    def mapped_ids(self, cigar):
        header = "chr \t ID \t\t str_pos \t end_pos \t clone_name \t type"
        sig_ids = open("./sig_read_ids.tsv", "w")
        sig_ids.write(header + "\n")
        # unmapped_ids = {"chr": [], "ID": [], "clone_name": []}
        bam = []
        for f in self.bam_folder:
            dir_list = os.listdir(f)
            print("\n\n", f)
            for el in dir_list:
                if el.endswith(".bam"):
                    bam.append(f + el)

                    # reads_bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength)
                    # update_bar = 0
        for file in bam:
            bam_file = pysam.AlignmentFile(file, "rb")
            clone_name = file[file.rfind("/") + 1:file.find(".bam")]
            print("\n", clone_name)
            for ch, bin_str in zip(self.bins.keys(), self.bins.values()):
                for ref in bam_file.references:
                    if ch in ref:
                        for read in bam_file.fetch(contig=ref):
                            if str(read.flag) in self.flags and \
                                    cigar and \
                                    read.cigarstring is not None and \
                                    bin_str <= int(read.reference_start) < (bin_str + self.bin_size):
                                if not any(filters in read.cigarstring for filters in self.cigar_filter):
                                    sig_ids.write(ref + " \t " +
                                                  read.query_name + " \t " +
                                                  clone_name + " \t " +
                                                  str(read.reference_start) + " \t " +
                                                  str(int(read.reference_start) + len(read.query_sequence) - 1) +
                                                  " \t " + "properly_mapped" + "\n")
                                else:
                                    sig_ids.write(ref + " \t " +
                                                  read.query_name + " \t\t " +
                                                  clone_name + " \t " +
                                                  str(read.reference_start) + " \t " +
                                                  str(int(read.reference_start) + len(read.query_sequence) - 1) +
                                                  " \t " + "clipped" + "\n")
                            elif str(read.flag) in self.flags and \
                                    not cigar and \
                                    read.cigarstring is not None and \
                                    bin_str <= int(read.reference_start) < (bin_str + self.bin_size):
                                sig_ids.write(ref + " \t " +
                                              read.query_name + " \t " +
                                              clone_name + " \t " +
                                              str(read.reference_start) + " \t " +
                                              str(int(read.reference_start) + len(read.query_sequence) - 1) +
                                              " \t " + "mapped" + "\n")

        sig_ids.close()
        #                 for ref in bam_file.references:
        #                     if ch in ref:
        #                         for read in bam_file.fetch(contig=ref):
        #                                                    # start=self.bins[ch],
        #                                                    # end=self.bins[ch] + self.bin_size - 1):
        #
        #                             if str(read.flag) in self.flags:
        #                                 if cigar and \
        #                                         read.cigarstring is not None and \
        #                                         self.bins[ch] <= int(read.reference_start) < (self.bins[ch] + self.bin_size):
        #                                     if not any(filt in read.cigarstring for filt in self.cigar_filter):
        #                                         # always the left-most position in POS field in SAM
        #                                         ids["chr"].append(ref)
        #                                         ids["ID"].append(read.query_name)
        #                                         ids["clone_name"].append(clone_name)
        #                                         ids["str_pos"].append(int(read.reference_start))
        #                                         ids["end_pos"].append(int(read.reference_start) +
        #                                                               len(read.query_sequence) - 1)
        #                                         ids["type"].append("properly_mapped")
        #                                     else:
        #                                         ids["chr"].append(ref)
        #                                         ids["ID"].append(read.query_name)
        #                                         ids["clone_name"].append(clone_name)
        #                                         ids["str_pos"].append(int(read.reference_start))
        #                                         ids["end_pos"].append(int(read.reference_start) +
        #                                                               len(read.query_sequence) - 1)
        #                                         ids["type"].append("clipped")
        #
        #                                 elif not cigar:
        #                                     ids["chr"].append(ref)
        #                                     ids["ID"].append(read.query_name)
        #                                     ids["clone_name"].append(clone_name)
        #                                     ids["str_pos"].append(int(read.reference_start))
        #                                     ids["end_pos"].append(int(read.reference_start) +
        #                                                           len(read.query_sequence) - 1)
        #                                     ids["type"].append("no_cig_filter")
        #                             else:
        #                                 continue
        #
        #             bam_file.close()
        # # print(len(ids["chr"]))
        # # print(len(ids["ID"]))
        # # print(len(ids["clone_name"]))
        # # print(len(ids["str_pos"]))
        # # print(len(ids["end_pos"]))
        # # print(len(ids["type"]))
        # ids_df = pd.DataFrame(ids)
        # print(ids_df)
        # # print(len(ids_df[ids_df["type"] == "clipped"]))
        # # print(ids_df[ids_df["str_pos"] == 777316])
        # with open("sig_read_ids.tsv", "w") as file:
        #     ids_df.to_csv(path_or_buf=file, sep="\t", index=False)

        #             for chrom in self.bins:
        #                 for line in header:
        #                     chr_name = line["SN"]
        #                     if chrom in chr_name:
        #                         chr_list.append(chrom)
        #                         for bin_str in self.bins[chrom]:
        #                             for read in bam_file.fetch(chr_name):
        #                                 # progress bar updating
        #                                 update_bar += 1
        #                                 reads_bar.update(update_bar)
        #
        #                                 # retrieve unmapped read id
        #                                 bit_flag = bin(int(read.flag))
        #                                 if bit_flag[-3] == "1":
        #                                     # -3 because in the order of bitwise FLAGS, the bit for which identify
        #                                     # the unmapped reads is the third; transforming the FLAG into binary
        #                                     # number, the order is respected starting from the wright going
        #                                     # to the left
        #                                     unmapped_ids["chr"].append(chr_name)
        #                                     unmapped_ids["ID"].append(read.query_name)
        #                                     unmapped_ids["clone_name"].append(clone)
        #
        #                                 # retrieve properly mapped reads
        #                                 elif str(read.flag) in self.flags:
        #                                     if cigar and \
        #                                             read.cigarstring is not None and \
        #                                             bin_str <= int(read.reference_start) < (bin_str + self.bin_size):
        #
        #                                         if not any(filt in read.cigarstring for filt in self.cigar_filter):
        #                                             # always the left-most position in POS field in SAM
        #                                             ids["chr"].append(chr_name)
        #                                             ids["ID"].append(read.query_name)
        #                                             ids["str_pos"].append(int(read.reference_start))
        #                                             ids["end_pos"].append(int(read.reference_start) +
        #                                                                   len(read.query_sequence) - 1)
        #                                             ids["type"].append("properly_mapped")
        #                                         else:
        #                                             ids["chr"].append(chr_name)
        #                                             ids["ID"].append(read.query_name)
        #                                             ids["str_pos"].append(int(read.reference_start))
        #                                             ids["end_pos"].append(int(read.reference_start) +
        #                                                                   len(read.query_sequence) - 1)
        #                                             ids["type"].append("clipped")
        #
        #                                     elif not cigar:
        #                                         ids["chr"].append(chrom)
        #                                         ids["ID"].append(read.query_name)
        #                                         ids["str_pos"].append(int(read.reference_start))
        #                                         ids["end_pos"].append(int(read.reference_start) + len(read.query_sequence) - 1)
        #                                 else:
        #                                     continue
        #
        #             bam_file.close()
        #
        # ids_df = pd.DataFrame(ids)
        # unmapped_ids_df = pd.DataFrame(unmapped_ids)
        # print(ids_df)
        # print(unmapped_ids_df)
        # return ids_df, unmapped_ids_df

    # def unmapped_ids(self):
    # print("unmapped_reads:", bam_file.unmapped)


#     if read.is_unmapped:
#         update_bar += 1
#         unmapped_ids["chr"].append(ref)
#         unmapped_ids["ID"].append(read.query_name)
#         unmapped_ids["clone_name"].append(clone_name)
if __name__ == "__main__":
    # bin_size = 30000
    # flags = ["99", "147", "163", "83"]
    # folder = "./"
    # bins = {"chr6": [7770000]}
    # cigar_filter = ["S", "H"]
    print("IDENTIFIER")
    ide = TestingBinReadIdentifier(bin_size=30000,
                                   flags=["99", "147", "163", "83"],
                                   folder=["./"],
                                   bins={"chr6": 7770000, "chr3": 7770000},
                                   cigar_filter=["S", "H"])

    ide.mapped_ids(cigar=True)
