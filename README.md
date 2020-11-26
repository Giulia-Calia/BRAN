# BRAN: a novel tool to explore subclonal genetic variants

BRAN is a Python implemented software that paves the way to the investigation of possible subclonal genetic variants through Illumina whole-genome sequencing data coming from _in vitro_ regenerated plants.
The development of this tool was guided by the previous work on _Solanum tuberosum_ conducted by [Fossi _et.al_, 2019] [1], in which the genomic instability caused by plant regeneration techniques was investigated.  


[1]: Fossi, M., Amundson, K., Kuppu, S., Britt, A. & Comai, L. Regeneration of Solanum tuberosum Plants from Protoplasts Induces Widespread Genome Instability. Plant Physiology 180, 78–86 (2019).

## Implementation 
BRAN is constituted by 3 interrelated classes all converged into a 4th one, the BinReadAnalyzer. This classes perform a chromosome binning parsing BAM formatted files and dividing each chromosome into equal, consecutive and non overlapping segment of a choosen size (bins). Subsequently BRAN counts how many reads (properly filtered), fall into every single bin. After the counting process, BRAN normalizes the counts and calculates a log2 fold-change, to try to capture the differences between control plant(s) and sample plant(s) coming form the rigeneration process. Finally it displays different informative plots as interactive HTML pages and in the meanwhile it locally saves them in a static format. 
Additionally it outputs two files containing information about all bins resulted respectivelly significantly different and not significantly different from the control plant(s), based on the log2 fold-change calculation. 

## Getting Started 

BRAN has many parameters that can be setted, depending on the scope of it usage, and this makes it more user-friendly than other complementary tools.

```
python3 /path_to_BRAN/BinReadAnalyzer.py -h


  -bs BIN_SIZE, --bin_size BIN_SIZE
                        The length in bp, of the segments that split the
                        chromosomes equally.
  -fl FLAG_LIST [FLAG_LIST ...], --flag_list FLAG_LIST [FLAG_LIST ...]
                        A list of the bitwise-flags in BAM/SAM format that
                        identify the reads to be counted during analyses; if
                        different flags wants to be added, add them as strings
                        (e.g. "177" "129").
  -fc FOLD_CHANGE, --fold_change FOLD_CHANGE
                        A number to set as fold-change cut-off to detect 
                        significantly different bins. 
  -co CONTROL_NAME, --control_name CONTROL_NAME
                        The name of the control sample(s) for the fold_change
                        analysis, it should be the same name as the column
                        name in the read_counts data structure, so the name of
                        the alignment file used as a baseline, without the
                        '.bam'.
  -pw, --pairwise       If specified, the fold change is calculated pairwise
                        between each sample and the control_sample.
  -c, --cigar           If specified, it allows the application of filters on
                        cigar_string, per read.
  -cf CIGAR_FILTER [CIGAR_FILTER ...], --cigar_filter CIGAR_FILTER [CIGAR_FILTER ...]
                        If specified, the reads mapped with soft and hard
                        clipping (S and H) by default, are taken out form the
                        read counts; it returns a data frame with same
                        structure of the default one but with 2 different
                        columns for sample, one without counting the reads
                        having these filters and the other counting only the
                        reads having these filters. (Specify other filters
                        using e.g. "I" "D").
  -u, --unmapped        If specified, also a .txt file is created, with all
                        the unmapped reads and, as last raw, the counts for
                        each sample.
  -ch CHROMOSOME, --chromosome CHROMOSOME
                        The number of the chromosome of interest to obtain a
                        restricted view of the data.
  -s SAMPLE, --sample SAMPLE
                        The name of the clone of interest to obtain a
                        restricted view of the data.
  -bc BIN_CHROMOSOMES [BIN_CHROMOSOMES ...], --bin_chromosomes BIN_CHROMOSOMES [BIN_CHROMOSOMES ...]
                        The name of the chromosome for each interesting bin
                        (no repetitions) from which retrieve IDs information. 
                        Specify it as e.g. -bc chr1 chr2 chr4.
  -bp BIN_POSITIONS, --bin_positions BIN_POSITIONS
                        The bin position on the corresponding chromosome, be
                        careful that for each position/list of positions there
                        is one and only one chromosome. 
                        Specify it like e.g. -bp 1000000 -bp 1200000 -bp [1340000, 6260000].
  -id, --identifier     If specified also the retrieval of read IDs will be done.
  -f FOLDER [FOLDER ...], --folder FOLDER [FOLDER ...]
                        The path to the folder in which are located the files
                        to be analyzed (.bam).
  -op OUTPUT_PICKLE, --output_pickle OUTPUT_PICKLE
                        The path to the folder in which search the pickle file
                        already created.
  -sf SAVING_FOLDER, --saving_folder SAVING_FOLDER
                        Path to the directory in which save all outputs.
                        ATTENTION: for plots, new directory is automatically
                        created with default name: /plots/[bin_size]/.
  -fo SAVING_FORMAT, --saving_format SAVING_FORMAT
                        file format for saved plot images, the choice is
                        between: ["svg", "jpeg", "pdf", "png"]. Default format
                        is: svg.
  -vb, --violin_bar     If specified BRAN returns and save only the violin and
                        the bar plots.
  -sc, --scatter        If specified BRAN returns and save only the scatter
                        distribution plots.
  -cp, --fold_change_pl
                        If specified BRAN returns and save only the fold change
                        plots.
  -i, --read_info       If specified, a data-frame with information on the
                        read ID, and if the read and its mate map in the same
                        bin in the same chromosome, is created.
```
