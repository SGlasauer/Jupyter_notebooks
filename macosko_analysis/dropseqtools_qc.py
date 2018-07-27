import pysam
from tqdm import tnrange, tqdm_notebook
from collections import Counter
import pandas as pd
import numpy as np
import pickle
import os


def get_cell_barcode_counts(bamfile, savefile):
    """
    Generates a dictionary of counts for each barcode observed in the XC tag
    :param bamfile: bamfile output from dropseqtools (star_gene_exon_tagged)
    :param savefile: full path and name of file to save results in pickle format
    :return: Dictionary of counts for each barcode detected
    """

    if os.path.exists(savefile):
        cell_barcode_counts = load_pickle(savefile)

    else:

        final_bam = pysam.AlignmentFile(bamfile, "rb")
        total_reads = int(final_bam.mapped)
        cell_barcode_counts = Counter()

        progress = tnrange(total_reads)

        for read in final_bam.fetch():
            barcode = dict(read.tags)['XC']
            cell_barcode_counts[barcode] += 1
            progress.update(1)

        save_dict_as_pickle(cell_barcode_counts, savefile)

    return cell_barcode_counts


def save_dict_as_pickle(cell_barcode_counts, savefile):
    '''
    saves the counts of cell barcodes in a pickle file
    :param cell_barcode_counts: counter output from get_cell_barcode_counts
    :param savefile: string to full path and name of picklefile to save (extension .pickle)
    :return:saves the counter as a pickle file
    '''
    with open(savefile, 'wb') as outputfile:
        pickle.dump(cell_barcode_counts, outputfile)


def load_pickle(savefile):
    """
    loads the pickle file saved from save_cell_barcode_counts
    :param savefile: full path to pickle file that was saved previously
    """
    cell_barcode_counts = pickle.load( open( savefile, "rb" ) )
    return cell_barcode_counts


def get_cell_barcodes_to_analyze(cell_barcodes_dict, num_cells=None):
    '''
    Makes a list of cell barcodes to interrogate further based on the number of reads mapping to them
    :param cell_barcodes_dict: dictionary output from get_cell_barcode_counts
    :num_cells: optional, if left empty, the top 1% of barcodes will be kept
    :return: list of cell barcodes sequences to analyze
    '''

    sorted_counts = dict(sorted(cell_barcodes_dict.items(), key=lambda x: x[1], reverse=True))

    if num_cells is None:
        num_to_keep = int(np.floor(0.01*len(sorted_counts)))

    else:
        num_to_keep = num_cells

    barcodes_to_analyze = list(sorted_counts.keys())[:num_to_keep]
    return barcodes_to_analyze


def summarize_counter_to_dataframe (counter):
    """

    :param counter: result of count_umis_per_barcode
    :return: dataframe of cell_barcode, umi, and gene_name counts
    """
    df = pd.DataFrame(list(counter.keys()))
    df_2 = pd.DataFrame(list(counter.values()))
    results = df.join(df_2, rsuffix="_count")
    old = ['0', 1, 2, '0_count']
    new = ['cell_barcode', 'umi', 'gene_name', 'count']
    renamer = dict(zip(old, new))
    results.rename(columns = renamer, inplace=True)
    results.fillna("unmapped", inplace=True)

    return results


def get_umi_counts_per_barcode(dataframe, require_mapped=False):
    """

    :param dataframe: result of summarize_counter_to_dataframe
    :return: umi_counts mapped to genes per cell barcode
    """

    if require_mapped is True:
        dataframe = dataframe.loc[dataframe['gene_name'] != "unmapped"]

    df = pd.DataFrame(dataframe.groupby(by=['cell_barcode', 'umi'])['gene_name'].count())
    cell_barcode_umi_counts = df.reset_index().groupby(by='cell_barcode')['gene_name'].sum()

    umi_counts = pd.DataFrame(cell_barcode_umi_counts)
    umi_counts.rename(columns = {'gene_name' : 'umi_count'}, inplace=True)

    return umi_counts


def count_umis_per_barcode(bamfile, barcodes, savefile):
    '''
    Make a dictionary of all unique, cell barcode, molecular barcode, genename combos and count totals
    :param bamfile:output of dropseqtools star_exon_tagged.bam
    :param barcodes:list of ~1000 barcodes to analyze
    :return:dictionary of counts for that barcode
    '''

    if os.path.isfile(savefile):
        results_counter = load_pickle(savefile)

    else:

        final_bam = pysam.AlignmentFile(bamfile, "rb")
        total_reads = int(final_bam.mapped)

        results_counter = Counter()

        progress = tnrange(total_reads)

        barcodes = set(barcodes)

        for read in final_bam.fetch():
            barcode = dict(read.tags)['XC']
            if barcode in barcodes:
                umi = dict(read.tags)['XM']
                if read.has_tag("GE"):
                    gene = dict(read.tags)['GE']
                else:
                    gene = None
                results_counter[(barcode,umi,gene)] += 1
                progress.update(1)
            else:
                progress.update(1)

        save_dict_as_pickle(results_counter, savefile)

    return results_counter


def make_df_for_kneeplot(bamfile, barcodes, savefile):
    """

    :param bamfile: output from dropseqtools
    :param barcodes: list of cell_barcodes_to_analyze
    :param savefile: pickle file to save umi and gene counts for each barcode
    :return: umi_df_for_plotting
    """

    results_counter = count_umis_per_barcode(bamfile, barcodes, savefile)

    results = summarize_counter_to_dataframe(results_counter)
    umi_df = get_umi_counts_per_barcode(results)

    umi_df.sort_values(by='umi_count', ascending=False, inplace=True)
    umi_df['percent_of_total'] = umi_df['umi_count']/umi_df['umi_count'].sum()
    umi_df['cumulative'] = np.cumsum(umi_df['percent_of_total'])

    return umi_df


def subset_cell_barcode_counts(cell_barcode_counts_dict, cell_barcodes):
    """

    :param cell_barcode_counts_dict: dictionary of cell barcodes (keys) and read count (values)
    :param cell_barcodes: list of cell barcodes to keep
    :return: subset of original dict only containing cell barcodes in the provided list
    """
    new = dict()
    for barcode in cell_barcodes:

        new[barcode] = cell_barcode_counts_dict[barcode]
    return new


def make_tidy_for_boxplot(cell_barcode_counts, all_umis, mapped_umis):
    """

    :param cell_barcode_counts: dictionary of cell barcode counts
    :param all_umis: result of get_umi_counts_per_barcode with require_mapped=False
    :param mapped_umis: result of get_umi_counts_per_barcode with require_mapped=True
    :return: dataframe in tidy format that can be used for plotting
    """

    df = pd.DataFrame.from_dict(cell_barcode_counts, orient='index')
    df.rename(columns = {0:'cell_barcode_count'}, inplace=True)
    df = df.join(all_umis)
    df.rename(columns = {'umi_count':'all_umi_count'}, inplace=True)
    df = df.join(mapped_umis)
    df.rename(columns = {'umi_count':'mapped_umi_count'}, inplace=True)

    tidy = pd.melt(df)
    tidy['log2_value'] = np.log2(tidy['value'])

    return tidy

