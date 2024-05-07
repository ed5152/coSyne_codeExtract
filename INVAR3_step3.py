#!/usr/bin/env python
# coding: utf-8

# Now making the whole dataframe of reads


# Is mismatches per read or per fragment? FRAG. does frag database have only one F/R read?
import logging
import pysam
import pandas as pd
import numpy as np
import sys

from argparse import ArgumentParser
from pandas import DataFrame
from pandas.errors import EmptyDataError
from pathlib import Path
from pysam import AlignmentFile, FastaFile

from .pickleFunctions import *
from .trinucelotideContext import get_trinucleotide_context
from numpy import number

class INVAR3_step3:

    def __init__(self):
        self.logger = logging.getLogger('INVAR3_step3')

        self.parser = ArgumentParser(prog = 'INVAR3_step3.py',
                                     description = "INVAR3 step three.")
        self.parser.add_argument('-p', '--patient', metavar = '<id>', dest = 'patient', required = True,
                                 help = "The patient identifier.")
        self.parser.add_argument('-s', '--sample', metavar = '<id>', dest = 'sample_ID', required = True,
                                 help = "The sample identifier.")
        self.parser.add_argument('-c', '--case-or-control', metavar = '<str>', dest = 'case_or_control', required = True,
                                 help = "Whether this process is dealing with the case or a control.")
        self.parser.add_argument('-i', '--iteration', metavar = '<num>', dest = 'iteration', required = True, type = int,
                                 help = "The control iteration. Should be 0 for case.")
        self.parser.add_argument('-l', '--loci', metavar = '<file>', dest = 'loci_path', required = True, type = Path,
                                 help = "The selected loci CSV file.")
        self.parser.add_argument('-r', '--reference', metavar = '<str>', dest = 'reference_filepath', required = True, type = Path,
                                 help = "The reference genome FASTA file.")
        self.parser.add_argument('-b', '--bam', metavar = '<file>', dest = 'bam_filepath', required = True, type = Path,
                                 help = "The aligned reads in BAM format.")
        self.parser.add_argument('-o', '--out', metavar = '<file>', dest = 'outFile_path', required = True, type = Path,
                                 help = "The output file.")
        self.parser.add_argument('-e', '--export', dest = 'export', action = "store_true",
                                 help = 'Export the data frames as CSV as well.')

        self.export = False

    def parse(self, args = None):
        self.parser.parse_args(args, self)

        if not self.case_or_control in [ 'case', 'control' ]:
            raise ValueError(f"--case-or-control (-c) must be either \"case\" or \"control\", not \"{self.case_or_control}\".")

        return self


    # Might be wrong if query_position is 1-based. String index is 0-based.
    def calculate_mismatch_rate(self, fragment_sequence:str, reference_sequence:str, query_position:int) -> np.number:
        # Add check they are the same length
        if len(fragment_sequence) == len(reference_sequence):
            mismatches = 0
            # Loop and count is faster than making a list and summing.
            for i, base in enumerate(fragment_sequence):
                if i != query_position and base != reference_sequence[i]:
                    mismatches += 1
        else:
            mismatches = np.nan
        return mismatches

    def export_frame(self, df):
        try:
            export_file = self.outFile_path.parent / f"{self.outFile_path.stem}.csv"
            df.to_csv(export_file, index = False)
        except BaseException as e:
            self.logger.error(f"Error while exporting selected loci: {e}")

    def create_fragment_dataframe(self):
        bam = AlignmentFile(self.bam_filepath, "rb")

        try:
            next(bam.head(1))
        except StopIteration:
            raise EmptyDataError(f"No reads in {self.bam_filepath}.")

        # genomic positions to collect reads from
        genomic_positions = load_pickle(self.loci_path, DataFrame)

        if genomic_positions.shape[0] == 0:
            raise EmptyDataError(f"{self.loci_path} has no rows.")

        reference = FastaFile(self.reference_filepath)

        cols = ['PATIENT', 'SAMPLE_ID', 'READNAME', 'CHROM', 'POS', 'REF', 'ALT', 'READ_BASE', 'TRINUCLEOTIDE_CONTEXT', \
                'MUTANT', 'ORIG_MUT_POS',  'TUMOUR_AF', 'FRAG_LENGTH', 'BKGD_ERROR', 'READ_LEN', 'CASE_OR_CONTROL', 'ITERATION']
        analysis_df = DataFrame(columns = cols)

        for index, row in genomic_positions.iterrows():
            #sampleID = bam_file
            chromosome = row['CHROM']
            # This position refers to the now control regions of the genome
            position = row['POS']
            alt_allele = row['ALT']
            ref_allele = row['REF']
            tumour_af = row['TUMOUR_AF']
            # This position refers to the original mutation position with matched trinucleotide to the 'position' position
            origMut = chromosome + "_" + str(row['ORIG_MUT_POS'])

            for pileup_column in bam.pileup(chromosome, position-1, position, min_base_quality=20, min_mapping_quality=40, stepper="all"):
                if pileup_column.pos == position - 1:
                    for pileup_read in pileup_column.pileups:
                        if pileup_read.alignment.is_paired and not pileup_read.alignment.is_duplicate and not pileup_read.is_del \
                                and not pileup_read.is_refskip and pileup_read.alignment.is_proper_pair:
                            read = pileup_read.alignment
                            if ("D" in read.cigarstring) or ('I' in read.cigarstring):
                                #self.logger.debug("Indel found, skipping this read")
                                continue

                            try:
                                #self.logger.debug(f"{chromosome}={position}")
                                read_name = read.query_name
                                length = abs(read.template_length)
                                # For background error calc
                                read_len = abs(read.query_alignment_length)
                                start_coord = read.reference_start
                                # Position of loci of interest wrt aligned read
                                query_index_pos = position - start_coord - 1

                                refBase = read.get_reference_sequence()[query_index_pos]
                                readBase = read.query_alignment_sequence[query_index_pos]
                                #self.logger.debug(f"REF : {refBase}, READ_BASE : {readBase}")

                                # Throw Error if refBase != ref_allele
                                if refBase.upper() != ref_allele:
                                    raise ValueError(f'Indexing issue somewhere, ref base different to ref_allele in CHROM: {chromosome}, POS:{position} for file {self.loci_path}')

                            except ValueError as e:
                                self.logger.warning(f"{self.patient} {self.sample_ID}: {e}")
                                continue

                            # If no error the continue
                            # Removed this because of string issues - will add in later
                            mutant = readBase == alt_allele

                            # check the index on this too
                            mismatch_N = self.calculate_mismatch_rate(read.query_alignment_sequence, read.get_reference_sequence(), query_index_pos)
                            #context = read.query_alignment_sequence[query_index_pos - 1: query_index_pos + 2]
                            context = get_trinucleotide_context(reference, chromosome, position)
                            mergeframe = DataFrame([[self.patient, self.sample_ID, read_name, chromosome, position, ref_allele, alt_allele, \
                                                     readBase, context, mutant, origMut, tumour_af, length, mismatch_N, read_len,
                                                     self.case_or_control, self.iteration]],
                                                   columns = cols)
                            analysis_df = pd.concat([analysis_df, mergeframe], ignore_index=True)
        return analysis_df

    def main(self) -> DataFrame:
        control_df = self.create_fragment_dataframe()
        dump_pickle(self.outFile_path, control_df)

        if self.export:
            self.export_frame(control_df)

        return control_df

## Main Code :
if __name__ == "__main__":
    logging.basicConfig(level = logging.INFO)
    INVAR3_step3().parse().main()
