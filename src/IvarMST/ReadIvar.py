"""
Read all Ivar data and generate mutation keys for matrix creation.


Really tired when writing this and waiting to go home...

2022-06-09: Matthew Wells
"""

from dataclasses import dataclass
import os
import glob
from random import sample
import pandas as pd
import numpy as np
import sys
import datetime

class IvarMatrix:
    """
    Read ivar files into the program, generate the matrix of sample presence or abscence
    Header: Mutations
    Rows: samples 1 or 0
    """

    mutations = set()
    samples = set()

    def __init__(self, glob_dir, outfile) -> None:
        self.glob_dir = glob_dir
        self.outfile = outfile
        self.find_and_read_files()
    
    def find_and_read_files(self):
        """
        Take in all of the ivar files read their data in,
        need to populate a dicitionary of mutations found as well.

        Doing two reads of files because I am tired and this can be cached later if it is insightful
        """
        st = datetime.datetime.now()
        files = glob.glob(os.path.join(self.glob_dir, "*.tsv")) # ivar spits out a tsv
        for i in files:
            fp = os.path.join(self.glob_dir, i)
            if os.path.isfile(fp):
                self.samples.add(self.get_sample_name(fp))
                with open(fp, 'r') as vcf_file:
                    lines = vcf_file.readlines()
                    for line in lines[1:]: # 1 = pos, 2 = ref, 3 = alt
                        self.update_mutation_set(line)
        df = self.create_matrix_muts()
        for i in files:
            fp = os.path.join(self.glob_dir, i)
            sample_name = self.get_sample_name(fp)
            with open(fp, 'r') as vcf_file:
                lines = vcf_file.readlines() # so slow i am sorry
                for line in lines[1:]:
                    mut = self.prep_mutations(line)
                    df.at[sample_name, mut] = 2
        df1 = df.sum(axis=0)
        percent_sum = int(len(self.samples) * 1.05)
        df1 = df1.where(df1 > percent_sum)
        df = df[[i for i in df1.index if df1[i] > percent_sum]]
        #df = df.loc[df["sum"] >= percent_sum]
        df.to_csv(self.outfile)
        end = datetime.datetime.now()
        print(f"Program took {end - st} seconds")
    
    def create_matrix_muts(self):
        """
        """
        df = pd.DataFrame(1, index=self.samples, columns=self.mutations)
        return df

    def get_sample_name(self, i):
        """
        get the sample name of the file
        """
        fn = os.path.basename(i)
        return fn[:fn.index(".")]

    def update_mutation_set(self, line):
        """
        Update the mutations set
        """
        self.mutations.add(self.prep_mutations(line))
    
    def prep_mutations(self, line_of_vcf):
        """
        Add mutations to a set
        """
        split_line = line_of_vcf.split("\t")
        mut = split_line[1] + split_line[3]
        return mut




if __name__=="__main__":
    IvarMatrix(sys.argv[1], sys.argv[2])
