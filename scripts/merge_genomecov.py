#!/usr/bin/env python

import pandas as pd
import numpy as np
from multiprocessing import Pool
from functools import partial


# FUNCTIONS
def applyParallel(grouped, func, cores):
    # Credit: https://stackoverflow.com/a/29281494/570918
    with Pool(cores) as p:
        results = p.map(func, [group for name, group in grouped])
    return pd.concat(results)


def merge_to_mode(df, epsilon):
    idx = 0
    while idx < len(df):
        cur = df.iloc[idx]
        neighborhood = df.loc[(df.position >= (cur.position - epsilon)) & (df.position <= (cur.position + epsilon)), 'reads']
        max_name = neighborhood.idxmax()
        if max_name == cur.name:
            idx += 1
        else:
            df.loc[max_name, 'reads'] = cur.reads + df.loc[max_name, 'reads']
            df.drop(cur.name, inplace=True)

    return df


# IMPORT FILE
counts = pd.read_table(snakemake.input.cov, header=None,
                       names=['seqname', 'position', 'reads'],
                       dtype={'seqname': str, 'position': np.int32, 'reads': np.int32})
print("Imported %d genomic locations" % len(counts))

merged_counts = applyParallel(counts.groupby('seqname'),
                              partial(merge_to_mode, epsilon=int(snakemake.wildcards.epsilon)),
                              int(snakemake.threads))
print("Exporting %d genomic locations" % len(merged_counts))

merged_counts.to_csv(snakemake.output.cov, sep='\t', header=False, index=False, compression='gzip')
