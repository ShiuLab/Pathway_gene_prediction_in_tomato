import os,sys
import pandas as pd
import numpy
from numpy import random
from numpy.random import shuffle
file = sys.argv[1]
df = pd.read_csv(file, sep='\t', index_col = 0, header = 0)
# rowname = df.index.tolist()
# row = shuffle(rowname)
# df.index = rowname
shuffle(df.Class)
df.to_csv(file + '_random', index=True, header=True,sep="\t")
