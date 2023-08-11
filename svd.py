from scipy.io import loadmat
import numpy as np
import pickle as pkl
import pandas as pd
import time
import subprocess
import sys
#
df = pd.read_csv( 'data.csv', sep=' ', header=None)
U, S, Vh = np.linalg.svd(df.values, full_matrices=True)
print(U.shape)
print(S.shape)
print(Vh.shape)
#
df_U=pd.DataFrame(U)
df_U.to_csv('U.cvs',sep=' ',header=False,index=False)
#
df_S=pd.DataFrame(S)
df_S.to_csv('S.cvs',sep=' ',header=False,index=False)
#
df_Vh=pd.DataFrame(Vh)
df_Vh.to_csv('Vh.cvs',sep=' ',header=False,index=False)
