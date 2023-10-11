import sys
import pandas as pd 


df1 = pd.read_csv(sys.argv[1], sep=',')
df2 = pd.read_csv(sys.argv[2], sep=',')

keys = list(set(df1.columns).intersection(set(df2.columns)))
df = pd.merge(df2, df1, on=keys, how='inner')
df.to_csv(sys.stdout, sep=',', index=False)
