import pandas as pd
import glob
import sys

output = sys.argv[1]

files = glob.glob('results/*.g.vcf')
files = sorted(files)
names = [x[8:14] for x in files]
d = {'col1': names, 'col2': files}
df = pd.DataFrame(d)
df.to_csv(output, sep='\t', header=False, index=False)
