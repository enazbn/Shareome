import sys
import pandas as pd

input_csv = sys.argv[1]
output_parquet = sys.argv[2]

df = pd.read_csv(input_csv)

df.to_parquet(output_parquet, index=False)
