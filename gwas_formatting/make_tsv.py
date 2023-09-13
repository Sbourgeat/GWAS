import sys
import pandas as pd

def modify_dataframe(df):
    """
    Modifies the dataframe by adding a prefix to the 'DGRP' column and replacing values in the 'sex' column.
    
    Returns:
    df (pandas.DataFrame): The modified dataframe.
    """
    # Modify the column DGRP
    df['DGRP'] = df['DGRP'].apply(lambda x: 'DGRP_0' + x if len(x) == 2 else 'DGRP_' + x)
    # Replace values in the column sex
    df['sex'] = df['sex'].replace({'male': 'M', 'female': 'F'})

    return df
 
# Read the input file
print(sys.argv[1],"reading...")
df = pd.read_csv(sys.argv[1])
# Modify the dataframe
df = modify_dataframe(df)

# Save the modified dataframe as a new file with the same name as the input file but with .tsv at the end
output_file = sys.argv[1][:-4] + '.tsv'
df.to_csv(output_file, sep='\t', index=False)
