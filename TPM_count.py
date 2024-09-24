import pandas as pd
import gffutils

# Load GTF file and create a database
gtf_file = 'Homo_sapiens.gtf'
db = gffutils.create_db(gtf_file, dbfn='annotations.db', force=True, keep_order=True, merge_strategy='merge')

# Extract gene lengths
genes = db.features_of_type('gene')
lengths = {gene.id: gene.end - gene.start + 1 for gene in genes}
lengths_df = pd.DataFrame.from_dict(lengths, orient='index', columns=['length'])

# Load HTSeq count data
counts = pd.read_csv('counts.csv', index_col=0)  # Assuming columns are named for samples

# Merge lengths with counts
counts = counts.join(lengths_df)

# Calculate TPM for each sample
tpm_results = pd.DataFrame(index=counts.index)

for sample in counts.columns:
    # Calculate RPK
    counts['RPK'] = counts[sample] / (counts['length'] / 1000)
    
    # Calculate TPM
    total_rpk = counts['RPK'].sum()
    tpm_results[sample] = (counts['RPK'] / total_rpk) * 1e6

# Save TPM results
tpm_results.to_csv('tpm_results.csv')
