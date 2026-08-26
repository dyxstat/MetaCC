import numpy as np
import pandas as pd
import statsmodels.api as sm

def normcc(contig_file):
    # Updated column names to match the CSV headers
    names = ['contig_name', 'sites', 'length', 'covcc', 'signal']
    
    # Read the CSV, skipping the header row
    df = pd.read_csv(contig_file, header=0, names=names)
    
    # Add a small constant to avoid log(0) issues
    df['sample_site'] = np.log(df['sites'] + 1e-9)
    df['sample_len'] = np.log(df['length'] + 1e-9)
    df['sample_covcc'] = np.log(df['covcc'] + 1e-9)
    
    # Replace any infinite values with NaN
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    
    # Drop rows with NaN or zero values in the signal column
    df.dropna(subset=['sample_site', 'sample_len', 'sample_covcc', 'signal'], inplace=True)
    df = df[df['signal'] > 0]
    
    if df.empty:
        raise ValueError("Dataframe is empty after preprocessing or no valid rows with non-zero signal found")

    # Define exogenous and endogenous variables
    exog = df[['sample_site', 'sample_len', 'sample_covcc']]
    endog = df['signal']

    exog = sm.add_constant(exog)
    
    try:
        glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial(alpha=1))
        res = glm_nb.fit(method="lbfgs")
        norm_result = res.params.tolist()
        return norm_result
    except Exception as e:
        print(f"An error occurred during model fitting: {e}")
        return None
