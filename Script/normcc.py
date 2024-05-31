import numpy as np
import pandas as pd
import statsmodels.api as sm

def normcc(contig_file):
    names = ['contig_name', 'site', 'length', 'covcc', 'signal']
    df = pd.read_csv(contig_file, header=None, names=names)
    df['sample_site'] = np.log(df['site'])
    df['sample_len'] = np.log(df['length'])
    df['sample_covcc'] = np.log(df['covcc'])
    exog = df[['sample_site', 'sample_len', 'sample_covcc']]
    endog = df[['signal']]
    exog = sm.add_constant(exog)
    glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial(alpha=1))
    res = glm_nb.fit(method="lbfgs")
    norm_result = res.params.to_list()
    return norm_result


def normcc_LC(contig_file):
    names = ['contig_name', 'length', 'covcc', 'signal']
    df = pd.read_csv(contig_file, header=None, names=names)
    df['sample_len'] = np.log(df['length'])
    df['sample_covcc'] = np.log(df['covcc'])
    exog = df[['sample_len', 'sample_covcc']]
    endog = df[['signal']]
    exog = sm.add_constant(exog)
    glm_nb = sm.GLM(endog, exog, family=sm.families.NegativeBinomial(alpha=1))
    res = glm_nb.fit(method="lbfgs")
    norm_result = res.params.to_list()
    return norm_result