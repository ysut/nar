# conda create -n sss python=3.8 -y && conda activate sss
# conda install -y -c bioconda gffutils jupyter tqdm cyvcf2 pathlib2 pandarallel pysam liftover pybedtools

import os
import numpy as np
import pandas as pd
# from Bio.Seq import Seq
# from liftover import get_lifter
from pathlib2 import Path
from pandarallel import pandarallel
from tqdm import tqdm
import gffutils
import pysam
from cyvcf2 import VCF

### Logging setup
from logging import getLogger, config
import yaml
parent_directory = os.path.dirname(os.path.dirname('__file__'))
config_path: str = os.path.join(parent_directory, '../../../config/logging.yaml')
with open(config_path, 'r') as f:
    config.dictConfig(yaml.safe_load(f))
logger = getLogger(__name__)

########   Initialize and setup pandas methods   ########
os.environ['JOBLIB_TEMP_FOLDER'] = '/tmp' 
pandarallel.initialize(nb_workers=os.cpu_count()-1, progress_bar=False, verbose=2, use_memory_fs=False) 
tqdm.pandas()
import sys
sys.path.append(os.path.join(Path().resolve(), '../../../'))
# try: 
#     __file__
#     sys.path.append(os.path.join(os.path.dirname('__file__')))
# except NameError:
#     Path().resolve()
#     sys.path.append(os.path.join(Path().resolve(), '../../../'))

from libs import utils, preprocess, variantfilter, posparser, splaiparser
from libs import predeffect, scoring
from libs import anno_spliceai, anno_clinvar
# from libs.deco import print_filtering_count
# from libs import predeffect
from libs.scoring import Scoring
from libs import predeffect


gencode_gff = '../../../Resources/05_GENCODE_v43lift37/gencode.v43lift37.annotation.sort.gff3.gz'

try:
    db_anno_gencode = '../../../Resources/06_gffutilsdb/gencode.v43lift37.annotation.gtf.db'
    db_anno_intron = '../../../Resources/06_gffutilsdb/gencode.v43lift37.annotation.intron.gtf.db'
    db = gffutils.FeatureDB(db_anno_gencode)
    db_intron = gffutils.FeatureDB(db_anno_intron)
except ValueError:
    db_anno_gencode = '/resources/DBs/gencode.v43lift37.annotation.gtf.db'
    db_anno_intron = '/resources/DBs/gencode.v43lift37.annotation.intron.gtf.db'
    db = gffutils.FeatureDB(db_anno_gencode)
    db_intron = gffutils.FeatureDB(db_anno_intron)

## Thresholds configuration
thresholds_SpliceAI_parser: dict = {
    'TH_min_sALDL': 0.02, 'TH_max_sALDL': 0.2, 
    'TH_min_sAGDG': 0.01, 'TH_max_sAGDG': 0.05,
    'TH_min_GExon': 25, 'TH_max_GExon': 500,
    'TH_sAG': 0.2, 'TH_sDG': 0.2
    }

from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
import matplotlib.pyplot as plt
import seaborn as sns
import scipy

import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

try: 
    __file__
    sys.path.append(os.path.join(os.path.dirname('__file__')))
except NameError:
    Path().resolve()
    sys.path.append(os.path.join(Path().resolve(), '../../'))

from libs.scoring import Scoring

import warnings
warnings.simplefilter('ignore')


# Loading allmut variants from pickle
allmut_dm_maf0_snv_hg19 = pd.read_pickle('allmut_dm_maf0_snv_liftover.pkl')

# Rename POS_hg19 to POS
allmut_dm_maf0_snv_hg19.rename(columns={'POS_hg19': 'POS'}, inplace=True)

# Drop unknown positions in 'POS' column and assign integer type
allmut_dm_maf0_snv_hg19.dropna(subset=['POS'], inplace=True)
allmut_dm_maf0_snv_hg19 = allmut_dm_maf0_snv_hg19.astype({'POS': int})

# Change object name to allmut
allmut = allmut_dm_maf0_snv_hg19

# Generate ID column
allmut['ID'] = allmut['CHROM'].astype(str) + '-' + allmut['POS'].astype(str) + '-' + allmut['hgvs']

# Extract useful columns
allmut = allmut[['ID', 'mutype', 'clinvar_clnsig', 'tag', 'deletion', 'insertion', 'expected_inheritance', 'gnomad_AF']]

# Load VCF file annoteted by analysis pipeline

# df = pd.read_pickle('splai_vep_vcfs/hgmd_dm/allchr.DM.splai.vep.nondel.enst.prescore.hgnconly.v2.pkl')
df = pd.read_pickle("variant_data_set_vcfs/hgmd_all.prescore.onlyhgnc.pkl")

df['HGVSc'] = df['HGVSc'].str.replace('c.', '')
df['ID'] = df['CHROM'].astype(str) + '-' + df['POS'].astype(str) + '-' + df['HGVSc']

# merge df and allmut on 'ID' column with inner join
# logger.info(len(df))
df = pd.merge(df, allmut, on='ID', how='inner')
# logger.info(len(df))

exclude_csq = {
    '3_prime_UTR_variant', '5_prime_UTR_variant', 'mature_miRNA_variant',
    'mature_miRNA_variant', 'downstream_gene_variant', 'upstream_gene_variant',
    'non_coding_transcript_exon_variant'
}

exclude_non_spl_dm: set = {'splice_region_variant'}

def is_orf_variants(row):
    csqs: list = row['Consequence'].split('&')
    if set(csqs).isdisjoint(exclude_csq):
        return True
    else:
        return False
    
def is_non_spl_tn(row):
    csqs: list = row['Consequence'].split('&')
    if set(csqs).isdisjoint(exclude_non_spl_dm):
        return True
    else:
        return False

def is_gnomad_tn(row):
    if row['gnomad_AF'] == 0:
        return True
    else:
        return False


df = df[df.apply(is_orf_variants, axis=1)]

# df_spl contains splicing mutations (splice, canonical-splice, exonic-splice)
df_spl = df[df['mutype'].str.contains('splice')].copy()

# df_non_spl contains non-splicing mutations
df_non_spl = df[df['mutype'].str.contains('missense|nonsense|synonymous')].copy()
# df_non_spl = df[df['mutype'].str.contains('missense|synonymous')].copy()
# df_non_spl = df[df['mutype'].str.contains('missense|nonsense')].copy()
# df_non_spl = df[df['mutype'].str.contains('missense')].copy()
# df_non_spl = df[df['mutype'].str.contains('synonymous')].copy()

df_non_spl = df_non_spl[df_non_spl.apply(is_non_spl_tn, axis=1)]

logger.info(f"Splicing: {len(df_spl)}, Non-splicing: {len(df_non_spl)}, total: {len(df_spl) + len(df_non_spl)}")

# Annotating the label and variant_id (CHROM-POS-REF-ALT)
# When mutype is splice, the label is 1, otherwise 0
# df_gnomad = pd.read_pickle('splai_vep_vcfs/gnomadv211/allchr.gnomad.splai.vep.vcf.enst.prescore.hgnconly.v2.pkl')
df_gnomad = pd.read_pickle('variant_data_set_vcfs/gnomad_all.prescore.onlyhgnc.pkl')
df_gnomad = df_gnomad[df_gnomad.apply(is_orf_variants, axis=1)]

df_spl['LABEL'] = 1
df_non_spl['LABEL'] = 0
df_gnomad['LABEL'] = 0

df_spl['variant_id'] = df_spl['CHROM'].astype(str) + '-' + df_spl['POS'].astype(str) + '-' + df_spl['REF'] + '-' + df_spl['ALT']
df_non_spl['variant_id'] = df_non_spl['CHROM'].astype(str) + '-' + df_non_spl['POS'].astype(str) + '-' + df_non_spl['REF'] + '-' + df_non_spl['ALT']
df_gnomad['variant_id'] = df_gnomad['CHROM'].astype(str) + '-' + df_gnomad['POS'].astype(str) + '-' + df_gnomad['REF'] + '-' + df_gnomad['ALT']

# Create a dataframe tp (True Positive)
tp = df_spl.copy()

tn = df_gnomad.copy()
tn = df_non_spl.copy()


## Exclude non-ORF variants
logger.info(f"TN: {len(tn)}")
tn['is_ORF'] = tn.apply(is_orf_variants, axis=1)
tn = tn[tn['is_ORF']]
# tn = pd.read_pickle('../../TN/tn_prescore.pkl')
logger.info(f"TN: {len(tn)}")

## Summary of the dataset
logger.info(f"TP: {len(tp)}, TN: {len(tn)}")

def specificity_sensitivity_plotly(data):
    thresholds = np.arange(0, 11, 1)
    results = []

    for threshold in thresholds:
        tp = data[(data['PriorityScore'] >= threshold) & (data['LABEL'] == 1)].shape[0]
        fn = data[(data['PriorityScore'] < threshold) & (data['LABEL'] == 1)].shape[0]
        tn = data[(data['PriorityScore'] < threshold) & (data['LABEL'] == 0)].shape[0]
        fp = data[(data['PriorityScore'] >= threshold) & (data['LABEL'] == 0)].shape[0]
        specificity = tn / (tn + fp) if (tn + fp) else 0
        sensitivity = tp / (tp + fn) if (tp + fn) else 0
        # print(f"Threshold: {threshold}, TP: {tp}, FN: {fn}, TN: {tn}, FP: {fp}")
        # print(f"Threshold: {threshold}, Specificity: {specificity:.6f}, Sensitivity: {sensitivity:.6f}")
        results.append({'Threshold': threshold, 'Metric': 'Specificity', 'Value': specificity})
        results.append({'Threshold': threshold, 'Metric': 'Sensitivity', 'Value': sensitivity})

    results_df = pd.DataFrame(results)
    return results_df

def plot_sensitivity_specificity_plotly(
        results_df: pd.DataFrame, w: int, h: int):
    # Separate the dataframes for specificity and sensitivity
    specificity_df = results_df[results_df['Metric'] == 'Specificity']
    sensitivity_df = results_df[results_df['Metric'] == 'Sensitivity']

    # Plotly Graph Objectsを使用してプロット
    fig = go.Figure()

    # 特異性
    fig.add_trace(go.Scatter(
        x=specificity_df['Threshold'],
        y=specificity_df['Value'],
        marker=dict(color='#665990'),
        mode='lines+markers',
        name='Specificity',
        text=[f'Threshold: {th}, Specificity: {val:.3f}' for th, val in zip(specificity_df['Threshold'], specificity_df['Value'])],
        hoverinfo='text'
    ))
    
    # 感度
    fig.add_trace(go.Scatter(
        x=sensitivity_df['Threshold'],
        y=sensitivity_df['Value'],
        marker=dict(color='#F8ACAC'),
        mode='lines+markers',
        name='Sensitivity',
        text=[f'Threshold: {th}, Sensitivity: {val:.3f}' for th, val in zip(sensitivity_df['Threshold'], sensitivity_df['Value'])],
        hoverinfo='text'
    ))

    # Y軸のフォーマット設定
    fig.update_yaxes(tickformat=".1f")

    # グラフのレイアウト設定
    fig.update_layout(title='Sensitivity and Specificity for each threshold',
                      xaxis_title='Threshold',
                      yaxis_title='Sensitivity/Specificity',
                      plot_bgcolor='rgba(243, 243, 243, 1)',
                      paper_bgcolor='rgba(243, 243, 243, 0)',
                      font=dict(family="Arial, sans-serif", size=12, color="black"),
                      legend=dict(y=0.075, x=0.75, xanchor='right', yanchor='bottom', 
                              bgcolor='rgba(243, 243, 243, 1)',
                              font=dict(family="Arial, sans-serif", size=12, color="black"))
                              )

    # グラフサイズの調整
    fig.update_layout(width=w, height=h)
    fig.write_html("sensitivity_specificity_plot.html")

    # fig.show()
    return fig

def plot_sensitivity_specificity_plotly_without_legened(results_df):
    # Separate the dataframes for specificity and sensitivity
    specificity_df = results_df[results_df['Metric'] == 'Specificity']
    sensitivity_df = results_df[results_df['Metric'] == 'Sensitivity']

    # Plotly Graph Objectsを使用してプロット
    fig = go.Figure()

    # 特異性
    fig.add_trace(go.Scatter(
        x=specificity_df['Threshold'],
        y=specificity_df['Value'],
        marker=dict(color='green'),
        mode='lines+markers',
        name='Specificity',
        text=[f'Threshold: {th}, Specificity: {val:.8f}' for th, val in zip(specificity_df['Threshold'], specificity_df['Value'])],
        hoverinfo='text',
        showlegend=False 
    ))
    
    # 感度
    fig.add_trace(go.Scatter(
        x=sensitivity_df['Threshold'],
        y=sensitivity_df['Value'],
        marker=dict(color='orange'),
        mode='lines+markers',
        name='Sensitivity',
        text=[f'Threshold: {th}, Sensitivity: {val:.8f}' for th, val in zip(sensitivity_df['Threshold'], sensitivity_df['Value'])],
        hoverinfo='text',
        showlegend=False
    ))

    # Y軸のフォーマット設定
    fig.update_yaxes(tickformat=".1f")

    # グラフのレイアウト設定
    fig.update_layout(title='Sensitivity and Specificity for each threshold',
                      xaxis_title='Threshold',
                      yaxis_title='Sensitivity/Specificity',
                      plot_bgcolor='rgba(243, 243, 243, 1)',
                      paper_bgcolor='rgba(243, 243, 243, 0)',
                      font=dict(family="Arial, sans-serif", size=12, color="black"),
                      legend=dict(y=0.075, x=0.75, xanchor='right', yanchor='bottom', 
                              bgcolor='rgba(243, 243, 243, 1)',
                              font=dict(family="Arial, sans-serif", size=12, color="black"))
                              )

    # グラフサイズの調整
    fig.update_layout(width=600, height=600)
    fig.write_html("sensitivity_specificity_plot.html")

    # fig.show()
    return fig


# Code below is adapted from Netflix's VMAF and BesenbacherLab's ROC-utils
# https://github.com/Netflix/vmaf/
# https://github.com/BesenbacherLab/ROC-utils
# Modifications: np.float -> np.float64

def compute_midrank(x):
    """Computes midranks.
    Args:
       x - a 1D numpy array
    Returns:
       array of midranks
    """
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N, dtype=np.float64)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5*(i + j - 1)
        i = j
    T2 = np.empty(N, dtype=np.float64)
    # Note(kazeevn) +1 is due to Python using 0-based indexing
    # instead of 1-based in the AUC formula in the paper
    T2[J] = T + 1
    return T2


def fastDeLong(predictions_sorted_transposed, label_1_count):
    """
    The fast version of DeLong's method for computing the covariance of
    unadjusted AUC.
    Args:
       predictions_sorted_transposed: a 2D numpy.array[n_classifiers, n_examples]
          sorted such as the examples with label "1" are first
    Returns:
       (AUC value, DeLong covariance)
    Reference:
     @article{sun2014fast,
       title={Fast Implementation of DeLong's Algorithm for
              Comparing the Areas Under Correlated Receiver Operating Characteristic Curves},
       author={Xu Sun and Weichao Xu},
       journal={IEEE Signal Processing Letters},
       volume={21},
       number={11},
       pages={1389--1393},
       year={2014},
       publisher={IEEE}
     }
    """
    # Short variables are named as they are in the paper
    m = label_1_count
    n = predictions_sorted_transposed.shape[1] - m
    positive_examples = predictions_sorted_transposed[:, :m]
    negative_examples = predictions_sorted_transposed[:, m:]
    k = predictions_sorted_transposed.shape[0]

    tx = np.empty([k, m], dtype=np.float64)
    ty = np.empty([k, n], dtype=np.float64)
    tz = np.empty([k, m + n], dtype=np.float64)
    for r in range(k):
        tx[r, :] = compute_midrank(positive_examples[r, :])
        ty[r, :] = compute_midrank(negative_examples[r, :])
        tz[r, :] = compute_midrank(predictions_sorted_transposed[r, :])
    aucs = tz[:, :m].sum(axis=1) / m / n - float(m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx[:, :]) / n
    v10 = 1.0 - (tz[:, m:] - ty[:, :]) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    delongcov = sx / m + sy / n
    return aucs, delongcov


def calc_pvalue(aucs, sigma):
    """Computes log(10) of p-values.
    Args:
       aucs: 1D array of AUCs
       sigma: AUC DeLong covariances
    Returns:
       log10(pvalue)
    """
    l = np.array([[1, -1]])
    z = np.abs(np.diff(aucs)) / np.sqrt(np.dot(np.dot(l, sigma), l.T))
    return np.log10(2) + scipy.stats.norm.logsf(z, loc=0, scale=1) / np.log(10)


def compute_ground_truth_statistics(ground_truth):
    assert np.array_equal(np.unique(ground_truth), [0, 1])
    order = (-ground_truth).argsort()
    label_1_count = int(ground_truth.sum())
    return order, label_1_count


def delong_roc_variance(ground_truth, predictions):
    """
    Computes ROC AUC variance for a single set of predictions
    Args:
       ground_truth: np.array of 0 and 1
       predictions: np.array of floats of the probability of being class 1
    """
    order, label_1_count = compute_ground_truth_statistics(ground_truth)
    predictions_sorted_transposed = predictions[np.newaxis, order]
    aucs, delongcov = fastDeLong(predictions_sorted_transposed, label_1_count)
    assert len(aucs) == 1, "There is a bug in the code, please forward this to the developers"
    return aucs[0], delongcov


def delong_roc_test(ground_truth, predictions_one, predictions_two):
    """
    Computes log(p-value) for hypothesis that two ROC AUCs are different
    Args:
       ground_truth: np.array of 0 and 1
       predictions_one: predictions of the first model,
          np.array of floats of the probability of being class 1
       predictions_two: predictions of the second model,
          np.array of floats of the probability of being class 1
    """
    order, label_1_count = compute_ground_truth_statistics(ground_truth)
    predictions_sorted_transposed = np.vstack((predictions_one, predictions_two))[:, order]
    aucs, delongcov = fastDeLong(predictions_sorted_transposed, label_1_count)
    return calc_pvalue(aucs, delongcov)

# Calculate AUC confidence interval (95%)
def compute_auc_confidence_interval(auc, var, confidence_level=0.95):
    alpha = 1 - confidence_level
    z_score = scipy.stats.norm.ppf(1 - alpha/2)  # 2-tailed z score
    se = np.sqrt(var)  # Calculate SE from variance
    lower_bound = auc - z_score * se
    upper_bound = auc + z_score * se
    return lower_bound, upper_bound


import pickle
all_solutions = pickle.load(open('all_solutions.pkl', 'rb'))

frac = 0.8
random_state = 13

results = []
buf: float = 0.980

for i, solution in enumerate(all_solutions):
    # if i < 2000:
    #     continue

    ths_scores = {'clinvar_same_pos': solution['s1'],
            'clinvar_same_motif': solution['s2'],
            'clinvar_else': solution['s3'],
            'non_canon_splai_lte_0.1_outside': solution['s4'],    
            'non_canon_splai_lte_0.1_other': solution['s5'],
            'non_canon_splai_bet_0.1_0.2': solution['s6'],
            'non_canon_splai_gte_0.2': solution['s7'],
            'canon_strong': solution['s8'], 
            'canon_moderate': solution['s9'], 
            'frameshift_nmd_eloF': solution['s10'], 
            'frameshift_nmd_not_eloF': solution['s11'],
            'canon_splai_lte_0.1': solution['s12'],
            'canon_splai_bet_0.1_0.2': solution['s13'],
            'canon_splai_gte_0.2': solution['s14'],
            'clinvar_blb': solution['s15']
            }
    
    scoring = Scoring(ths=ths_scores)

    tp_train = pd.read_pickle(f'train_test_pkls/tp_prescore_train_{random_state}.pkl')
    tn_train = pd.read_pickle(f'train_test_pkls/tn_prescore_train_{random_state}.pkl')

    tp_train['insilico_screening'] = tp_train.parallel_apply(scoring.insilico_screening, axis=1)
    tp_train['clinvar_screening'] = tp_train.parallel_apply(scoring.clinvar_screening, axis=1)
    # tp_train['PriorityScore'] = tp_train.parallel_apply(scoring.calc_priority_score, axis=1)
    # tp_train = scoring.calc_priority_score(tp_train)
    tp_train = tp_train[tp_train['insilico_screening'] != 'Not available']
    tp_train['PriorityScore'] = tp_train['insilico_screening'] + tp_train['clinvar_screening']
    

    tn_train['insilico_screening'] = tn_train.parallel_apply(scoring.insilico_screening, axis=1)
    tn_train['clinvar_screening'] = tn_train.parallel_apply(scoring.clinvar_screening, axis=1)
    # tn_train['PriorityScore'] = tn_train.parallel_apply(scoring.calc_priority_score, axis=1)
    # tn_train = scoring.calc_priority_score(tn_train)
    tn_train = tn_train[tn_train['insilico_screening'] != 'Not available']
    tn_train['PriorityScore'] = tn_train['insilico_screening'] + tn_train['clinvar_screening']

    # Extract the columns needed
    tp_train = tp_train[['variant_id', 'LABEL', 'PriorityScore', 'maxsplai']]
    tn_train = tn_train[['variant_id', 'LABEL', 'PriorityScore', 'maxsplai']]

    ### ========================================================== ##
    data = pd.concat([tp_train, tn_train], ignore_index=True)
    data.drop_duplicates(subset='variant_id', keep=False, inplace=True)

    # Cast the columns to float type
    data['LABEL'] = data['LABEL'].astype(int)
    # Extract rows with PriorityScore not 'Not available'
    data = data[data['PriorityScore'] != 'Not available']
    data['PriorityScore'] = data['PriorityScore'].astype(float)
    data['maxsplai'] = data['maxsplai'].astype(float)

    ## DeLong test and AUC confidence interval
    ground_truth = np.array(data['LABEL'])
    predictions_fw = np.array(data['PriorityScore'])

    auc1, var1 = delong_roc_variance(ground_truth, predictions_fw)
    cilower1, ciupper1 = compute_auc_confidence_interval(auc1, var1)

    results.append(
        {'index': i+1, 's1': solution['s1'], 's2': solution['s2'], 
         's3': solution['s3'], 's4': solution['s4'], 's5': solution['s5'], 
         's6': solution['s6'], 's7': solution['s7'], 's8': solution['s8'], 
         's9': solution['s9'], 's10': solution['s10'], 's11': solution['s11'], 
         's12': solution['s12'], 's13': solution['s13'], 's14': solution['s14'],
         'auROC': f"{auc1:.10f}, '95% Confidence Interval': {cilower1:.12f}-{ciupper1:.12f}"
        }
    )
    if i % 50 == 0:
        logger.info(f"###  Processed {i} solutions  ###")
    logger.info(f"Processed solution {i+1}: AUC: {auc1:.10f}")
    
    if auc1 > buf:
        buf = auc1
        logger.info(f"\n===== New best AUC: {auc1:.10f} with solution {i+1} =======")
        logger.info(f"New best solution {i}: {solution} \n")
        predictions_sp = np.array(data['maxsplai'])
        auc2, var2 = delong_roc_variance(ground_truth, predictions_sp)
        cilower2, ciupper2 = compute_auc_confidence_interval(auc2, var2)
        p_value_log = delong_roc_test(ground_truth, predictions_fw, predictions_sp)
        logger.info(f"AUC - Framework (95%CI): {auc1:.3f} [{cilower1:.4f}-{ciupper1:.4f}]")
        logger.info(f"AUC - SpliceAI (95%CI) : {auc2:.3f} [{cilower2:.4f}-{ciupper2:.4f}]")
        logger.info(f"p-value (DeLong Test)  : {10**p_value_log[0][0]:.2e}\n")
        logger.info("===========================================================")
        
