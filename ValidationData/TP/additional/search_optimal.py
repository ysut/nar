
import os
import numpy as np
import pandas as pd
from pathlib2 import Path
from pandarallel import pandarallel
from tqdm import tqdm
import gffutils
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
# pandarallel.initialize(nb_workers=os.cpu_count()-8, progress_bar=False, verbose=2, use_memory_fs=False) 
pandarallel.initialize(nb_workers=5, progress_bar=False, verbose=2, use_memory_fs=False) 

tqdm.pandas()

import sys
sys.path.append(os.path.join(Path().resolve(), '../../../'))


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


# 関数の設定
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

def map_and_calc_score(row, score_map: dict) -> int:
    """
    Map the score to the solution
    s1, s2, s3, and s15 are clinvar_screening
    s4, s5, s6, s7, s8, s9, s10 and s11 are insilico_screening
    s12, s13 and s14 are recalibrated_splai
    PriortiyScore is the sum of the "clinvar_screening", "insilico_screening", and "recalibrated_splai"
    """
    return int(score_map[row['recalibrated_splai']]) + int(score_map[row['insilico_screening']]) + int(score_map[row['clinvar_screening']])


def search_optimal(all_solutions, start, tp_train, tn_train, auc2):
    results = []
    buf: float = auc2
    max_auc: float = 0

    for i, solution in enumerate(all_solutions):
        if i <= start:
            continue

        solution.update({"s0": 0})

        tp_train['PriorityScore'] = tp_train.parallel_apply(
            map_and_calc_score, args=(solution,), axis=1)
        tn_train['PriorityScore'] = tn_train.parallel_apply(
            map_and_calc_score, args=(solution,), axis=1)

        data = pd.concat([tp_train, tn_train], ignore_index=True)
        data.drop_duplicates(subset='variant_id', keep=False, inplace=True)
        data = data[data['PriorityScore'] != 'Not available']
        data['LABEL'] = data['LABEL'].astype(int)
        data['PriorityScore'] = data['PriorityScore'].astype(float)
        data['maxsplai'] = data['maxsplai'].astype(float)
        ground_truth = np.array(data['LABEL'])

        predictions_fw = np.array(data['PriorityScore'])
        auc1, var1 = delong_roc_variance(ground_truth, predictions_fw)
        cilower1, ciupper1 = compute_auc_confidence_interval(auc1, var1)

        results.append({'index': i+1, 'solution': solution, 
                        'AUC': auc1, 'CI_lower': cilower1, 'CI_upper': ciupper1})

        max_auc = max(max_auc, auc1)
        
        if auc1 > buf:
            buf = auc1
            print(f"\n===== New best AUC: {auc1:.10f} with solution {i+1} =======")
            print(f"New best solution {i}: {solution} \n")
            predictions_sp = np.array(data['maxsplai'])
            auc2, var2 = delong_roc_variance(ground_truth, predictions_sp)
            cilower2, ciupper2 = compute_auc_confidence_interval(auc2, var2)
            p_value_log = delong_roc_test(ground_truth, predictions_fw, predictions_sp)
            print(f"AUC - Framework (95%CI): {auc1:.3f} [{cilower1:.4f}-{ciupper1:.4f}]")
            print(f"AUC - SpliceAI (95%CI) : {auc2:.3f} [{cilower2:.4f}-{ciupper2:.4f}]")
            print(f"p-value (DeLong Test)  : {10**p_value_log[0][0]:.2e}\n")

        if i % 10 == 0:
            print(f"###  Processed {i+1} solutions. Max AUC: {max_auc:.8f}  ###")

        if i > start + 4051:
            break
    return results


if __name__ == '__main__':
    start: int = int(sys.argv[1])
    
    auc2 = 0.983
    random_state = 12
    all_solutions_pkl = 'all_solutions_-5-9.pkl'
    
    # Load the results of ortools
    import pickle
    all_solutions = pickle.load(open(all_solutions_pkl, 'rb'))

    # Load train data from pickle
    tp_train = pd.read_pickle(f'train_test_pkls/tp_prescore_train_{random_state}_scored.pkl')
    tn_train = pd.read_pickle(f'train_test_pkls/tn_prescore_train_{random_state}_scored.pkl')

    results = search_optimal(all_solutions, start=start, tp_train=tp_train, tn_train=tn_train, auc2=auc2)
    # Save the results with start index
    results_pkl = f'results/results_{start}_to_end.pkl'
    with open(results_pkl, 'wb') as f:
        pickle.dump(results, f)
