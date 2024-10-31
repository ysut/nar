# conda create -n sss python=3.8 -y && conda activate sss
# conda install -y -c bioconda gffutils jupyter tqdm cyvcf2 pathlib2 pandarallel pysam liftover pybedtools

import os
import re
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
pandarallel.initialize(nb_workers=5, progress_bar=True, verbose=1, use_memory_fs=False) 
tqdm.pandas()

import sys
try: 
    __file__
    sys.path.append(os.path.join(os.path.dirname('__file__')))
except NameError:
    Path().resolve()
    sys.path.append(os.path.join(Path().resolve(), '../../../'))

from libs import utils, preprocess, variantfilter, posparser, splaiparser
# from libs import predeffect, scoring
from libs import anno_spliceai, anno_clinvar
from libs.deco import print_filtering_count
# from libs import predeffect
from libs.scoring import Scoring

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

## Parse VCF to simple input table
chr: str = "22"
raw_vcf: str = f"splai_vep_vcfs/hgmd_dm/all_DM_chr{chr}.splai.vep.nondel.vcf"

vcf = VCF(raw_vcf)
header = vcf.header_iter()
for h in header:
	try:
		h['ID']
	except KeyError:
		continue
	else:
		if h['ID'] == 'CSQ':
			vep_cols_list = h['Description'].split('Format: ')[1].rstrip('"').split('|')
		elif h['ID'] == 'SpliceAI':
			splai_cols_list = h['Description'].split('Format: ')[1].rstrip('"').split('|')
		else:
			pass

vepidx: dict = {col: i for i, col in enumerate(vep_cols_list)}
splaidx: dict = {col: i for i, col in enumerate(splai_cols_list)}

cols = [
	'CHROM', 'POS', 'REF', 'ALT', 'GeneSymbol', 'SymbolSource', 'HGNC_ID', 
	'ENST', 'HGVSc', 'Consequence', 'EXON', 'INTRON', 'Strand',
	'DS_AG', 'DS_AL', 'DS_DG', 'DS_DL', 'DP_AG', 'DP_AL', 'DP_DG', 'DP_DL', 'MaxSpliceAI'
]

# print(vepidx)

df: pd.DataFrame = pd.DataFrame(columns=cols)
for v in VCF(raw_vcf):
	vep: list = v.INFO.get('CSQ').split('|')

	# Get SpliceAI scores
	if v.INFO.get('SpliceAI'):
		splai: list = v.INFO.get('SpliceAI').split(',')[0].split('|')
	else:
		splai = ['NA'] * len(splai_cols_list)

	# Get HGVSc from VEP
	try:
		hgvsc = re.search('(?<=:).*',vep[vepidx['HGVSc']])[0]
	except TypeError:
		hgvsc = "NA"

	# Convert strand to +/- 
	strand = lambda s: '+' if s == '1' else '-'

	# Get max SpliceAI scores
	ds_ag: float = splai[splaidx['DS_AG']]
	ds_al: float = splai[splaidx['DS_AL']]
	ds_dg: float = splai[splaidx['DS_DG']]
	ds_dl: float = splai[splaidx['DS_DL']]
	if splai[splaidx['DP_AG']] == 'NA':
		maxsplai: str = "NA"
	maxsplai: float = max(ds_ag, ds_al, ds_dg, ds_dl)

	# Add df row
	df = pd.concat([df, pd.DataFrame([[
		v.CHROM, v.POS, v.REF, v.ALT[0], 
		vep[vepidx['SYMBOL']], vep[vepidx['SYMBOL_SOURCE']], vep[vepidx['HGNC_ID']], 
		vep[vepidx['Feature']], hgvsc, vep[vepidx['Consequence']], 
		vep[vepidx['EXON']], vep[vepidx['INTRON']],
		strand(vep[vepidx['STRAND']]), 
		ds_ag, ds_al, ds_dg, ds_dl,
		splai[splaidx['DP_AG']], splai[splaidx['DP_AL']], 
		splai[splaidx['DP_DG']], splai[splaidx['DP_DL']],
		maxsplai
	]], columns=cols)], ignore_index=True)

	# if hgvsc == "NA":
	#     logger.warning(f"[{v.CHROM}:{v.POS}] HGVSc not found")
	# if maxsplai == "NA":
	#     logger.warning(f"[{v.CHROM}:{v.POS}] SpliceAI scores not found")

# ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
# CHROM, POS, REF, ALT, GeneSymbol, NCBI_ID, ENST, ExonIntronNumbers, FLAGS, SYMBOL_SOURCE|


#### Very slow process ####
# Annotate ENST Full ID for fetching variant information from GENCODE database

chr_list = [str(i) for i in range(1, 23)] + ['X', 'Y']

for chr in chr_list:
	logger.info(f"chr{chr}")
	df = pd.read_pickle(f'splai_vep_vcfs/hgmd_dm/all_DM_chr{chr}.splai.vep.nondel.vcf.pkl')
	print(len(df))
	df.drop_duplicates(inplace=True)
	df['ENST_Full'] = df.progress_apply(posparser.fetch_enst_full, db=db, axis=1)
	df.to_pickle(f'splai_vep_vcfs/hgmd_dm/all_DM_chr{chr}.splai.vep.nondel.vcf.enst.pkl')
	

for chr in chr_list:
    df = pd.read_pickle(f'splai_vep_vcfs/hgmd_dm/all_DM_chr{chr}.splai.vep.nondel.vcf.enst.pkl')

    logger.info('Classify "Canonical" splice site or "Non-canonical" splice site...')
    df = posparser.classifying_canonical(df)

    logger.info('Calculate the distance to the nearest splice site in intron variant...')
    df['IntronDist'] = df.progress_apply(
        posparser.signed_distance_to_exon_boundary, 
        db=db, db_intron=db_intron, axis=1)

    tbx_anno = pysam.TabixFile(gencode_gff)
    df['exon_loc'] = df.progress_apply(
        posparser.calc_exon_loc, tabixfile=tbx_anno, enstcolname='ENST', axis=1)
    df = pd.concat([df, df['exon_loc'].str.split(':', expand=True)], axis=1)
    df.rename(columns={0: 'ex_up_dist', 1: 'ex_down_dist'}, inplace=True)
    df.drop(columns=['exon_loc'], inplace=True)

    #2-2. Select minimum distance from upstream distance and downstream distance
    df['exon_pos'] = df.parallel_apply(posparser.select_exon_pos, axis=1)
    #2-3. Relative exon location
    df['prc_exon_loc'] = df.parallel_apply(posparser.calc_prc_exon_loc, axis=1)

    #2-4. Decision exonic splice sites (1 nt in acceptor site or 3 nts on Donor site)
    df['exon_splice_site'] = df.parallel_apply(posparser.extract_splicing_region, axis=1)

    #3.   Additional Splicing information
    logger.info('Annotating splicing information...')
    #3-1. Annotate splicing type ('Exonic Acceptor' etc.)
    df['SpliceType'] = df.parallel_apply(posparser.select_donor_acceptor, axis=1)

    #5.   Annotate ClinVar varaints interpretations
    logger.info('Annotating ClinVar varaints interpretations...')
    clinvar_file = '../../../Resources/03_ClinVar/variant_summary.snv.grch37.germline.criteria.sort.bed.gz'
    tbx_clinvar = pysam.TabixFile(clinvar_file)
    df['clinvar_same_pos'] = df.progress_apply(
        anno_clinvar.anno_same_pos_vars, tabixfile=tbx_clinvar, axis=1)
    df['clinvar_same_motif'] = df.progress_apply(
        anno_clinvar.anno_same_motif_vars, tabixfile=tbx_clinvar, axis=1)

    logger.info('Parsing SpliceAI results...')
    logger.info('Annotating Exon/Intron position information...')
    df['ExInt_INFO'] = df.progress_apply(
        splaiparser.calc_exint_info, db=db, db_intron=db_intron, axis=1)

    #6-3. Predict splicing effects
    df['Pseudoexon'] = df.progress_apply(
        splaiparser.pseudoexon_activation,
        thresholds=thresholds_SpliceAI_parser, 
        db_intron=db_intron,
        axis=1)

    df['Part_IntRet'] = df.parallel_apply(
        splaiparser.partial_intron_retention,
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)

    df['Part_ExDel'] = df.parallel_apply(
        splaiparser.partial_exon_deletion,
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)

    df['Exon_skipping'] = df.parallel_apply(
        splaiparser.exon_skipping, 
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)
                                            
    df['Int_Retention'] = df.parallel_apply(
        splaiparser.intron_retention, 
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)

    df['multiexs'] = df.parallel_apply(
        splaiparser.multi_exon_skipping, 
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)

    # df = pd.read_pickle('splai_vep_vcfs/hgmd_dm/all_DM_chr1.splai.vep.enst.introndist.exintinfo.splicing2.pkl')
    #7.   Annotate aberrant splicing size (bp)
    logger.info('Annotating aberrant splicing size (bp)...')
    #7-1. Annotate size of 
    df['Size_Part_ExDel'] = df.parallel_apply(
        splaiparser.anno_partial_exon_del_size, 
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)

    #7-3. Annotate size of partial intron retention
    df['Size_Part_IntRet'] = df.parallel_apply(
        splaiparser.anno_partial_intron_retention_size, 
        thresholds=thresholds_SpliceAI_parser,
        axis=1)

    #7-2. Annotate size of pseudoexon
    df['Size_pseudoexon'] = df.parallel_apply(
        splaiparser.anno_gained_exon_size, 
        thresholds=thresholds_SpliceAI_parser, 
        axis=1)

    #7-4. Annotate size of intron retention
    df['Size_IntRet'] = df.parallel_apply(
        splaiparser.anno_intron_retention_size, 
        thresholds=thresholds_SpliceAI_parser,
        axis=1)

    #7-5. Annotate size of exon skipping
    df['Size_skipped_exon'] = df.parallel_apply(
        splaiparser.anno_skipped_exon_size, 
        thresholds=thresholds_SpliceAI_parser,
        axis=1)

    df['variant_id'] = df['CHROM'].astype(str) + '-' \
        + df['POS'].astype(str) + '-' + df['REF'] + '-' + df['ALT']

    #8.   Evaluate splicing effects
    logger.info('Predicting CDS change...')
    #8-1. Predict CDS change
    df['CDS_Length'] = df.progress_apply(predeffect.calc_cds_len, db=db, axis=1)
    df['is_10%_truncation'] = df.progress_apply(predeffect.calc_cds_len_shorten, axis=1)

    #8-2. Determine if the gene is included in eLoFs genes
    df['is_eLoF'] = df.parallel_apply(predeffect.elofs_judge, axis=1)

    #8-3. Determine causing NMD or not
    df['is_NMD_at_Canon'] = df.parallel_apply(predeffect.nmd_judge, axis=1)

    #8-4. Frame check
    # Covert to str (Cannot predict splicing event) to np.nan


    cannot_predict: str = 'Cannot predict splicing event'
    df['Size_Part_ExDel'] = df['Size_Part_ExDel'].replace(cannot_predict, np.nan)
    df['Size_Part_IntRet'] = df['Size_Part_IntRet'].replace(cannot_predict, np.nan)
    df['Size_pseudoexon'] = df['Size_pseudoexon'].replace(cannot_predict, np.nan)
    df['Size_IntRet'] = df['Size_IntRet'].replace(cannot_predict, np.nan)
    df['Size_skipped_exon'] = df['Size_skipped_exon'].replace(cannot_predict, np.nan)

    df['is_Frameshift_Part_ExDel'] = df['Size_Part_ExDel'].parallel_apply(
        predeffect.frame_check)
    df['is_Frameshift_Part_IntRet'] = df['Size_Part_IntRet'].parallel_apply(
        predeffect.frame_check)
    df['is_Frameshift_pseudoexon'] = df['Size_pseudoexon'].parallel_apply(
        predeffect.frame_check)
    df['is_Frameshift_IntRet'] = df['Size_IntRet'].parallel_apply(
        predeffect.frame_check)
    df['is_Frameshift_skipped_exon'] = df['Size_skipped_exon'].parallel_apply(
        predeffect.frame_check)
    df['is_Frameshift'] = df[['is_Frameshift_Part_ExDel', 
                            'is_Frameshift_Part_IntRet', 
                            'is_Frameshift_pseudoexon', 
                            'is_Frameshift_IntRet', 
                            'is_Frameshift_skipped_exon'
                            ]].any(axis=1)

    #9.   CCRs
    logger.info('Annotating CCRs info...')
    #9-1. Annotate truncated regions 
    df['skipped_region'] = df.parallel_apply(
        splaiparser.anno_skipped_regions, axis=1)
    df['deleted_region'] = df.parallel_apply(
        splaiparser.anno_deleted_regions, 
        thresholds=thresholds_SpliceAI_parser, axis=1)

    #9-2. Intersect with CCRs
    logger.info('Annotate CCR score')
    df = predeffect.anno_ccr_score(df)

    df.to_pickle(f'splai_vep_vcfs/hgmd_dm/all_DM_chr{chr}.splai.vep.nondel.vcf.enst.prescore.pkl')
    logger.info(f"Comleted chr{chr}")
	

    df: pd.DataFrame = pd.concat([pd.read_pickle(f'splai_vep_vcfs/hgmd_dm/all_DM_chr{chr}.splai.vep.nondel.vcf.enst.prescore.pkl') for chr in chr_list])

#10.  Scoring
sccore_ths = {'clinvar_same_pos': 2,     
                'clinvar_same_motif': 1,
                'clinvar_else': 0,
                'non_canon_splai_lte_0.1_outside': -3,
                'non_canon_splai_lte_0.1_other': -2,
                'non_canon_splai_bet_0.1_0.2': 1,
                'non_canon_splai_gte_0.2': 2,
                'canon_strong': 6, 
                'canon_moderate': 5, 
                'frameshift_nmd_eloF': 7, 
                'frameshift_nmd_not_eloF': 3,
                'canon_splai_lte_0.1': -3,
                'canon_splai_bet_0.1_0.2': -1,
                'canon_splai_gte_0.2': 0}
scoring = Scoring(ths=sccore_ths)

logger.info('Annotating Screening scores...')
df.rename(columns={'MaxSpliceAI': 'maxsplai'}, inplace=True)
df['insilico_screening'] = df.parallel_apply(scoring.insilico_screening, axis=1)
df['clinvar_screening'] = df.parallel_apply(scoring.clinvar_screening, axis=1)
df['PriorityScore'] = df.parallel_apply(scoring.calc_priority_score, axis=1)



### FINALIZE ###
"""
Create an output table with the following columns:
#CHROM POS ID REF ALT QUAL FILTER INFO
ID, QUAL FILTER columns are not used in this table, then set them as ".".
INFO column should contain the following information:
Priority Score

"""