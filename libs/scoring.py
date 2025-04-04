import pandas as pd

bp7_csq: set = {'intron_variant', 'synonymous_variant'}

# scores = {'clinvar_same_pos': 2,                    #1
#           'clinvar_same_motif': 1,                  #2
#           'clinvar_else': 0,                        #3 
#           'non_canon_splai_gte_0.2': 3,             #4
#           'non_canon_splai_bet_0.1_0.2': 0,         #5
#           'non_canon_splai_lte_0.1_outside': -2,    #6
#           'non_canon_splai_lte_0.1_other': -1,      #7
#           'frameshift_nmd_eloF': 7,                 #8
#           'frameshift_nmd_not_eloF': 3,             #9
#           'canon_strong': 6,                        #10
#           'canon_moderate': 5,                      #11
#           'canon_splai_lte_0.1': -2,                #12
#           'canon_splai_bet_0.1_0.2': -1,            #13
#           'canon_splai_gte_0.2': 1,                 #14
#           }

class Scoring:
    def __init__(self) -> None: 
        self.scores: dict = {}

    def recal_scores_in_canon(self, row) -> str:
        # if float(row['maxsplai']) < 0.1:
        if row['is_Canonical'] == "Yes":
            if float(row['maxsplai']) <= 0.1:
                return "s12"
            elif float(row['maxsplai']) < 0.2:
                return "s13"
            else:
                return "s14"
        else:
            return "s0"

    def insilico_screening(self, row) -> str:
        #0. No score
        try:
            maxsplai = float(row['maxsplai'])
        except ValueError:
            return "Not available"

        #1. Canonical
        if row['is_Canonical'] == "Yes":
            pre_score = self.recal_scores_in_canon(row)
            # Frameshift variants
            if row['is_Frameshift']:
                # print(f"Frameshift: {row['is_Frameshift']}")
                if ((row['is_NMD_at_Canon'] == 'Possibly_NMD') 
                    | (row['loftee'] == 'HC')
                    | (row['loftee'] == 'OS')):
                    if row['is_eLoF']:
                        raw_score = "s10"
                    else:
                        raw_score = "s11"
                else:
                    if ((float(row['skipped_ccrs']) >= 95) | (float(row['deleted_ccrs']) >= 95)):
                        raw_score = "s8"
                    else:
                        if row['is_10%_truncation']:
                            raw_score = "s8"
                        else:
                            raw_score = "s9"
            # In-frame
            else:
                if ((float(row['skipped_ccrs']) >= 95) | (float(row['deleted_ccrs']) >= 95)):
                    raw_score = "s8"
                else:
                    if row['is_10%_truncation']:
                        # print('≥10% Truncation')
                        raw_score = "s8"
                    else:
                        # print(f"≤10% Truncation {self.scores['canon_moderate']}")
                        raw_score = "s9"
        
        #2. Non-canonical
        else:
            if maxsplai >= 0.2:
                return "s7"
            elif maxsplai <= 0.1:
                if ((row['SpliceType'] == 'Acceptor_int') | (row['SpliceType'] == 'Donor_int')):
                    # if ((int(row['Int_loc']) <= -21) | (int(row['Int_loc']) >= 7)):
                    if ((int(row['IntronDist']) <= -21) | (int(row['IntronDist']) >= 7)):
                        raw_score = "s4"
                    else:
                        raw_score = "s5"
                elif ((row['SpliceType'] == 'Acceptor_ex') | (row['SpliceType'] == 'Donor_ex')):
                    csqs: list = row['Consequence'].split('&')
                    if not set(csqs).isdisjoint(bp7_csq):
                        if ((int(row['ex_up_dist']) > 1) & (int(row['ex_down_dist']) > 3)):
                            raw_score = "s4"
                        else:
                            raw_score = "s5"
                    else:
                        # Not in bp7_csq
                        raw_score = "s5"
                else:
                    raw_score = "s5"
            else:
                raw_score = "s6"
    
        # Calibrate minus scores to 0

        # if raw_score < 0:
        #     raw_score = 0

        return raw_score


    # def clinvar_screening(self, row) -> str:
    #     if str(row['insilico_screening']).item() == "Not available":
    #         return "Not available"
    #     else:
    #         if row['insilico_screening'] >= 0:
    #             if row['clinvar_same_pos']:
    #                 return self.scores['clinvar_same_pos']
    #             else:
    #                 if row['clinvar_same_motif']:
    #                     return self.scores['clinvar_same_motif']
    #                 else:
    #                     return self.scores['clinvar_else']
    #         else:
    #             return self.scores['clinvar_else']


    # def calc_priority_score(self, row):
    #     if row['insilico_screening'].item() == "Not available":
    #         return "Not available"
    #     else:
    #         return row['insilico_screening'] + row['clinvar_screening']

    # def clinvar_screening(self, row) -> str:
    #     if row['insilico_screening'] >= 0:
    #         if row['clinvar_same_pos']:
    #             return self.scores['clinvar_same_pos']
    #         else:
    #             if row['clinvar_same_motif']:
    #                 return self.scores['clinvar_same_motif']
    #             else:
    #                 return self.scores['clinvar_else']
    #     else:
    #         return self.scores['clinvar_else']
    # def clinvar_screening(self, row) -> str:
    #     if row['clinvar_same_pos']:
    #         return self.scores['clinvar_same_pos']
    #     else:
    #         if row['clinvar_same_motif']:
    #             return self.scores['clinvar_same_motif']
    #         else:
    #             return self.scores['clinvar_else']

    def clinvar_screening(self, row) -> str:
        cln_same_pos = row['clinvar_same_pos'].replace("'", "")
        if cln_same_pos in ['Benign', 'Likely_benign', 'Benign/Likely_benign']:
            return "s15"
        else:
            if cln_same_pos in ['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic']:
                return "s1"
            else:
                if 'Pathogenic' in row['same_motif_clinsigs']:
                    return "s2"
                elif 'pathogenic' in row['same_motif_clinsigs']:
                    return "S2"
                else:
                    return "s3"
    
    def calc_priority_score(self, row):
        # print(row['insilico_screening'] + row['clinvar_screening'])
        if row['insilico_screening'] == "Not available":
            return "Not available"
        else:
            total: int = int(row['insilico_screening'] + row['clinvar_screening'])
            if total < 0:
                return 0
            else:
                return total
        

    def calc_priority_score2(self, df: pd.DataFrame) -> pd.DataFrame:
        df['PriorityScore'] = df['insilico_screening'] + df['clinvar_screening']
        return df
    
    