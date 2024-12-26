import pysam

def __remove_square_brackets(s: str) -> str:
    return s.replace("[", "").replace("]", "")

def __remove_quotations(s: str) -> str:
    return s.replace("'", "").replace('"', '')

def anno_same_pos_vars(row, cln_bcf: pysam.VariantFile) -> str:
    query_chr: str = f"{row['CHROM']}"
    query_pos: int = int(row['POS'])
    query_ref: str = row['REF']
    query_alt: str = row['ALT']
    query_variant: str = f"{query_chr}-{query_pos}-{query_ref}-{query_alt}"
    query_start: int = query_pos - 1
    query_end: int = query_pos
    recs = cln_bcf.fetch(query_chr, query_start, query_end)

    # main loop 
    samepos = []
    while 1:
        try:
            rec = next(recs)
        except StopIteration:
            break
        else:
            if rec.alts is None:
                rec_alt: str = "."
            else:
                rec_alt: str = rec.alts[0]
            rec_id = f"{rec.contig}-{rec.pos}-{rec.ref}-{rec_alt}"
            if query_variant == rec_id:
                clnsigs: list = [x for x in rec.info["CLNSIG"]]
                samepos.append(clnsigs)
            else:
                pass
    
    # Return the results
    if samepos == []:
        return "No_ClinVar_info_found"
    else:
        return __remove_square_brackets(str(samepos))


def _generate_query_pos(row) -> tuple:
    query_pos = int(row['POS'])

    if row['SpliceType'] == 'Donor_int':
        if row['Strand'] == '+':
            query_start: int = query_pos - int(row['IntronDist']) - 2 - 1
            query_end: int = query_pos - int(row['IntronDist']) + 6
        elif row['Strand'] == '-':
            query_start: int = query_pos + int(row['IntronDist']) - 6 - 1
            query_end: int = query_pos + int(row['IntronDist']) + 2
        else:
            return 'unk_Strand',

    elif row['SpliceType'] == 'Donor_ex':
        if row['Strand'] == '+':
            query_start: int = query_pos + int(row['exon_pos']) - 3 - 1
            query_end: int = query_pos + int(row['exon_pos']) + 5
        elif row['Strand'] == '-':
            query_start: int = query_pos - int(row['exon_pos']) - 5 - 1
            query_end: int = query_pos - int(row['exon_pos']) + 3
        else:
            return 'unk_Strand',
    
    elif row['SpliceType'] == 'Acceptor_int':
        if row['Strand'] == '+':
            query_start: int = query_pos + (- int(row['IntronDist'])) - 20 - 1
            query_end: int = query_pos + (- int(row['IntronDist'])) + 0
        elif row['Strand'] == '-':
            query_start: int = query_pos - (- int(row['IntronDist'])) - 0 - 1
            query_end: int = query_pos - (- int(row['IntronDist'])) + 20
        else:
            return 'unk_Strand',
    
    elif row['SpliceType'] == 'Acceptor_ex':
        if row['Strand'] == '+':
            query_start: int = query_pos - int(row['exon_pos']) - 19 - 1
            query_end: int = query_pos - int(row['exon_pos']) + 1
        elif row['Strand'] == '-':
            query_start: int = query_pos + int(row['exon_pos']) - 1 - 1
            query_end: int = query_pos + int(row['exon_pos']) + 19
        else:
            return 'unk_Strand',
        
    else:
        return 'unk_SpliceType',
    
    return str(row['CHROM']), query_start, query_end

def anno_same_motif_vars(row, cln_bcf: pysam.VariantFile) -> str:
    region: tuple = _generate_query_pos(row)
    # print(f"Query region: {region[0]}")

    if region[0] == 'unk_Strand':
        return 'unk_Strand'
    elif region[0] == 'unk_SpliceType':
        return 'unk_SpliceType'
    else:
        recs = cln_bcf.fetch(*region)

    # Main loop
    samemotifs = []
    while 1:
        try:
            rec = next(recs)
        except StopIteration:
            break
        else:
            if rec.alts is None:
                rec_alt: str = "."
            else:
                rec_alt: str = rec.alts[0]

            registered_var = f"{rec.contig}-{rec.pos}-{rec.ref}-{rec_alt}"
            clnsigs: list = [x for x in rec.info["CLNSIG"]]
            samemotifs.append(f"{registered_var}:{clnsigs}")

    if samemotifs == []:
        return "No_ClinVar_info_found"
    else:
        return __remove_quotations(__remove_square_brackets(str(samemotifs)))


def extract_same_motif_clinsigs(row) -> list:
    if row == "No_ClinVar_info_found":
        return ["No_ClinVar_info_found"]
    elif row.startswith('unk'):
        return ["unk_Strand_or_SpliceType"]
    else:
        # print(row)
        return [var_clinsig.split(':')[1] for var_clinsig in row.split(', ')]