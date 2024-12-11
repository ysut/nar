from pysam import VariantFile

def remove_square_brackets(s: str) -> str:
    return s.replace("[", "").replace("]", "")

def anno_same_pos_vars(row, cln_bcf: VariantFile) -> str:
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
        return remove_square_brackets(str(samepos))


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
        return 'unk_SpliceType'
    
    return str(row['CHROM']), query_start, query_end

def anno_same_motif_vars(row, cln_bcf: VariantFile) -> str:
    region: tuple = _generate_query_pos(row)

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
            clinvar = next(recs)
        except StopIteration:
            break
        else:
            registered_var = f'{clinvar[7]}_{clinvar[0]}:{clinvar[1]}'
            samemotifs.append(registered_var)

    if samemotifs == []:
        return "No_ClinVar_info_found"
    else:
        return remove_square_brackets(str(samemotifs))