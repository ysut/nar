import pysam

def anno_same_pos_vars(row, tabixfile: pysam.pysam.libctabix.TabixFile):
    samepos = []
    query_chr: str = f"{row['CHROM']}"
    query_pos: int = row['POS']
    query_ref: int = row['REF']
    query_alt: int = row['ALT']
    query_variant: str = f"{query_chr}-{query_pos}-{query_ref}-{query_alt}"
    query_start: int = int(query_pos) - 1
    query_end: int = int(query_pos) + 1
    clinvars = tabixfile.fetch(
        query_chr, query_start, query_end, parser=pysam.asBed())

    # main loop 
    while 1:
        try:
            clinvar = next(clinvars)
        except StopIteration:
            break
        else:
            clinvar_variant = f"{clinvar.contig}-{clinvar.start}-{clinvar[3]}-{clinvar[4]}"
            if query_variant == clinvar_variant:
                samepos.append(f'{clinvar[7]}_{clinvar[0]}:{clinvar[1]}')
            else:
                pass

    return samepos

def anno_same_motif_vars(row, tabixfile: pysam.pysam.libctabix.TabixFile):
    """ Fetch the variants within the same splicing motif from ClinVar data

    Args:
        row (_type_): _description_

    Returns:
        list: Return a list of variants registered in ClinVar 
                within the same splicing motif (Prediction_CHROM:POS)
    """
    
    # Initialize result list & set common variables
    samemotifs = []
    query_chr: str = f'{row["CHROM"]}'
    query_pos: int = int(row['POS'])

    # Generate query positions 
    if row['SpliceType'] == 'Donor_int':
        if row['Strand'] == '+':
            query_start: int = query_pos - int(row['IntronDist']) - 2 - 1
            query_end: int = query_pos - int(row['IntronDist']) + 6
        elif row['Strand'] == '-':
            query_start: int = query_pos + int(row['IntronDist']) - 6 - 1
            query_end: int = query_pos + int(row['IntronDist']) + 2
        else:
            return 'unk_Strand'

    elif row['SpliceType'] == 'Donor_ex':
        if row['Strand'] == '+':
            query_start: int = query_pos + int(row['exon_pos']) - 3 - 1
            query_end: int = query_pos + int(row['exon_pos']) + 5
        elif row['Strand'] == '-':
            query_start: int = query_pos - int(row['exon_pos']) - 5 - 1
            query_end: int = query_pos - int(row['exon_pos']) + 3
        else:
            return 'unk_Strand'
    
    elif row['SpliceType'] == 'Acceptor_int':
        if row['Strand'] == '+':
            query_start: int = query_pos + (- int(row['IntronDist'])) - 20 - 1
            query_end: int = query_pos + (- int(row['IntronDist'])) + 0
        elif row['Strand'] == '-':
            query_start: int = query_pos - (- int(row['IntronDist'])) - 0 - 1
            query_end: int = query_pos - (- int(row['IntronDist'])) + 20
        else:
            return 'unk_Strand'
    
    elif row['SpliceType'] == 'Acceptor_ex':
        if row['Strand'] == '+':
            query_start: int = query_pos - int(row['exon_pos']) - 19 - 1
            query_end: int = query_pos - int(row['exon_pos']) + 1
        elif row['Strand'] == '-':
            query_start: int = query_pos + int(row['exon_pos']) - 1 - 1
            query_end: int = query_pos + int(row['exon_pos']) + 19
        else:
            return 'unk_Strand'
        
    else:
        return 'unk_SpliceType'

    clinvars = tabixfile.fetch(
        query_chr, query_start, query_end, parser=pysam.asBed())
        
    # Main loop
    while 1:
        try:
            clinvar = next(clinvars)
        except StopIteration:
            break
        else:
            registered_var = f'{clinvar[7]}_{clinvar[0]}:{clinvar[1]}'
            samemotifs.append(registered_var)

    return samemotifs