import csv
from collections import defaultdict

infile  = {{ in.infile | str | repr }}
outfile = {{ out.outfile | str | repr }}
features = {{ args.features | repr }}
captured = {{ args.captured | ? int ! repr }}
samplecol = {{ args.samplecol | isinstance: str ?
                                lambda x: int(x) if x.isdigit() else x !
                              | repr }}
classcol = {{ args.classcol | isinstance: str ?
                                lambda x: int(x) if x.isdigit() else x !
                            | repr }}
genecol = {{ args.genecol | isinstance: str ?
                                lambda x: int(x) if x.isdigit() else x !
                            | repr }}

PPINGS = dict(
    # oncotator MAFs
    INTRON = "Intron",
    FIVE_PRIME_UTR = "5'UTR",
    THREE_PRIME_UTR = "3'UTR",
    IGR = "IGR",
    FIVE_PRIME_PRIME_FLANK = "5'Flank",
    THREE_PRIME_PRIME_FLANK = "3'Flank",
    MISSENSE = "Missense_Mutation",
    NONSENSE = "Nonsense_Mutation",
    NONSTOP = "Nonstop_Mutation",
    SILENT = "Silent",
    SPLICE_SITE = "Splice_Site",
    IN_FRAME_DEL = "In_Frame_Del",
    IN_FRAME_INS = "In_Frame_Ins",
    FRAME_SHIFT_INS = "Frame_Shift_Ins",
    FRAME_SHIFT_DEL = "Frame_Shift_Del",
    START_CODON_SNP = "Start_Codon_SNP",
    START_CODON_INS = "Start_Codon_Ins",
    START_CODON_DEL = "Start_Codon_Del",
    STOP_CODON_INS = "Stop_Codon_Ins",
    STOP_CODON_DEL = "Stop_Codon_Del",
    # Note: A STOP_CODON_SNP is a nonstop mutation (or Silent),
    DE_NOVO_START_IN_FRAME = "De_novo_Start_InFrame",
    DE_NOVO_START_OUT_FRAME = "De_novo_Start_OutOfFrame",
    RNA = "RNA",
    LINCRNA = "lincRNA",
    # vcf2maf,
    TRANSLATION_START_SITE = 'Translation_Start_Site',
    SPLICE_REGION = "Splice_Region",
    TARGETED_REGION = "Targeted_Region",
)

# See: https://github.com/PoisonAlien/maftools/blob/6af2fe77d87bdf078e81d6d143bbf626825804d8/R/read_maf_dt.R#L97
MAF_NONSYN = set([
    'FRAME_SHIFT_DEL',
    'FRAME_SHIFT_INS',
    'SPLICE_SITE',
    'TRANSLATION_START_SITE',
    'NONSENSE',
    'NONSTOP',
    'IN_FRAME_DEL',
    'IN_FRAME_INS',
    'MISSENSE'
])

IS_MAF = infile.endswith('.maf')

def get_colnames(samcol, clscol, gcol):
    if isinstance(samcol, str) and isinstance(clscol, str) and isinstance(gcol, str):
        return samcol, clscol, gcol

    if IS_MAF and samcol is None:
        samcol = 'Tumor_Sample_Barcode'
    elif samcol is None:
        samcol = 0

    if IS_MAF and clscol is None:
        clscol = 'Variant_Classification'

    if IS_MAF and gcol is None:
        gcol = 'Hugo_Symbol'

    with open(infile) as fp:
        reader = csv.DictReader(filter(lambda row: row[0]!='#', fp),
                                delimiter='\t')
        if isinstance(samcol, int):
            samcol = reader.fieldnames[samcol]
        if isinstance(clscol, int):
            clscol = reader.fieldnames[clscol]
        if isinstance(gcol, int):
            clscol = reader.fieldnames[gcol]

    return samcol, clscol, gcol

samplecol, classcol, genecol = get_colnames(samplecol, classcol, genecol)

def get_samples(samcol):
    samples = []

    with open(infile) as fp:
        reader = csv.DictReader(filter(lambda row: row[0]!='#', fp),
                                delimiter='\t')
        for row in reader:
            if row[samcol] not in samples:
                samples.append(row[samcol])
    return samples

def tmb(feat):
    samcol, clscol = samplecol, classcol
    ret = defaultdict(lambda: 0)

    with open(infile) as fp:
        reader = csv.DictReader(filter(lambda row: row[0]!='#', fp),
                                delimiter='\t')

        CLASS_COL_FLAG = None
        for row in reader:
            sample = row[samcol]
            if CLASS_COL_FLAG is None:
                if clscol is None:
                    ret[sample] += 1
                    CLASS_COL_FLAG = 'ALL'
                else:
                    klass = row[clscol]
                    if klass in ('0', '1', 0, 1):
                        klass = int(klass)
                        CLASS_COL_FLAG = 'BIN'
                        if klass:
                            ret[sample] += 1
                    else:
                        CLASS_COL_FLAG = 'NAME'
                        if klass.lower == 'nonsyn':
                            ret[sample] += 1
            elif CLASS_COL_FLAG == 'NAME':
                if row[clscol].lower() == 'nonsyn':
                    ret[sample] += 1
            elif CLASS_COL_FLAG == 'BIN':
                if int(row[clscol]):
                    ret[sample] += 1
            else:
                ret[sample] += 1
    if captured:
        ret = {key: float(val)/float(captured) for key, val in ret.items()}
    return ret

def gmut(feat):
    gene = feat[5:]
    samcol, gcol = samplecol, genecol
    ret = defaultdict(lambda: lambda: 0)

    with open(infile) as fp:
        reader = csv.DictReader(filter(lambda row: row[0]!='#', fp),
                                delimiter='\t')

        for row in reader:
            sample = row[samcol]
            if row[gcol] == gene:
                ret[sample] += 1
    return ret

def bgmut(feat):
    ret = gmut(feat[1:])
    return {sample: str(bool(count)).upper()
            for sample, count in ret.items()}

FEAT_FUNCS = {
    'tmb': tmb,
    'gmut': gmut,
    'bgmut': bgmut
}

if __name__ == "__main__":
    samples = get_samples(samplecol)

    with open(outfile, 'w') as fout:
        fout.write('\t'.join([""] + [
            f'{feat[6:]}_mut'
            if feat.startswith('bgmut-')
            else f'{feat[5:]}_mut'
            if feat.startswith('gmut-')
            else feat.upper().replace('-', '_')
            for feat in features
        ]) + '\n')
        for sample in samples:
            featvals = []
            for feat in features:
                featvals.append(
                    str(FEAT_FUNCS[feat](feat).get(sample, 'NA'))
                    if 'gmut' not in feat
                    else FEAT_FUNCS['gmut'](feat)
                    if feat.startswith('gmut-')
                    else FEAT_FUNCS['bgmut'](feat)
                )
            fout.write('\t'.join([sample] + featvals) + '\n')
