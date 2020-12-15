import sys
from pathlib import Path
from pipen import Proc
from imtherapy.modules import FTModule, ft_modules
from imtherapy.envs import envs

HERE = Path(__file__).parent.resolve()

class FTMut(Proc):
    """Mutation transformers for imtherapy"""
    lang = sys.executable
    input_keys = 'infile:file'
    output = 'outfile:file:mut.txt'
    script = f'file://{HERE}/ft_mut.py'
    envs = envs
    args = {'features': [], 'captured': None}

class FeatureTransformMut(FTModule):
    """Transform mutation related features for imtherapy"""
    name = 'mut'
    process = FTMut
    start_process = end_process = process

    @ft_modules.impl
    def on_args_init(self, params):
        params.add_param(
            'mut',
            type='ns',
            show=False,
            desc='Options for mut feature transform module'
        )
        params.add_param(
            'mut.mutfile',
            type='path',
            show=False,
            argname_shorten=False,
            desc='A tab-delimited mutation file.',
            callback=lambda val, all_vals: (
                ValueError('A mutfile is required '
                           'for related feature transformation.')
                if not val and self.name in all_vals.t
                else val
            )
        )
        params.add_param(
            'mut.captured',
            type=str,
            show=False,
            argname_shorten=False,
            desc=('The length of the captured regions. '
                  'If provided, tumor mutation burden will be stratified. '
                  '`K/M` is supported for kilo or mega bases.'),
            callback=lambda val: (
                None
                if not val
                else int(val[:-1]) * 1_000
                if val[-1].upper() == 'K'
                else int(val[:-1]) * 1_000_000
                if val[-1].upper() == 'M'
                else int(val)
            )
        )
        params.add_param(
            'mut.samplecol',
            type=str,
            show=False,
            argname_shorten=False,
            desc=('The column of the samples. If mutfile is a MAF file, '
                  '`Tumor_Sample_Barcode` will be used. Otherwise, first '
                  'column will be used. Could be 0-based index.')
        )
        params.add_param(
            'mut.classcol',
            type=str,
            show=False,
            argname_shorten=False,
            desc=('The column for the class of the mutation, which is, in '
                  'most cases, used to check if a mutation is nonsynonymous '
                  'or not. If input file is MAF, will be inferred from '
                  '`Variant_Classification`. Could be 0-based index.',
                  '- If this column is binarized (0 or 1, telling from the '
                  'first record), then 1 is treated as nonsyn',
                  '- Otherwise, it should have `nonsyn` to indicate nonsyn '
                  'mutations.'
                  '- If this is not specified, all mutations will be treated '
                  'as nonsyn mutations.')
        )
        params.add_param(
            'mut.feats',
            show=False,
            argname_shorten=False,
            type=list,
            default=['tmb'],
            desc=('The features to be calculated from the mutfile.')
        )

    @ft_modules.impl
    def on_args_parsed(self, args):
        self.process.input = [args.mut.mutfile]
        self.process.args['features'] = args.mut.feats
        self.process.args['captured'] = args.mut.captured
        self.process.args['samplecol'] = args.mut.samplecol
        self.process.args['classcol'] = args.mut.classcol
