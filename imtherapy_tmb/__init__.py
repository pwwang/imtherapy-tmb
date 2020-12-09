from pathlib import Path
from pipen import Proc
from imtherapy.modules import FTModule, ft_modules
from imtherapy.envs import envs

HERE = Path(__file__).parent.resolve()

class FTTmb(Proc):

    input_keys = 'infile:file'
    output = 'outfile:file:tmb.txt'
    script = 'file://{HERE}/ft_tmb.py'
    envs = envs
    args = {}

class FeatureTransformTmb(FTModule):
    """Transform a MAF file into Tumor mutation burden for each sample"""
    name = 'tmb'
    process = FTTmb
    start_process = end_process = process

    @ft_modules.impl
    def on_args_init(self, params):
        params.add_param(
            'tmb',
            type='ns',
            show=False,
            desc='Options for tmb feature transform module'
        )
        params.add_param(
            'tmb.maf',
            type='path',
            show=False,
            argname_shorten=False,
            desc='The MAF file to calculate the TMB.',
            callback=lambda val, all_vals: (
                ValueError('A .maf file is required to calculate the TMB.')
                if not val and self.name in all_vals.t
                else val
            )
        )
        params.add_param(
            'tmb.captured',
            type=str,
            show=False,
            argname_shorten=False,
            desc=('Either the length the captured regions. '
                  'Required for stratified TMB. `K/M` is supported '
                  'for kilo or mega bases.'),
            callback=lambda val, all_vals: (
                ValueError('Required for stratified TMB.')
                if not val and all_vals.tmb.method == 'stratified'
                else int(val[:-1]) * 1_000
                if val[-1].upper() == 'K'
                else int(val[:-1]) * 1_000_000
                if val[-1].upper() == 'M'
                else int(val)
            )
        )
        params.add_param(
            'tmb.method',
            show=False,
            argname_shorten=False,
            type='choice',
            default='count',
            choices=['count', 'stratified'],
            desc=('The method to calculate the TMB. '
                  'Stratified TMB requires the capture region and will be '
                  'stratified per megabase.')
        )

    @ft_modules.impl
    def on_args_parsed(self, args):
        self.process.input = [args.tmb.maf]
        self.process.args['method'] = args.tmb.method
