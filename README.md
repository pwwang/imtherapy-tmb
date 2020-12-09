# imtherapy-tmb

Tumor mutation burden feature transform module for `imtherapy`

## Installation

```shell
pip install -U imtherapy-tmb
```

## Usage

```shell
imtherapy --tmb.maf <maf file> ... [other imtherapy options]
# or get stratified TMB / megabase
imtherapy --tmb.maf <maf file> \
    --tmb.method stratified \
    --tmb.captured 30M \
    ... [other imtherapy options]
```
