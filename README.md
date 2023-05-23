# mtSVG

Convert a mtDNA GFF returned by MITOS2 to a linear SVG representation

## Installation

Require `Python3`, the script `mtSVG.py` and the library `drawsvg` which can be installed as follows. 

```
python3 -m pip install "drawsvg~=2.0"
```

## Usage

Example:
```
usage: mtSVG.py [-h] [--gff GFF] [--svg SVG] [--length LENGTH]
                [--cox2cob COX2COB]
                
arguments:
  -h, --help         show this help message and exit
  --gff <string>          The path to the gff file [REQUIRED]
  --svg <string>          The path of the SVG to create [OPTIONAL, default=linear_mtdna.svg]
  --length <int>          The mtDNA length in bp to scale the output [REQUIRED]
  --start <string>        Gene to use at the start of the ribbon [OPTIONAL, default=cox1]
  --reverse <string>      If true, the gene order in the gff is reversed [OPTIONAL, default=false]
```

## Example

```
./mtSVG.py --gff example.gff --length 14934 --cox2cob true
```

## Output

A linear representation of the mtDNA genes present in the GFF exported as a SVG image.

![alt text](https://rehost.diberie.com/Picture/Get/f/177361)

