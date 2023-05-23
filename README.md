# mtSVG
Convert a mtDNA GFF returned by MITOS2 to a linear SVG representation

# Installation

Require `Python3`, the script `mtSVG.py` and the library `drawsvg` which can be installed as follows. 
```
python3 -m pip install "drawsvg~=2.0"
```

# Usage

Example:
```
./mtSVG.py --gff example.gff --length 14934 --cox2cob true
```

# Output

A linear representation of the mtDNA genes present in the GFF exported as a SVG image.

![alt text](https://rehost.diberie.com/Picture/Get/f/177361)

