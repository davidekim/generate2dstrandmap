# *generate2dstrandmap*
This script generates a 2d beta strand map in SVG. The default script configuration colors glycines lime green, beta-bulges red, and 3-10 helices yellow.

![Example SVG](6d0tA.pdb.2dstrandmap.svg)

## Dependencies

pyrosetta

## Usage

```
./generate2dstrandmap.py 6d0tA.pdb
```

## Usage details

```
./generate2dstrandmap.py -h
```

## Usage notes

This script will often produce incorrect beta strand maps due to missing or incorrect beta bulge and/or beta strand (E) secondary structure assignments. If this occurs H-bonds may be displayed in a non-perpendicular orientation, the shear number may be an odd number, the strand count may be off, and/or other peculiarities may occur. You may fix these issues by manually assigning correct beta bulges and/or strand secondary structure which may require manual inspection of your protein structure and using the following command options:

```
--add_E <residue number(s), comma separated>
--rm_E <residue number(s), comma separated>
--add_bulges <residue number(s), comma separated>
--rm_bulges <residue number(s), comma separated>
```
