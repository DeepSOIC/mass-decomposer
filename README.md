# mass-decomposer
python script for processing mass-spectrum mass signals into molecule signals.

Takes as input:
* list of molecules and peak magnitudes for each molecule
* datalogged file with columns corresponding to recorded magnitude of each peak

Outputs:
* magnitudes of molecules to match the input peaks (least-squares fit, independently for each data row)
* residual data
* inverse matrix
* matrix condition value

## how to use
instructions are for Windows.

0. First of all, install python with numpy. You must be able to type `python` to cmd to launch python interpreter (i.e., python.exe must be in a directory listed in PATH environment variable).

1. Edit script.py to contain your list of molecules, and mass-spectrum peak magnitudes.

```
molecs = [
    ('NO', {
        'm29Int': 0.021892103,
        'm30Int': 0.072973677,
        'm31Int': 0.905134219,
    }),
    ('CO', {
        'm28Int': 0.436291522,
        'm29Int': 0.498264750,
        'm30Int': 0.058254834,
        'm31Int': 0.007188894,
    }),
    <and so on>
]
```
`m29Int` and the like must match column names in the input file. They are case-sensitive. 
Column names are parsed out from the first line of text in the input file that has 3 or more tab-delimited substrings. The remaining such strings are treated as data (numbers).

Magnitudes for peaks of one molecule should add up to around 1 (exactness not required, but if they add up to, say, 10, condition value may be useless) 

You can also specify weights for masses (also in script.py). You may want to decrease weights of "buggy" masses, such as 28, which includes CO outgassing of cathode, that may be affected by other gases being analyzed). Or you may want to increase some, to force calculation of some molecules from some peaks preferentially.

2. Apply some preprocessing to the input file. Particularly, subtract background

3. drag-drop the file to be processed onto `MassDecomposer.bat`.
Input file must be tab-delimited. Comma and point are both treated as decimal point; group separators are not supported.

-> two files should appear in the directory your original file is stored: one with "_proc" and one with "_proc_matrix" added to the name.

4. analysis:

  4.1 Inspect condition value. This value is approximately by how much the errors are being inflated. If >10, the decomposition is very ambiguous.

  4.2 Inspect residuals. Substantial magnitudes there mean you have some other signals in the mass-spectrum, and the magnitudes of molecules are inaccurate.

5. use the result for science!

...

#. PROFIT
