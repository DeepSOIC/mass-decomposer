__doc__ = """The engine of mass-decomposer
The main function to be used is processFile()"""

import numpy as np

def processFile(filepath, molecs, weights = {}, outpath = None, decimal_separator = '.'):
    """processFile(filepath, molecs, weights = {}, outpath = None, decimal_separator = '.'):
    molecs: list of tuples, each tuple containing a name and a dict of peak values and 
        amplitudes [('molecule_name', {'col_name':0.35,...}),...]
    weights: optional. Dict of weights for masses {'col_name':0.1}. Weights of columns not listed in
        this list are assumed to be 1.0.
    outpath: optional, file name for output file.
    decimal_separator: optional, sets which separator is used for writing output file. Does 
        not affect reading - there, both ',' and '.' are recognized as a decimal separator."""
    if outpath is None:
        outpath = appendedToFileName(filepath, '_proc')
            
    output = []
    
    print('reading {f}'.format(f= outpath))
    with open(filepath) as infile:
        in_header = True
        for line_number,line in enumerate(infile):
            split = line[:-1].split('\t')
            if len(split)>3:
                if in_header:
                    if val(split[0]) is None:
                        print('    header found: line #{num}'.format(num= line_number))
                        in_header = False
                        masses = split
                        matrix, mass_index, molec_index, mass_is_used = makeMatrix(molecs, masses)
                        weights_vec = makeWeightsVec(masses, mass_index, weights)
                        analyzeMatrix(molecs, matrix, weights_vec, outpath)
                        molec_list = [molec[0] for molec in molecs]
                        used_mass_list = ['rd' + mass_id for (i, mass_id) in enumerate(masses) if mass_is_used[i]]
                        output.append('\t'.join(masses + molec_list + used_mass_list))
                else: #in data
                    if len(split) >= len(masses):
                        vals = [val(it) for it in split]
                        molv, residuals = solvePoint(matrix, vals, weights_vec)
                        outputvals = vals + list(molv) + [v for (i, v) in enumerate(residuals) if mass_is_used[i]]
                        output.append(
                            '\t'.join([toString(v, decimal_separator) for v in outputvals])
                        )
                    else:
                        print('    line #{num} has too few values, skipped.'.format(num= line_number))
    print('writing to {f}'.format(f= outpath))
    with open(outpath, 'w') as outfile:
        for line in output:
            outfile.write(line)
            outfile.write('\n')


# ------------helper routines--------------

if not 'Err' in vars():
    def Err(*args, **kwargs):
        import sys
        print(*args, file=sys.stderr, **kwargs)

def makeMatrix(molecs, masses):
    """makeMatrix(molecs, masses): 
    molecs: list of tuples : [(molec_name, peak_dict),..]
    masses: list of strings (mass ids) : ['me32int', ...]
    returns: (matrix, mass_index, molec_index)"""
    #matrix converts molecs into masses
    # m29      ## ## ##                      
    # m30      ## ## ##       NO               
    # m31  = ( ## ## ## ) * ( CO  )                                 
    # m32      ## ## ##       N2O      
    # ...      ## ## ##             
    matrix = np.zeros( (len(masses), len(molecs)) )
    mass_index = {mass_id:i for (i, mass_id) in enumerate(masses)}
    molec_index = {molec[0]:i for (i, molec) in enumerate(molecs)}
    mass_is_used = [False]*len(masses)
    for molec_name, molec_masses in molecs:
        for mass_id, frac in molec_masses.items():
            matrix[mass_index[mass_id], molec_index[molec_name]] = frac
            mass_is_used[mass_index[mass_id]] = True
    
    #test rank
    rank = np.linalg.matrix_rank(matrix)
    print("    computing {n_molec} values from {n_used_masses} values with a matrix of rank {rank}"
        .format(n_molec= len(molecs), n_used_masses= sum(mass_is_used), rank= rank))
    if rank < len(molecs):
        raise RuntimeError("Rank of matrix is less than the number of variables. Solution is ambiguous.")
    
    return matrix, mass_index, molec_index, mass_is_used

def analyzeMatrix(molecs, matrix, weights, outpath = None):
    """analyzeMatrix(molecs, matrix, weights, outpath = None): analyzes matrix, prints basic
    info to stdout, and extended info to a file."""
    
    #generate solution matrix. The solution matrix is actually here just for reference, 
    #it is not used for actual calculations, as those should support dynamically changing weights 
    nmasses, nmolecs = matrix.shape
    vals = np.identity(nmasses)
    solvematrix = matrix #solvematrix is a modified matrix with weights burned into it
    solvevals = vals
    if weights is not None:
        solvematrix = (matrix.T * weights).T
        solvevals = vals.dot(weights)
    solution_matrix, residual, rank, singular_vals = np.linalg.lstsq(solvematrix, solvevals, rcond= None)
    
    #print condition value
    condition_value = singular_vals[0]/singular_vals[-1]
    print('    matrix condition value (more is worse): {cond}'.format(cond= condition_value))
    
    #write extended info to file
    if outpath:
        with open(appendedToFileName(outpath, '_matrix'), 'w') as f:
            import pprint
            f.write('molecs = ')
            pprint.pprint(molecs, stream= f)
            print('weights', file= f)
            np.savetxt(f, weights, delimiter='\t')
            print('condition value = {cond}'.format(cond= condition_value), file= f)
            print('matrix (unweighted)', file= f)
            np.savetxt(f, matrix, delimiter='\t')
            print('inverse matrix', file= f)
            np.savetxt(f, solution_matrix, delimiter='\t')
    
def appendedToFileName(path, suffix):
    i_dot = path.rfind('.')
    if i_dot == -1:
        i_dot = len(path)
    return path[0:i_dot] + suffix + path[i_dot:]

def solvePoint(matrix, vals, weights = None):
    """converts a single set of mass values to molecules. Returns result and residuals as arrays of values."""
    solvematrix = matrix
    solvevals = vals
    if weights is not None:
        solvematrix = (matrix.T * weights).T
        solvevals = np.array(vals) * weights
    x, residual, rank, singular_vals = np.linalg.lstsq(solvematrix, solvevals, rcond= None)
    residuals = np.array(vals) - matrix.dot(x)
    return x, residuals

def val(txt):
    """val(txt): converts a string into a float number."""
    txt = txt.replace(",", ".")
    try:
        return float(txt)
    except ValueError:
        return None

def toString(number, dec_sep):
    """toString(number, dec_sep): converts a number to a string"""
    st = str(number)
    if dec_sep != '.':
        return st.replace('.', dec_sep, 1)
    else:
        return st

def makeWeightsVec(masses, mass_index, weights_dict):
    """converts weights in a dict to an array of numbers"""
    weights_vec = np.ones(len(masses))
    for mass_id, wgt in weights_dict.items():
        weights_vec[mass_index[mass_id]] = wgt
    return weights_vec

