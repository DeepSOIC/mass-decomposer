#example!
molecs = [
    ('NO', {
        'm29Int': 0.021892103,
        'm30Int': 0.072973677,
        'm31Int': 0.905134219,
    }),
    ('CO', {
        #'m28Int': 0.436291522,
        'm29Int': 0.498264750,
        'm30Int': 0.058254834,
        'm31Int': 0.007188894,
    }),
    ('N2O', {
        #'m28Int': 0.000321429,
        'm29Int': 0.012214286,
        'm30Int': 0.123892857,
        'm31Int': 0.149285714,
        'm44Int': 0.001785714,
        'm45Int': 0.067857143,
        'm46Int': 0.644642857,
    }),
    ('N2', {
        #'m28Int': 0.0025,
        'm29Int': 0.095,
        'm30Int': 0.9025,
    }),    
    ('CO2', {
        'm44Int': 1.0,
    }),    
]

import numpy as np

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
    rank = np.linalg.matrix_rank(matrix)
    print("   computing {n_molec} values from {n_used_masses} values with a matrix of rank {rank}"
        .format(n_molec= len(molecs), n_used_masses= sum(mass_is_used), rank= rank))
    if rank < len(molecs):
        raise RuntimeError("Rank of matrix is less than the number of variables. Unsolvable.")
    return matrix, mass_index, molec_index, mass_is_used

def analyzeMatrix(molecs, matrix, outpath = None):
    nmasses, nmolecs = matrix.shape
    b = np.identity(nmasses)
    solution_matrix, residual, rank, singular_vals = np.linalg.lstsq(matrix, b, rcond= None)
    condition_value = singular_vals[0]/singular_vals[-1]
    print('    matrix condition value (more is worse): {cond}'.format(cond= condition_value))
    if outpath:
        with open(appendedToFileName(outpath, '_matrix'), 'w') as f:
            f.write('molecs = ')
            import pprint
            pprint.pprint(molecs, stream= f)
            print('condition value = {cond}'.format(cond= condition_value), file= f)
            print('matrix', file= f)
            np.savetxt(f, matrix, delimiter='\t')
            print('inverse matrix', file= f)
            np.savetxt(f, solution_matrix, delimiter='\t')
    
def appendedToFileName(path, suffix):
    i_dot = path.rfind('.')
    if i_dot == -1:
        i_dot = len(path)
    return path[0:i_dot] + suffix + path[i_dot:]

def solvePoint(matrix, vals):
    x, residual, rank, singular_vals = np.linalg.lstsq(matrix, vals, rcond= None)
    residuals = np.array(vals) - matrix.dot(x)
    return x, residuals

def val(txt):
    txt = txt.replace(",", ".")
    try:
        return float(txt)
    except ValueError:
        return None

def processFile(filepath, molecs, outpath = None):
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
                        analyzeMatrix(molecs, matrix, outpath)
                        molec_list = [molec[0] for molec in molecs]
                        used_mass_list = ['rd' + mass_id for (i, mass_id) in enumerate(masses) if mass_is_used[i]]
                        output.append('\t'.join(masses + molec_list + used_mass_list))
                else: #in data
                    if len(split) >= len(masses):
                        vals = [val(it) for it in split]
                        molv, residuals = solvePoint(matrix, vals)
                        outputvals = vals + list(molv) + [v for (i, v) in enumerate(residuals) if mass_is_used[i]]
                        output.append(
                            '\t'.join([str(v) for v in outputvals])
                        )
                    else:
                        print('    line #{num} has too few values, skipped.'.format(num= line_number))
    print('writing to {f}'.format(f= outpath))
    with open(outpath, 'w') as outfile:
        for line in output:
            outfile.write(line)
            outfile.write('\n')