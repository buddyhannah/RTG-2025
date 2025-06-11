import sympy as sp
from sympy import ImmutableMatrix
from functools import reduce
import os


'''
    Gets the full path for an output file in the quantum_output directory 
    within the home directory. Creates quantum_output it doesn't exist.

    Input:
    filename - Name of the output file

    Output:
    Full path to the output file
'''
def get_data_file_path(filename):
    home_dir = os.path.expanduser("~")
    output_dir = os.path.join(home_dir, "quantum_output")
    os.makedirs(output_dir, exist_ok=True)
    return os.path.join(output_dir, filename)


'''
    Writes output to the console and file
    
    Input:
        text: Text to output
        prt: Boolean indicating whether to print to console
        output_file: File to write to
'''
def write_output(text, prt, output_file):
    if prt:
        print(text)
        if output_file:
            output_file.write(text + "\n")



"""
    Computes R^n
    
    Input:
    R: Dictionary of matrices and their labels
    n: Positive integer indicating the power to compute
    verbose=False: Whether to show the computational steps
    output_path=None: File to write output to, or None if not outputing to file
        
    Output:
    Dictionary of basis matrices and their labels
   
    """
def compute_R_n(R, n, verbose=False, output_path=None):
    output_file = None
    if output_path:
        output_file = open(output_path, 'a')


    
    # Print parameters
    write_output("*"*30, True, output_file)
    write_output("INPUT PARAMETERS:", True, output_file)
    write_output("*"*30, True, output_file)
    
    write_output("\nR=\n{", True, output_file)
    for i, (mat, name) in enumerate(R.items(), 1):
        write_output(f"{name}=", True, output_file)
        write_output(f"{sp.pretty(mat)}\n", True, output_file)
    write_output("}", True, output_file)
    write_output(f"Find R^{n}", True, output_file)


    write_output("\n"+"*"*30, verbose,  output_file)
    write_output("COMPUTING:", verbose,  output_file)
    write_output("*"*30, verbose, output_file)
    
    try:
        # Track all matrices
        basis = R.copy()  

        # Track matrices in the current step
        current_span = set(R.keys())  
        
        for step in range(1, n):
            if verbose:
                write_output(f"\nStep {step}: Computing R^{step+1}", verbose, output_file)
            
            new_span = set()

            for mat1 in  current_span:
                for mat2 in R:
                    # Compute product
                    product = mat1 * mat2
                    name = f"{basis[mat1]}{basis[mat2]}"
                    
                    if verbose:
                        write_output(f"\nProduct ({basis[mat1]})({basis[mat2]}):", verbose,  output_file)
                        write_output(sp.pretty(product), verbose,  output_file)
                    
                    # Simplify matrix
                    linearCombo = simplify_matrix(product, basis)
                    
                    # Matrix is linear combination
                    if linearCombo:
                        if verbose:
                            write_output(f"= {linearCombo}", verbose, output_file)
                    # New matrix found
                    else:
                        basis[product] = name
                        new_span.add(product)
                        if verbose:
                            write_output("New matrix", verbose, output_file)

        
            # Stop is last spanning set equals current spanning set
            if not new_span:
                write_output(f"\nSTABILIZED: R^{step} = R^{step+1} = ...", True, output_file)
                break

            current_span = new_span
        
        # Final spanning set
        write_output("\n"+"*"*30, True,  output_file)
        write_output("SOLUTION:", True,  output_file)
        write_output("*"*30, True,  output_file)

        write_output(f"\nFinal spanning set for R^{n}:", True, output_file)
        for i, (mat, name) in enumerate(basis.items(), 1):
            write_output(f'{i}. {name} ', True,  output_file)
        
        
        return basis
    
    finally:
        # Done writing to file
        write_output(f"Wrote output to {output_path}", output_path!=None,  None)
        write_output("\n" + "═"*50 + "\n", True, output_file)
        if output_file:
            output_file.close()
            

'''
    Tries to output_path a matrix as a linear combination of basis matrices.

    Input:
    mat: Matrix to decompose 
    basis: Dictionary of {matrix: label} pairs representing the basis

    Output:
    Returns the linear combination as a string of form "c1*matrix1 + c2*matrix2 + ..."
    if a solution exists, and None otherwise
'''
def simplify_matrix(mat, basis):

    # Check if the matrix is the zero matrix
    zero_mat = sp.zeros(*mat.shape)
    if mat == zero_mat:
        return "zero matrix"


    # Convert matrix-label dictionary to list of matrices
    basis_list = list(basis)

    # Create symbolic coefficients 
    coeffs = sp.symbols(f'c0:{len(basis_list)}')
    
    # Add each (coefficient, matrix) pair to lin_combo
    lin_comb = sp.zeros(*mat.shape)
    for c, b in zip(coeffs, basis_list):
        lin_comb += c * b

    # Set up the system of equations
    equations = []
    # For every entry in mat, create the equation mat[i,j] == lin_comb[i,j]
    # Goal: Find c0, c1, etc., so that mat = lin_comb
    for i in range(mat.rows):
        for j in range(mat.cols):
            equations.append(sp.Eq(mat[i,j], lin_comb[i,j]))
    

    # Solve the system of linear equations
    # Return solutions as a list of dictionaries, like: [{'c0': 2, 'c1': 3}]
    solution = sp.solve(equations, coeffs, dict=True)

    # Solution found
    if solution:
        # Pick first solution
        sol = solution[0]

        # Only keep non-zero terms. Example, {c0: 2, c1: 0} --> {c0: 2}.
        non_zero = {c: val for c, val in sol.items() if val != 0}
    
        # Return linear combination as string
        terms = []
        for c, val in non_zero.items():
            b = basis_list[coeffs.index(c)]
            b_str = basis[b]
            terms.append(f"{val}*{b_str}")
        return " + ".join(terms)
    
    # No solution, return none
    return None


def main():


    I = ImmutableMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) 
    A = ImmutableMatrix([[0, -1*sp.I, 0], [sp.I, 0, 0], [0, 0, 0]]) 
    B = ImmutableMatrix([[0, 0, -1*sp.I], [0, 0, 0], [sp.I, 0, 0]]) 
    C = ImmutableMatrix([[0, 0, 0], [0, 0, -1*sp.I], [0, sp.I, 0]]) 
    R = {I:'I', B:'B', C:'C'}
    n = 5
    compute_R_n(R, n, verbose=True)


    I = ImmutableMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) 
    A = ImmutableMatrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) 
    B = ImmutableMatrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) 
    C = ImmutableMatrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) 
    R = {I:'I', A: 'A', B:'B', C:'C'}
    n = 2
    compute_R_n(R, n, verbose=True)



    I = ImmutableMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) 
    A = ImmutableMatrix([[0, -1*sp.I, 0], [sp.I, 0, 0], [0, 0, 0]]) 
    B = ImmutableMatrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) 
    C = ImmutableMatrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) 
    R = {I:'I', A: 'A', B:'B', C:'C'}
    n = 1
    compute_R_n(R, n, verbose=True)


    A = ImmutableMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) 
    B = ImmutableMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 2]]) 
    R = {A: 'A', B:'B'}
    n = 2
    compute_R_n(R, n, verbose=True)
        
    A = ImmutableMatrix([[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])  
    R = {A: 'A'}
    n = 5
    compute_R_n(R, n, verbose=True)

    I = ImmutableMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  
    A = ImmutableMatrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) 
    B = ImmutableMatrix([[0, -1*sp.I, 0], [sp.I, 0, 0], [0, 0, 0]])   
    R = {I: 'I', A: 'A', B: 'B'}
    n = 3
    compute_R_n(R, n, verbose=True)

    A = ImmutableMatrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])  
    I = ImmutableMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  
    R = {I: 'I', A: 'A'}
    n = 3
    compute_R_n(R, n, verbose=True)


    I = ImmutableMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  
    A = ImmutableMatrix([[1,0,0], [0,-1, 0], [0,0, 0]])
    omega = -sp.Rational(1, 2) - sp.Rational(1, 2) * sp.sqrt(3) * sp.I
    B = ImmutableMatrix([[0, 1-omega,1-omega**2],[1-omega**2,0,1-omega],[1-omega,1-omega**2,0]])
    R = {I: 'I', A: 'A', B: 'B'}
    n = 3
    compute_R_n(R, n, verbose=True)
        

    B = ImmutableMatrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]])  
    C = ImmutableMatrix([[0, 0, 0], [0, 0, 1], [1, 0, 0]])  
    R = {B: 'B', C:'C'}
    n = 2
    compute_R_n(R, n, verbose=True)
        


    A = ImmutableMatrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]])  
    B = ImmutableMatrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]])  
    C = ImmutableMatrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) 
    R = {A: 'A', B: 'B', C: 'C'}
    n = 2
    compute_R_n(R, n, verbose=True)
                

    A = ImmutableMatrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) 
    B = ImmutableMatrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]])  
    C = ImmutableMatrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) 
    R = {A: 'A', B: 'B', C: 'C'}
    n = 2
    compute_R_n(R, n, verbose=True)


    I = ImmutableMatrix([[1, 0], [0, 1]])  # Identity
    A = ImmutableMatrix([[0, 1], [1, 0]])  # X
    B = ImmutableMatrix([[0, 0], [0, 2]]) 
    R = {I: 'I', A: 'A', B: 'B'}
    n = 3
    compute_R_n(R, n, verbose=True)

    
    I = ImmutableMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  # Identity
    B_tilda = ImmutableMatrix([[0, 0, 0], [0, 2*sp.pi/3, 0], [0, 0, 4*sp.pi/3]]) 
    # compute A_tilda
    omega = sp.symbols('omega')
    A1 = ImmutableMatrix([[1, omega**2, omega], [1, omega, omega**2], [1,1,1]])
    A2 = ImmutableMatrix([[0,0,0], [1, 2*sp.pi/3, 0], [0, 0, 4*sp.pi/3]])
    A3 = ImmutableMatrix([[1, 1,1], [omega, omega**2,1], [omega**2, omega,1]])
    A_tilda = sp.Rational(1, 3) * ImmutableMatrix(reduce(lambda a, b: a * b, [A1, A2, A3]))
    R = {I: 'I', A_tilda: 'Ã', B_tilda: 'B̃'}
    n = 3
    compute_R_n(R, n, verbose=True, output_path=get_data_file_path("R_n_computation.txt"))


if __name__ == "__main__":
    main()