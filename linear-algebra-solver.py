import numpy as np
from fractions import Fraction
import sympy as sp
from sympy import Matrix, symbols, solve
import json
import os
from pathlib import Path

class LinearAlgebraSolver:
    """
    Linear algebra solver that reads problems from JSON files
    Handles: Gaussian elimination, homogeneous systems, parametric solutions
    """
    
    def __init__(self):
        self.steps = []  # Store solution steps for learning
        self.current_dir = Path(__file__).parent if '__file__' in globals() else Path.cwd()
    
    def print_matrix(self, matrix, title="Matrix"):
        """Pretty print matrix with fractions"""
        print(f"\n{title}:")
        if isinstance(matrix, np.ndarray):
            for row in matrix:
                formatted_row = [str(Fraction(x).limit_denominator()) for x in row]
                print("[" + "  ".join(f"{x:>8}" for x in formatted_row) + "]")
        else:  # sympy matrix
            print(matrix)
    
    def load_from_json(self, filename):
        """
        Load matrix problem from JSON file
        Expected format:
        {
            "problem_name": "Example 1.3.6",
            "problem_type": "homogeneous",  // or "system" or "parametric"
            "matrix_A": [[1, -3, 0, 2, 2], [-2, 6, 1, 2, -5], ...],
            "vector_b": [7, 1, 4],  // optional, for non-homogeneous systems
            "description": "Find basic solutions..."
        }
        """
        filepath = self.current_dir / filename
        
        try:
            with open(filepath, 'r') as file:
                data = json.load(file)
            
            print(f"📁 Loaded: {data.get('problem_name', 'Unknown Problem')}")
            print(f"📝 Description: {data.get('description', 'No description')}")
            print(f"🔧 Type: {data.get('problem_type', 'unknown')}")
            
            # Check if matrix contains symbolic expressions (strings)
            matrix_data = data['matrix_A']
            has_symbols = any(isinstance(element, str) for row in matrix_data for element in row)
            
            if has_symbols:
                print("⚠️  Warning: Matrix contains symbolic expressions (strings)")
                print("   This requires symbolic computation with SymPy")
                print("   Converting to numerical approximation for basic analysis...")
                
                # Try to convert symbolic expressions to numbers for basic analysis
                # This is a simplified approach - real parametric analysis needs SymPy
                A = self._convert_symbolic_matrix(matrix_data)
            else:
                A = np.array(matrix_data, dtype=float)
            
            b = None
            if 'vector_b' in data and data['vector_b'] is not None:
                vector_data = data['vector_b']
                has_symbols_b = any(isinstance(element, str) for element in vector_data)
                
                if has_symbols_b:
                    print("⚠️  Warning: Vector b contains symbolic expressions")
                    b = self._convert_symbolic_vector(vector_data)
                else:
                    b = np.array(vector_data, dtype=float)
            
            return A, b, data
            
        except FileNotFoundError:
            print(f"❌ Error: File '{filename}' not found in {self.current_dir}")
            return None, None, None
        except json.JSONDecodeError:
            print(f"❌ Error: Invalid JSON format in '{filename}'")
            return None, None, None
        except Exception as e:
            print(f"❌ Error loading file: {e}")
            return None, None, None
    
    def _convert_symbolic_matrix(self, matrix_data):
        """Convert matrix with symbolic expressions to numerical matrix"""
        try:
            # For demonstration, replace common symbolic expressions with numerical values
            # In practice, you'd want to use SymPy for proper symbolic computation
            converted = []
            for row in matrix_data:
                converted_row = []
                for element in row:
                    if isinstance(element, str):
                        # Simple symbolic substitutions for demo (k=1)
                        expr = element.replace('k', '1')
                        expr = expr.replace('^', '**')  # Python exponentiation
                        try:
                            value = eval(expr)  # Caution: eval can be dangerous
                            converted_row.append(float(value))
                        except:
                            converted_row.append(0.0)  # Default fallback
                    else:
                        converted_row.append(float(element))
                converted.append(converted_row)
            
            print("   Substituted k=1 for numerical analysis")
            return np.array(converted, dtype=float)
            
        except Exception as e:
            print(f"   Error converting symbolic matrix: {e}")
            # Return a default matrix if conversion fails
            return np.eye(len(matrix_data))
    
    def _convert_symbolic_vector(self, vector_data):
        """Convert vector with symbolic expressions to numerical vector"""
        try:
            converted = []
            for element in vector_data:
                if isinstance(element, str):
                    # Simple symbolic substitutions for demo (k=1)
                    expr = element.replace('k', '1')
                    expr = expr.replace('^', '**')
                    try:
                        value = eval(expr)
                        converted.append(float(value))
                    except:
                        converted.append(0.0)
                else:
                    converted.append(float(element))
            
            return np.array(converted, dtype=float)
            
        except Exception as e:
            print(f"   Error converting symbolic vector: {e}")
            return np.zeros(len(vector_data))
    
    def gaussian_elimination(self, A, b=None, show_steps=True):
        """
        Perform Gaussian elimination
        A: coefficient matrix
        b: constant vector (None for homogeneous systems)
        """
        if b is not None:
            # Create augmented matrix for non-homogeneous systems
            augmented = np.column_stack([A, b])
        else:
            # Homogeneous system
            augmented = A.copy()
        
        m, n = augmented.shape
        self.steps = []
        
        if show_steps:
            self.print_matrix(augmented, "Original Matrix")
        
        # Forward elimination
        pivot_row = 0
        for col in range(min(m, n - (1 if b is not None else 0))):
            # Find pivot
            max_row = pivot_row
            for i in range(pivot_row + 1, m):
                if abs(augmented[i, col]) > abs(augmented[max_row, col]):
                    max_row = i
            
            # Skip if column is all zeros
            if abs(augmented[max_row, col]) < 1e-10:
                continue
            
            # Swap rows if needed
            if max_row != pivot_row:
                augmented[[pivot_row, max_row]] = augmented[[max_row, pivot_row]]
                if show_steps:
                    print(f"R{pivot_row+1} ↔ R{max_row+1}")
                    self.print_matrix(augmented, f"After row swap")
            
            # Make pivot = 1
            pivot = augmented[pivot_row, col]
            augmented[pivot_row] = augmented[pivot_row] / pivot
            if show_steps:
                print(f"R{pivot_row+1} → R{pivot_row+1}/{pivot}")
                self.print_matrix(augmented, f"After scaling row {pivot_row+1}")
            
            # Eliminate below pivot
            for i in range(pivot_row + 1, m):
                if abs(augmented[i, col]) > 1e-10:
                    factor = augmented[i, col]
                    augmented[i] = augmented[i] - factor * augmented[pivot_row]
                    if show_steps:
                        print(f"R{i+1} → R{i+1} - {factor}*R{pivot_row+1}")
            
            if show_steps:
                self.print_matrix(augmented, f"After eliminating column {col+1}")
            
            pivot_row += 1
            if pivot_row >= m:
                break
        
        return augmented
    
    def reduced_row_echelon(self, matrix, show_steps=True):
        """Convert to reduced row-echelon form (RREF)"""
        rref_matrix = matrix.copy()
        m, n = rref_matrix.shape
        
        # Find leading 1s and eliminate above them
        for i in range(m-1, -1, -1):
            # Find leading 1 in this row
            leading_col = -1
            for j in range(n):
                if abs(rref_matrix[i, j] - 1) < 1e-10:
                    leading_col = j
                    break
            
            if leading_col == -1:
                continue  # No leading 1 in this row
            
            # Eliminate above the leading 1
            for k in range(i):
                if abs(rref_matrix[k, leading_col]) > 1e-10:
                    factor = rref_matrix[k, leading_col]
                    rref_matrix[k] = rref_matrix[k] - factor * rref_matrix[i]
                    if show_steps:
                        print(f"R{k+1} → R{k+1} - {factor}*R{i+1}")
        
        if show_steps:
            self.print_matrix(rref_matrix, "Reduced Row-Echelon Form")
        
        return rref_matrix
    
    def find_rank(self, matrix):
        """Find rank of matrix"""
        rref = self.gaussian_elimination(matrix, show_steps=False)
        rank = 0
        m, n = rref.shape
        
        for i in range(m):
            # Check if row is non-zero
            if any(abs(rref[i, j]) > 1e-10 for j in range(n)):
                rank += 1
        
        return rank
    
    def solve_system(self, A, b, show_steps=True):
        """
        Solve linear system Ax = b
        Returns solution type and solutions
        """
        print("="*50)
        print("SOLVING LINEAR SYSTEM Ax = b")
        print("="*50)
        
        augmented = self.gaussian_elimination(A, b, show_steps)
        rref = self.reduced_row_echelon(augmented, show_steps)
        
        m, n = A.shape
        rank_A = self.find_rank(A)
        rank_augmented = self.find_rank(augmented)
        
        print(f"\nRank of A: {rank_A}")
        print(f"Rank of [A|b]: {rank_augmented}")
        print(f"Number of variables: {n}")
        
        # Check consistency
        if rank_A != rank_augmented:
            print("\n🚫 SYSTEM IS INCONSISTENT (No solution)")
            return "inconsistent", None
        
        elif rank_A == n:
            print("\n✅ UNIQUE SOLUTION")
            solution = self.back_substitution(rref, show_steps)
            return "unique", solution
        
        else:
            print(f"\n♾️ INFINITE SOLUTIONS ({n - rank_A} free variables)")
            solution = self.parametric_solution(rref, n, show_steps)
            return "infinite", solution
    
    def solve_homogeneous_system(self, A, show_steps=True):
        """
        Solve homogeneous system Ax = 0
        Find basic solutions and general solution
        """
        print("="*50)
        print("SOLVING HOMOGENEOUS SYSTEM Ax = 0")
        print("="*50)
        
        rref = self.gaussian_elimination(A, show_steps=show_steps)
        rref = self.reduced_row_echelon(rref, show_steps)
        
        m, n = A.shape
        rank = self.find_rank(A)
        num_free = n - rank
        
        print(f"\nRank: {rank}")
        print(f"Number of variables: {n}")
        print(f"Number of free variables: {num_free}")
        print(f"Number of basic solutions: {num_free}")
        
        if num_free == 0:
            print("\n✅ ONLY TRIVIAL SOLUTION: x = 0")
            return [np.zeros(n)]
        
        # Find basic solutions
        basic_solutions = self.find_basic_solutions(rref, n, show_steps)
        
        return basic_solutions
    
    def find_basic_solutions(self, rref, n, show_steps=True):
        """Find basic solutions for homogeneous system"""
        m = rref.shape[0]
        
        # Identify leading variables (pivot columns)
        leading_vars = []
        free_vars = []
        
        for i in range(m):
            for j in range(n):
                if abs(rref[i, j] - 1) < 1e-10:
                    leading_vars.append(j)
                    break
        
        free_vars = [j for j in range(n) if j not in leading_vars]
        
        print(f"\nLeading variables (columns): {[x+1 for x in leading_vars]}")
        print(f"Free variables (columns): {[x+1 for x in free_vars]}")
        
        basic_solutions = []
        
        # Create one basic solution for each free variable
        for idx, free_var in enumerate(free_vars):
            solution = np.zeros(n)
            solution[free_var] = 1  # Set this free variable to 1
            
            # Express leading variables in terms of free variables
            for i, leading_var in enumerate(leading_vars):
                if i < m:  # Make sure we don't go out of bounds
                    # From row i: leading_var + (linear combination of free vars) = 0
                    value = 0
                    for free_var_col in free_vars:
                        if free_var_col < rref.shape[1]:
                            value -= rref[i, free_var_col] * solution[free_var_col]
                    solution[leading_var] = value
            
            basic_solutions.append(solution)
            
            if show_steps:
                print(f"\nBasic solution {idx+1} (free variable x_{free_var+1} = 1):")
                self.print_vector(solution)
        
        # Show general solution
        print(f"\n🎯 GENERAL SOLUTION:")
        print("x = ", end="")
        for i, sol in enumerate(basic_solutions):
            if i > 0:
                print(" + ", end="")
            print(f"t_{i+1} * ", end="")
            self.print_vector(sol, inline=True)
        print(f"\nwhere t_1, t_2, ..., t_{len(basic_solutions)} are arbitrary real parameters")
        
        return basic_solutions
    
    def print_vector(self, vector, inline=False):
        """Pretty print vector"""
        if inline:
            formatted = [str(Fraction(x).limit_denominator()) for x in vector]
            print("[" + ", ".join(formatted) + "]", end="")
        else:
            for x in vector:
                frac_str = str(Fraction(x).limit_denominator())
                print(f"[{frac_str:>8}]")
    
    def back_substitution(self, rref, show_steps=True):
        """Solve system in RREF using back substitution"""
        m, n_aug = rref.shape
        n = n_aug - 1  # Number of variables
        solution = np.zeros(n)
        
        if show_steps:
            print("\nBack substitution:")
        
        for i in range(m-1, -1, -1):
            # Find the leading variable in this row
            leading_col = -1
            for j in range(n):
                if abs(rref[i, j] - 1) < 1e-10:
                    leading_col = j
                    break
            
            if leading_col != -1:
                solution[leading_col] = rref[i, -1]  # The constant term
                if show_steps:
                    print(f"x_{leading_col+1} = {Fraction(solution[leading_col]).limit_denominator()}")
        
        return solution
    
    def parametric_solution(self, rref, n, show_steps=True):
        """Find parametric solution for systems with infinite solutions"""
        # Implementation for parametric solutions
        # This is complex and depends on the specific form
        print("\nParametric solution analysis needed...")
        print("Use basic solutions method for homogeneous systems")
        print("Or manual analysis for non-homogeneous systems")
        return None
    
    def solve_from_file(self, filename):
        """
        Solve linear algebra problem from JSON file
        """
        # Load data from file
        A, b, data = self.load_from_json(filename)
        
        if A is None:
            return None
        
        print("\n" + "="*60)
        
        # Solve based on problem type
        problem_type = data.get('problem_type', 'system')
        
        if problem_type == 'homogeneous':
            return self.solve_homogeneous_system(A)
        elif problem_type == 'system':
            return self.solve_system(A, b)
        elif problem_type == 'parametric':
            print("🔍 Parametric system detected. Use symbolic analysis.")
            return self.gaussian_elimination(A, b)
        else:
            print(f"🤔 Unknown problem type '{problem_type}'. Defaulting to system solver.")
            return self.solve_system(A, b)
    
    def create_sample_files(self):
        """Create sample input files for demonstration"""
        
        # Sample JSON file for Example 1.3.6
        sample_json = {
            "problem_name": "Example 1.3.6 - Homogeneous System",
            "problem_type": "homogeneous",
            "description": "Find basic solutions of the homogeneous system",
            "matrix_A": [
                [1, -3, 0, 2, 2],
                [-2, 6, 1, 2, -5],
                [3, -9, -1, 0, 7],
                [-3, 9, 2, 6, -8]
            ]
        }
        
        json_path = self.current_dir / "example_1_3_6.json"
        with open(json_path, 'w') as f:
            json.dump(sample_json, f, indent=4)
        
        # Sample regular system
        system_json = {
            "problem_name": "Regular 3x3 System",
            "problem_type": "system",
            "description": "Solve the linear system Ax = b",
            "matrix_A": [
                [2, 3, -1],
                [1, -1, 2],
                [3, 2, 1]
            ],
            "vector_b": [7, 1, 4]
        }
        
        system_path = self.current_dir / "sample_system.json"
        with open(system_path, 'w') as f:
            json.dump(system_json, f, indent=4)
        
        # Sample parametric system (simplified for numerical analysis)
        parametric_json = {
            "problem_name": "NUMBAS Parametric System (k=1 substitution)",
            "problem_type": "system",
            "description": "Parametric system with k=1 for numerical analysis",
            "matrix_A": [
                [2, 3, 2],
                [-2, 3, -4],
                [-2, -3, -11]
            ],
            "vector_b": [1, 2, -4],
            "note": "This uses k=1 substitution. For full parametric analysis, use SymPy directly"
        }
        
        param_path = self.current_dir / "parametric_system.json"
        with open(param_path, 'w') as f:
            json.dump(parametric_json, f, indent=4)
        
        print(f"✅ Created sample JSON files in {self.current_dir}:")
        print(f"   📄 example_1_3_6.json")
        print(f"   📄 sample_system.json")
        print(f"   📄 parametric_system.json")
    
    def list_available_files(self):
        """List all JSON files in current directory"""
        json_files = list(self.current_dir.glob("*.json"))
        
        print(f"📂 Available JSON files in {self.current_dir}:")
        for f in json_files:
            print(f"   📄 {f.name}")
        
        if not json_files:
            print("   (No .json files found)")
            print("   Use solver.create_sample_files() to create examples")
    
    def batch_solve(self, pattern="*.json"):
        """Solve all JSON files matching the pattern"""
        files = list(self.current_dir.glob(pattern))
        
        if not files:
            print(f"No JSON files found matching pattern: {pattern}")
            return
        
        print(f"🔄 Batch solving {len(files)} files:")
        results = {}
        
        for file_path in files:
            print(f"\n{'='*60}")
            print(f"🔍 Processing: {file_path.name}")
            print('='*60)
            
            try:
                result = self.solve_from_file(file_path.name)
                results[file_path.name] = result
            except Exception as e:
                print(f"❌ Error processing {file_path.name}: {e}")
                results[file_path.name] = None
        
        return results

def main():
    """Main function to demonstrate JSON-only functionality"""
    solver = LinearAlgebraSolver()
    
    print("🎯 LINEAR ALGEBRA SOLVER - JSON INPUT ONLY")
    print("="*60)
    
    # Check for existing files
    solver.list_available_files()
    
    # Create sample files if none exist
    json_files = list(solver.current_dir.glob("*.json"))
    if not json_files:
        print("\n📝 No JSON files found. Creating sample files...")
        solver.create_sample_files()
        print("\n📋 Now available:")
        solver.list_available_files()
    
    # Solve all available JSON files
    print("\n🔄 Solving all JSON files:")
    solver.batch_solve()

# Convenience functions for JSON file input
def solve_file(filename):
    """Quick function to solve from JSON file"""
    solver = LinearAlgebraSolver()
    return solver.solve_from_file(filename)

def create_samples():
    """Quick function to create sample JSON files"""
    solver = LinearAlgebraSolver()
    solver.create_sample_files()

def list_files():
    """Quick function to list available JSON files"""
    solver = LinearAlgebraSolver()
    solver.list_available_files()

def batch_solve_all():
    """Quick function to solve all JSON files"""
    solver = LinearAlgebraSolver()
    return solver.batch_solve()

if __name__ == "__main__":
    main()
    
    print("\n" + "="*70)
    print("📚 HOW TO USE JSON-ONLY SOLVER:")
    print("="*70)
    print("""
    🔧 SETUP:
    1. Create JSON files with your matrix problems
    2. Place them in the same directory as this script
    
    🎯 USAGE:
    # Create sample files
    python3 linear_solver.py
    
    # Or use convenience functions:
    solver = LinearAlgebraSolver()
    solver.create_sample_files()        # Create examples
    solver.list_available_files()       # List JSON files
    solver.solve_from_file("my.json")   # Solve single file
    solver.batch_solve()                # Solve all JSON files
    
    📄 JSON FORMAT:
    {
        "problem_name": "Your Problem Name",
        "problem_type": "homogeneous|system|parametric",
        "description": "Problem description",
        "matrix_A": [[1, 2, 3], [4, 5, 6]],
        "vector_b": [7, 8]  // only for non-homogeneous
    }
    
    🎛️ PROBLEM TYPES:
    - "homogeneous": Ax = 0 (finds basic solutions)
    - "system": Ax = b (finds unique/infinite/no solution)  
    - "parametric": Systems with parameters
    
    📋 REQUIREMENTS:
    - pip install numpy sympy
    """)