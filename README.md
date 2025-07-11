# Linear Algebra Solver

A comprehensive Python tool for solving linear algebra problems from JSON input files. Designed for educational purposes and coursework verification.

## Features

- 🔧 **Gaussian Elimination** with step-by-step output
- 🏠 **Homogeneous Systems** with basic solution finding
- 🔢 **General Linear Systems** (unique, infinite, no solutions)
- 📊 **Matrix Rank** calculation
- 🔄 **Batch Processing** of multiple problems
- 📱 **JSON-based Input** for organized problem management
- 🎓 **Educational Output** with detailed explanations

## Installation

### Prerequisites

- Python 3.7 or higher
- pip package manager

### Required Packages

```bash
pip install numpy sympy
```

### Installation Steps

1. Clone or download this repository
2. Navigate to the project directory
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
   Or manually:
   ```bash
   pip install numpy sympy
   ```

## Quick Start

1. **Run the solver** to create sample files:

   ```bash
   python3 linear_solver.py
   ```

2. **This creates sample JSON files:**

   - `example_1_3_6.json` - Homogeneous system example
   - `sample_system.json` - Regular linear system
   - `parametric_system.json` - Parametric system example

3. **The solver automatically processes all JSON files** in the directory

## Usage

### Basic Usage

```python
from linear_solver import LinearAlgebraSolver

# Create solver instance
solver = LinearAlgebraSolver()

# Solve single file
result = solver.solve_from_file("my_problem.json")

# Solve all JSON files
results = solver.batch_solve()

# List available files
solver.list_available_files()

# Create sample files
solver.create_sample_files()
```

### Command Line Usage

```bash
# Process all JSON files in directory
python3 linear_solver.py

# Or use convenience functions
python3 -c "from linear_solver import solve_file; solve_file('example_1_3_6.json')"
python3 -c "from linear_solver import list_files; list_files()"
python3 -c "from linear_solver import batch_solve_all; batch_solve_all()"
```

## JSON File Format

### Basic Structure

```json
{
  "problem_name": "Your Problem Name",
  "problem_type": "homogeneous|system|parametric",
  "description": "Problem description",
  "matrix_A": [
    [1, 2, 3],
    [4, 5, 6]
  ],
  "vector_b": [7, 8]
}
```

### Problem Types

#### Homogeneous Systems (`"problem_type": "homogeneous"`)

For systems of the form **Ax = 0**. Finds basic solutions and general solution.

```json
{
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
```

#### Linear Systems (`"problem_type": "system"`)

For systems of the form **Ax = b**. Determines if solution is unique, infinite, or nonexistent.

```json
{
  "problem_name": "3x3 Linear System",
  "problem_type": "system",
  "description": "Solve the linear system Ax = b",
  "matrix_A": [
    [2, 3, -1],
    [1, -1, 2],
    [3, 2, 1]
  ],
  "vector_b": [7, 1, 4]
}
```

#### Parametric Systems (`"problem_type": "parametric"`)

For systems with parameters. Currently handles simple numerical substitution.

```json
{
  "problem_name": "NUMBAS Parametric System",
  "problem_type": "system",
  "description": "Parametric system with k=1 substitution",
  "matrix_A": [
    [2, 3, 2],
    [-2, 3, -4],
    [-2, -3, -11]
  ],
  "vector_b": [1, 2, -4],
  "note": "Uses k=1 for numerical analysis"
}
```

## Examples

### Example Problems

The solver comes with three example problems:

1. **Example 1.3.6** - Homogeneous system with basic solutions
2. **Sample System** - 3×3 linear system with unique solution
3. **Parametric System** - Simplified parametric analysis

### Creating Your Own Problems

Create JSON files for your coursework:

**numbas_n1_3_1.json:**

```json
{
  "problem_name": "NUMBAS N1.3.1 - Gaussian Elimination",
  "problem_type": "system",
  "description": "Step-by-step Gaussian elimination practice",
  "matrix_A": [
    [3, -15, -15],
    [-1, 2, 11],
    [-3, 12, 18]
  ],
  "vector_b": [-9, 0, 6]
}
```

### Output Example

```
📁 Loaded: Example 1.3.6 - Homogeneous System
📝 Description: Find basic solutions of the homogeneous system
🔧 Type: homogeneous

==================================================
SOLVING HOMOGENEOUS SYSTEM Ax = 0
==================================================

Original Matrix:
[       1      -3       0       2       2]
[      -2       6       1       2      -5]
[       3      -9      -1       0       7]
[      -3       9       2       6      -8]

... [step-by-step Gaussian elimination] ...

Rank: 2
Number of variables: 5
Number of free variables: 3
Number of basic solutions: 3

Leading variables (columns): [1, 3]
Free variables (columns): [2, 4, 5]

🎯 GENERAL SOLUTION:
x = t_1 * [3, 1, 0, 0, 0] + t_2 * [-2, 0, -6, 1, 0] + t_3 * [-2, 0, 1, 0, 1]
```

## Features in Detail

### Gaussian Elimination

- Step-by-step row operations
- Automatic pivot selection
- Clear elimination process
- Fraction output for exact results

### Solution Classification

- **Unique Solution**: When rank(A) = rank([A|b]) = n
- **Infinite Solutions**: When rank(A) = rank([A|b]) < n
- **No Solution**: When rank(A) ≠ rank([A|b])

### Homogeneous Systems

- Automatic basic solution finding
- General solution in parametric form
- Linear combination representation

### Educational Output

- Step-by-step explanations
- Clear mathematical notation
- Learning-focused feedback

## Project Structure

```
linear-algebra-solver/
├── linear_solver.py          # Main solver code
├── README.md                 # This file
├── LICENSE                   # MIT License
├── .gitignore               # Git ignore rules
├── requirements.txt         # Python dependencies
├── example_1_3_6.json      # Sample homogeneous system
├── sample_system.json      # Sample linear system
└── parametric_system.json  # Sample parametric system
```

## Educational Use

This tool is designed for:

- **Coursework Verification**: Check your hand calculations
- **Learning Reinforcement**: See step-by-step solutions
- **Problem Organization**: Maintain a library of solved problems
- **Batch Processing**: Solve multiple problems efficiently

### For Students

- Create JSON files for textbook problems
- Verify NUMBAS exercise solutions
- Organize coursework systematically
- Learn through detailed step-by-step output

### For Portfolio Evidence

- Demonstrate systematic problem-solving approach
- Show technical implementation skills
- Organize coursework in professional format
- Prove deep understanding through automation

## Limitations

- **Parametric Systems**: Currently uses simple numerical substitution (k=1)
- **Symbolic Math**: For full symbolic computation, use SymPy directly
- **Large Systems**: Performance may decrease with very large matrices
- **Floating Point**: Uses standard floating-point arithmetic

## Future Enhancements

- [ ] Full symbolic parametric analysis
- [ ] LaTeX output for mathematical notation
- [ ] Graphical visualization of solutions
- [ ] Web interface
- [ ] Support for complex numbers
- [ ] Performance optimization for large systems

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Developed for linear algebra coursework verification
- Inspired by educational needs in SIT292 Linear Algebra
- Based on W. Keith Nicholson's "Linear Algebra with Applications"

## Support

For questions or issues:

1. Check existing examples in the repository
2. Review the JSON format specification
3. Ensure all dependencies are installed
4. Create an issue with your specific problem

---

**Happy Computing!** 🎓📊🔢
