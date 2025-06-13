
# Rⁿ Calculator

Computes R^n by multiplying matrices from a given set R up to n times. The output includes only linearly independent matrices.

---

## Getting Started

### 1. Clone the Repository

```
git clone https://github.com/yourusername/RTG-2025
cd RTG-2025
```

> Replace `yourusername` with your actual GitHub username.

---

### 2. Install Python

Make sure you have Python 3.8+ installed.

Download Python: https://www.python.org/downloads/

Check your Python version:

```
python --version
# or
python3 --version
```

---

### 3. Install Dependencies

Install all required packages using requirements.txt:

```
pip install -r requirements.txt
# or
pip3 install -r requirements.txt
```

---

### 4. Run the Program

Navigate to the `src` directory:

```
cd src
```

Run the program:

```
python computeRn.py
# or
python3 computeRn.py
```

---

## Running Your Own Examples

Here’s an example of how to construct a matrix set and compute R^n:

```python

# Identity matrix
I = ImmutableMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

B_tilda = ImmutableMatrix([[0, 0, 0], [0, 2*sp.pi/3, 0], [0, 0, 4*sp.pi/3]])

# Compute A_tilda
omega = sp.symbols('omega')
A1 = ImmutableMatrix([[1, omega**2, omega], [1, omega, omega**2], [1, 1, 1]])
A2 = ImmutableMatrix([[0, 0, 0], [1, 2*sp.pi/3, 0], [0, 0, 4*sp.pi/3]])
A3 = ImmutableMatrix([[1, 1, 1], [omega, omega**2, 1], [omega**2, omega, 1]])

A_tilda = sp.Rational(1, 3) * ImmutableMatrix(reduce(lambda a, b: a * b, [A1, A2, A3]))

R = {I: 'I', A_tilda: 'Ã', B_tilda: 'B̃'}
n = 3

# Run computation
compute_R_n(R, n, verbose=True, output_path=get_data_file_path("R_n_computation.txt"))
```

- Set `verbose=True` to print detailed output of the computations at each step.
- Set `output_path` to the path where you'd like the output to be saved, or don't include this parameter to print to the console only.