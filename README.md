# Complex Unit Circle Visualization

An interactive visualization of complex numbers and their relationship to the unit circle.
This educational tool helps understand the connection between complex numbers, trigonometric functions, and polar form.

## Features

- Interactive unit circle visualization
- Real-time conversion between Cartesian and polar forms
- Visual representation of complex number operations
- Educational tool for learning complex numbers

## Preview
![image](https://github.com/user-attachments/assets/e7f3ad34-6836-4160-be2b-092bd658a63c)
![image](https://github.com/user-attachments/assets/cbe8e421-254c-47aa-8c18-d0c732ee17eb)



## Installation

### Prerequisites
- Python 3.6+
- NumPy
- Matplotlib
- PyInstaller (for building executable only)

### Setup
```bash
# Clone the repository
git clone https://github.com/MaNafromSaar/complex_unit_circle

# Install dependencies
pip install -r requirements.txt
```

## Usage

Run the main script:
```bash
python complex_unit_circle.py
```

## Building Standalone Executable

A standalone executable can be built using the included `build_executable.py` script:

```bash
# Install required dependencies first
pip install -r requirements.txt

# Run the build script
python build_executable.py
```

This will create a standalone executable file named `ComplexUnitCircle.exe` in the `dist` folder. The executable can be run on any compatible Windows system without requiring Python or any dependencies to be installed.

## Project Evolution

This project evolved from several small command line programs located in the root directory:
- `angle_to_complex.py` - Converts angles to complex numbers
- `complex_to_angle.py` - Converts complex numbers to angles

These simpler programs were the starting point that led to the development of this comprehensive visualization tool.

## Educational Use

This software is free to use for educational purposes under the MIT License. 
Feel free to use it in classrooms, workshops, or online courses.
Only exception are schools in the Saarland state of Germany. Please reach out via GitHub for a license.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Links

- [GitHub Repository](https://github.com/MaNafromSaar/complex_unit_circle)
- [Khan Academy Forum Discussion](https://www.khanacademy.org/math/precalculus/x9e81a4f98389efdf:complex/x9e81a4f98389efdf:complex-abs-angle/a/complex-number-absolute-value-and-angle-review)
