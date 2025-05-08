import math
import cmath
import fractions

def get_exact_form(value):
    """Return the exact form of common trigonometric values in terms of square roots"""
    # Common angles with exact forms
    special_values = {
        0: "0",
        1: "1",
        -1: "-1",
        0.5: "1/2",
        -0.5: "-1/2",
        math.sqrt(2)/2: "√2/2",
        -math.sqrt(2)/2: "-√2/2",
        math.sqrt(3)/2: "√3/2",
        -math.sqrt(3)/2: "-√3/2",
        1/math.sqrt(2): "1/√2",
        -1/math.sqrt(2): "-1/√2",
        math.sqrt(3)/3: "√3/3",
        -math.sqrt(3)/3: "-√3/3"
    }
    
    # Check for close matches to account for floating-point precision
    for exact_val, exact_form in special_values.items():
        if abs(value - exact_val) < 1e-10:
            return exact_form
    
    # If not a special value, return the decimal form
    return f"{value:.3f}"

mag = float(input("What is the magnitude of your complex number? "))
angle = float(input("What angle does it express? "))
rad = 0
if (input("Is the angle given in degrees or radians? (r for radians) ")) == "r":
    rad = 1

# Convert angle to radians for calculations if it's in degrees
if rad == 0:
    anglerad = angle * (math.pi / 180)
else:
    anglerad = angle

# Calculate the real and imaginary parts
num1 = mag * math.cos(anglerad)
num2 = mag * math.sin(anglerad)

# Check if this is a common angle with exact form
exact_real = get_exact_form(num1/mag) # Normalize by magnitude
exact_imag = get_exact_form(num2/mag)

print(f"The complex number with magnitude {mag} and angle {angle:.3f} {'degrees' if rad == 0 else 'radians'} is: {num1:.3f} + {num2:.3f}i")

# Show exact form if it's a special angle
if exact_real != f"{num1/mag:.3f}" or exact_imag != f"{num2/mag:.3f}":
    if mag == 1:
        print(f"Exact form: {exact_real} + {exact_imag}i")
    else:
        print(f"Exact form: {mag} × ({exact_real} + {exact_imag}i)")