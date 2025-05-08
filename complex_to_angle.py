import math
import cmath

num = float(input("What is the real part of your complex number? "))
num2 = float(input("What is the imaginary part of your complex number? Enter without i. "))
rad = 0
if (input("Do you want to see the angle in degrees or radians? (r for radians) ")) == "r":
	rad = 1
if input("Is the sought angle confined to a certain range? (y/n) ") == "y":
	low = float(input("What is the lower limit? "))
	high = float(input("What is the upper limit? "))
	if low > high:
		print("The lower limit must be less than the upper limit. ")
		exit()
	if low < 0 and high > 0:
		print("The angle is not confined to a certain range. ")
		exit()
	if low < 0:
		low += 360
	if high < 0:
		high += 360
angle = math.atan2(num2, num)  # Calculate the angle using atan2
if rad == 0:
	# angle = math.degrees(angle)  # Convert to degrees if needed
	angle = angle * (180 / math.pi)  # Convert to degrees
print(f"The angle of the complex number {num} + {num2}i is: {angle:.3f} {'degrees' if rad == 0 else 'radians'}")

