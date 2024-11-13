import FreeCAD as App
import FreeCADGui as Gui
import Part
import math
import numpy as np
from scipy.special import fresnel

# Clothoid Spiral class
class ClothoidSpiral:
    def __init__(self, length, radius, direction, num_points=100):
        self.length = length
        self.radius = radius
        self.direction = direction
        self.num_points = num_points
        self.points = []

    def calculate_points(self, placement=App.Vector(0, 0, 0), angle=0):
        # Calculate Fresnel parameters
        A = math.sqrt(self.length * self.radius * math.pi)
        t_values = np.sqrt(np.linspace(0, 1, self.num_points)) * (self.length / A)

        S, C = fresnel(t_values)

        # Calculate coordinates with the direction and scaling factor
        x_coords = A * C
        y_coords = A * S * self.direction

        # Create a rotation object for the given angle around the Z-axis
        rotation = App.Rotation(App.Vector(0, 0, 1), Radian=angle)

        # Apply translation and rotation to the points
        self.points = [placement + rotation.multVec(App.Vector(x, y, 0))
                       for x, y in zip(x_coords, y_coords)]

    def toShape(self):
        if not self.points:
            return None
        bspline = Part.BSplineCurve()
        bspline.interpolate(self.points)
        return bspline.toShape()

def bearing_angle(v1, v2, axis="x"):
    # Function to calculate the bearing angle between two vectors
    # v1: Start vector
    # v2: End vector
    # axis: Rotation axis ("x" or "y")

    if axis == "x":
	    start = App.Vector(1, 0, 0)
    elif axis == "y":
        start = App.Vector(0, 1, 0)

    direction = v1.sub(v2)
    angle = direction.getAngle(start)
    return 2 * math.pi - angle if direction.y < 0 and direction.x < 0 else angle

# Function to find the angle difference and turn direction
def turn_direction(STin, STpi,STout):
    # Calculate initial and final angles
    angle_in = bearing_angle(STpi, STin, "y")
    angle_out = bearing_angle(STout, STpi, "y")

    deflection = angle_out - angle_in
    direction = 1 if deflection > 0 else -1
    return abs(deflection), direction

# Function to calculate spiral parameters
def spiral_parameters(length, radius):
    p = (length ** 2) / (24 * radius)
    X = length - (length ** 3 / (40 * radius ** 2) - length ** 5 / (3456 * radius ** 4))
    Y = X * math.tan((length * math.pi / 2) / (math.pi * radius) / 3)
    k = X - radius * math.sin((length * math.pi / 2) / (math.pi * radius)) 
    Is = (length * math.pi / 2) / (math.pi * radius) 
    SPchord = math.sqrt(X ** 2 + Y ** 2)
    return p, X, Y, k, Is, SPchord

# Function to create a tangent line
def makeTangent(start, end):
    tangent = Part.LineSegment(start, end)
    return tangent.toShape()

# Function to create a Clothoid spiral
def makeSpiral(length, radius, direction, placement, angle):
    spiral = ClothoidSpiral(length, radius, direction)
    spiral.calculate_points(placement, angle)
    return spiral.toShape(), spiral.points[-1]

# Function to create an arc
def makeCurve(SC, CS, radius, direction):
    # Calculate the midpoint and the distance between the points
    chord_middle = (SC + CS) / 2
    chord_length = SC.distanceToPoint(CS)

    # Calculate the distance from the midpoint to the center
    dist_to_center = math.sqrt(abs(radius**2 - (chord_length / 2)**2))

    # Calculate the vector perpendicular to the line segment between the points
    perp_vector = App.Vector(0,0,1).cross(CS.sub(SC)).normalize().multiply(dist_to_center)

    # Select the correct center based on direction
    center = chord_middle.add(perp_vector) if direction > 0 else chord_middle.sub(perp_vector)
    middle = chord_middle.sub(center).normalize().multiply(radius).add(center)
    curve = Part.Arc(SC, middle, CS)

    return curve.toShape()



# Data
radius = 3000000
length_in = 2000000
length_out = 2000000

previous = App.Vector(451298775, 1433498332, 0)
current = App.Vector(450028351, 1428116943, 0)
next = App.Vector(459061372, 1422718608, 0)

deflection, direction = turn_direction(previous, current, next)

angle_in = bearing_angle(current, previous)
angle_out = bearing_angle(current, next)

# Calculate spiral parameters
p1, X1, Y1, k1, Is1, SP1chord = spiral_parameters(length_in, radius)
p2, X2, Y2, k2, Is2, SP2chord = spiral_parameters(length_out, radius)

# Calculate transition points
T1 = (radius + p2) / math.sin(deflection) - \
     (radius + p1) / math.tan(deflection) + k1

T2 = (radius + p1) / math.sin(deflection) - \
     (radius + p2) / math.tan(deflection) + k2

TS = previous.sub(current).normalize().multiply(T1).add(current)
ST = next.sub(current).normalize().multiply(T2).add(current)

tangent_in = makeTangent(previous, TS)
spiral_in, SC = makeSpiral(length_in, radius, -direction, TS, angle_in)
spiral_out, CS = makeSpiral(length_out, radius, direction, ST, angle_out)
curve = makeCurve(SC, CS, radius, -direction)
tangent_out = makeTangent(ST, next)

shape = Part.makeCompound([tangent_in, spiral_in, curve, spiral_out, tangent_out])

if not Gui.activeDocument():
    App.Console.PrintWarning("FreeCAD GUI is not active.\n")

else:
    doc = FreeCAD.activeDocument()
    # Add the compound object to the document
    alignment = doc.addObject("Part::Feature", "Alignment")
    alignment.Shape = shape

    # Update the document
    doc.recompute()
