import FreeCAD as App
import Part
import math
import numpy as np
from scipy.special import fresnel

# Data
arc_radius = 3000000
spiral1_length = 2000000
spiral2_length = 2000000

start_station = App.Vector(451298775, 1433498332, 0)
pi_station = App.Vector(450028351, 1428116943, 0)
end_station = App.Vector(459061372, 1422718608, 0)

# Function to create a Clothoid spiral
def makeSpiral(length, radius, direction, placement, rotation):
    spiral = ClothoidSpiral(length, radius, direction)
    spiral.calculate_points(placement, rotation)
    return spiral.toShape()

# Clothoid Spiral class
class ClothoidSpiral:
    def __init__(self, length, radius, direction, num_points=100):
        self.length = length
        self.radius = radius
        self.direction = direction
        self.num_points = num_points
        self.points = []

    def calculate_points(self, placement=App.Vector(0, 0, 0), rotation=0):
        # Calculate Fresnel parameters
        A = np.sqrt(self.length * self.radius * np.pi)
        t_values = np.linspace(0, self.length / A, self.num_points)
        S, C = fresnel(t_values)
        
        # Calculate coordinates with the direction and scaling factor
        x_coords = A * C
        y_coords = A * S * self.direction

        # Create a rotation object for the given angle around the Z-axis
        rotation = App.Rotation(App.Vector(0, 0, 1), Radian=rotation)

        # Apply translation and rotation to the points
        self.points = [placement + rotation.multVec(App.Vector(x, y, 0))
                       for x, y in zip(x_coords, y_coords)]

    def toShape(self):
        if not self.points:
            return None
        bspline = Part.BSplineCurve()
        bspline.interpolate(self.points)
        return bspline.toShape()

# Function to calculate the angle between two vectors
def bearing_angle(v1, v2):
    bearing = math.atan2(v2.x - v1.x, v2.y - v1.y)
    return bearing + 2 * math.pi if bearing < 0 else bearing

# Function to determine the rotation angle between two vectors
def get_rotation(v1, v2):
    direction = v1.sub(v2)
    rotation = math.acos(direction.normalize().dot(App.Vector(1, 0, 0)))
    return 2 * math.pi - rotation if direction.y < 0 else rotation

# Function to find the angle difference and turn direction
def find_turn_direction(angle1, angle2):
    deflection = angle2 - angle1
    direction = 1 if deflection > 0 else -1
    return abs(deflection), direction

# Calculate initial and final angles
angle1 = bearing_angle(start_station, pi_station)
angle2 = bearing_angle(pi_station, end_station)

deflection, direction = find_turn_direction(angle1, angle2)
# Function to calculate spiral parameters
def calculate_spiral_params(length, radius):
    p = (length ** 2) / (24 * radius)
    X = length - (length ** 3 / (40 * radius ** 2) - length ** 5 / (3456 * radius ** 4))
    Y = X * math.tan((length * math.pi / 2) / (math.pi * radius) / 3)
    k = X - radius * math.sin((length * math.pi / 2) / (math.pi * radius)) 
    Is = (length * math.pi / 2) / (math.pi * radius) 
    SPchord = math.sqrt(X ** 2 + Y ** 2)
    return p, X, Y, k, Is, SPchord

# Calculate spiral parameters
p1, X1, Y1, k1, Is1, SP1chord = calculate_spiral_params(spiral1_length, arc_radius)
p2, X2, Y2, k2, Is2, SP2chord = calculate_spiral_params(spiral2_length, arc_radius)

# Calculate transition points
T1 = (arc_radius + p2) / math.sin(deflection) - \
     (arc_radius + p1) / math.tan(deflection) + k1

TS = pi_station + App.Vector(T1 * math.sin(angle1 + math.pi),
                             T1 * math.cos(angle1 + math.pi), 0)

SC = TS + App.Vector(SP1chord * math.sin(angle1 + (Is1 / 3) * direction),
                     SP1chord * math.cos(angle1 + (Is1 / 3) * direction), 0)

T2 = (arc_radius + p1) / math.sin(deflection) - \
     (arc_radius + p2) / math.tan(deflection) + k2

ST = pi_station + App.Vector(T2 * math.sin(angle1 + deflection * direction),
                             T2 * math.cos(angle1 + deflection * direction), 0)

CS = ST - App.Vector(math.sin(angle1 + (deflection * direction) - ((30 / math.pi / arc_radius / 
                             (spiral2_length + 0.0001) * spiral2_length ** 2) * direction)) * spiral2_length,
                     math.cos(angle1 + (deflection * direction) - ((30 / math.pi / arc_radius / 
                             (spiral2_length + 0.0001) * spiral2_length ** 2) * direction)) * spiral2_length, 0)

CS = ST - App.Vector(math.sin(angle1 + (deflection * direction) - (((math.pi / 6) / math.pi / arc_radius / 
                             (spiral2_length + 0.0001) * spiral2_length ** 2) * direction)) * spiral2_length,
                     math.cos(angle1 + (deflection * direction) - (((math.pi / 6) / math.pi / arc_radius / 
                             (spiral2_length + 0.0001) * spiral2_length ** 2) * direction)) * spiral2_length, 0)


# Function to create an arc
def makeCurve(SC, CS, radius, direction):
    # Calculate the midpoint and the distance between the points
    midpoint = (SC + CS) / 2
    dist_between_points = SC.distanceToPoint(CS)

    # Calculate the distance from the midpoint to the center
    dist_to_center = math.sqrt(radius**2 - (dist_between_points / 2)**2)

    # Calculate the vector perpendicular to the line segment between the points
    dx = CS.x - SC.x
    dy = CS.y - SC.y
    perp_vector = App.Vector(-dy, dx, 0).normalize() * dist_to_center

    # Select the correct center based on direction
    center = midpoint + perp_vector if direction > 0 else midpoint - perp_vector
    chord_mid = (SC + CS) / 2
    thirdPt = chord_mid.sub(center).normalize().multiply(radius).add(center)
    curve = Part.Arc(SC, thirdPt, CS)

    return curve.toShape()

# Function to create a tangent line
def makeTangent(start, end):
    tangent = Part.LineSegment(start, end)
    return tangent.toShape()

rotation1 = get_rotation(pi_station, start_station)
rotation2 = get_rotation(pi_station, end_station)

tangent1 = makeTangent(start_station, TS)
spiral1 = makeSpiral(spiral1_length, arc_radius, -direction, TS, rotation1)
curve = makeCurve(SC, CS, arc_radius, -direction)
spiral2 = makeSpiral(spiral2_length, arc_radius, direction, ST, rotation2)
tangent2 = makeTangent(ST, end_station)

shape = Part.makeCompound([tangent1, spiral1, curve, spiral2, tangent2])
doc = FreeCAD.activeDocument()
# Add the compound object to the document
alignment = doc.addObject("Part::Feature", "Alignment")
alignment.Shape = shape

# Update the document
doc.recompute()
