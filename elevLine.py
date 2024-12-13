import math
import matplotlib.pyplot as plt
import numpy as np
import sys

def read_elevation_map(file_path):
    with open(file_path, 'r') as file:
        file.readline()
        height, width = map(int, file.readline().split())
        height_distance, width_distance = map(float, file.readline().split())
        zero_pointX, zero_pointY = map(int, file.readline().split())
        elevation_map = np.zeros((height, width))
        for i in range(height):
            elevation_map[i] = list(map(float, file.readline().split()))
    return elevation_map, height_distance, width_distance, zero_pointX, zero_pointY

def linear_interpolate(a, b, x):
    res = a * x + b * (1 - x)
    return res

def bilinear_interpolate(a, b, c, d, x, y):
    res1 = linear_interpolate(a, b, x) 
    res2 = linear_interpolate(c, d, x)
    return linear_interpolate(res1, res2, y)

def interpolate_elevation(line_coordinates, elevation_map):
    finalAshVals = []
    for i in range(len(line_coordinates)):
        row, col = line_coordinates[i]
        x00 = elevation_map[math.floor(row)][math.floor(col)]
        x01 = elevation_map[math.floor(row)][math.ceil(col)]
        x10 = elevation_map[math.ceil(row)][math.floor(col)]
        x11 = elevation_map[math.ceil(row)][math.ceil(col)]
        x = row - math.floor(row)
        y = col - math.floor(col)
        xVal = bilinear_interpolate(x00, x01, x10, x11, x, y)
        finalAshVals.append(xVal)
    
    return finalAshVals

def interpolate_line(start, end):
    x1, y1 = start
    x2, y2 = end
    numPoints = math.ceil(math.sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2)))
    dx = (x2 - x1) / numPoints 
    dy = (y2 - y1) / numPoints
    line_coordinates = []

    for i in range(numPoints):
        coord = (x1 + i * dx, y1 + i * dy)
        line_coordinates.append(coord)

    return line_coordinates

def interpolate_distances(line_coordinates, zero_point, cellHeight, cellWidth):
    distances = []
    row1, col1 = zero_point
    minIndex = 0
    min = 800000
    for i in range(len(line_coordinates)):
        row2, col2 = line_coordinates[i]
        dist = math.sqrt(pow((row2 - row1) * cellHeight, 2) + pow((col2 - col1) * cellWidth, 2))
        if(dist < min):
            minIndex = i
            min = dist
    for i in range(len(line_coordinates)):
        row2, col2 = line_coordinates[i]
        dist = math.sqrt(pow((row2 - row1) * cellHeight, 2) + pow((col2 - col1) * cellWidth, 2))
        if(i < minIndex):
            dist = dist * -1
        # if row2 - row1 < 0 or col2 - col1 < 0:
            # dist = dist * -1
        distances.append(dist)
    return distances

def main():
    args = []
    if (len(sys.argv)) % 4 != 2 or len(sys.argv) < 6:
        print("Usage: python3 elevLine.py <file_path> <start_row1> <start_col1> <end_row1> <end_col1> <start_row2> <start_col2> <end_row2> <end_col2> ...")
        print("Running with default values: ImageOutput/MergeTest/AshRainBigScale/Terrain_Ash_Step_48000.elv 256 0 256 512 0 256 512 256 0 0 511 511")
        args = ["elevLine.py", "ImageOutput/MergeTest/AshRainBigScale/Terrain_Ash_Step_48000.elv", 256, 0, 256, 512]
    else:
        args = sys.argv
    
    file_path = args[1]

    elevation_map, height_distance, width_distance, zeroX, zeroY = read_elevation_map(file_path)
    zero_point = (zeroX, zeroY)

    num_lines = (len(args) - 2) // 4

    plt.figure(figsize=(8, 6))

    for i in range(num_lines):
        start_row, start_col, end_row, end_col = map(int, args[2 + i * 4: 6 + i * 4])
        start_point = (start_row, start_col)
        end_point = (end_row, end_col)

        line_coordinates = interpolate_line(start_point, end_point)
        elevation_values = interpolate_elevation(line_coordinates, elevation_map)
        distances = interpolate_distances(line_coordinates, zero_point, height_distance, width_distance)

        plt.plot(distances, elevation_values, label=f'Line {i + 1}')

    plt.xlabel('Distance from the plume (m)')
    plt.ylabel('Ash density (kg ⋅ m^-2)')
    plt.title('Ash density profiles as a function of distance from the crater')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
