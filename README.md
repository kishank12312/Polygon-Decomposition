# Decomposition of ploygon into convex polygons

Implementation of the algorithms in [Algorithms for the decomposition of a polygon into convex polygons](https://www.sciencedirect.com/science/article/abs/pii/S0377221799000338)

We implement the algorithm to decompose a given polygon into convex polygons, followed by the merging algorithm, which merges any excess convex polygons obtained in the first step.

### Usage

 - Input polygon should be in a text file, and the format of the input is as follows:
	 - First line contains a positive integer n, the number of nodes in the input polygon
	 - The next n lines contains two space separated numbers denoting the X and Y coordinates of the points in ***clockwise order***
 - Run `main.cpp`
	 - Usage: `g++ main.cpp -o main.exe` followed by `./main.exe <input_file>` where `<input_file>` is the name of the text file with the input polygon, in the above format.
 - The output after running the first decomposition algorithm is stored in `Unmerged Polygon.txt`. The output after merging extra polygons (final output) is stored in `Merged Polygon.txt`
 - Run `visual.py` to visualize both of these outputs in order.

The timing analysis (Time of execution vs Size of input polygon) is present in Timing Analysis.png

Further analysis on select test cases now available in Report.html