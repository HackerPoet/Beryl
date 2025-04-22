# Beryl
Beryl is a polyhedron solver designed to find embeddings of regular maps.
But it can also be used to find polyhedral embeddings from any list of triangles or at least attempt to minimize the intersections.

For a more detailed explanation about this topic, check out:
* [YouTube Video](todo)
* [Research Paper](https://www.mdpi.com/2073-8994/17/4/622)
* [Interactive Viewer](http://codeparade.net/embeddings/)

## Installation

The only required dependency is the Eigen library which can be downloaded [here](https://eigen.tuxfamily.org/index.php?title=Main_Page).
It is a header-only library, just extract it in the source code directory.

#### Windows
Use Visual Studio 2019 or above to open the project `Beryl.vcxproj` and compile the project.

#### Mac
Open the terminal in this directory and run:
```
clang++ *.cpp -arch x86_64 -o Beryl -Ieigen-3.4.0 -std=c++17
```

#### Linux
Open the terminal in this directory and run:
```
g++ *.cpp -O2 -o Beryl -Ieigen-3.4.0 -std=c++17
```

## Finding Regular Map Embeddings

### Symmetries
Almost all regular maps with genus less than 100 can be found in the 'Maps' folder. The files come in 3 parts:

* A list of triangles
* A list of automorphisms (optional)
* A list of geometric symmetries (for reference)

While the list of automorphisms are technically optional, they are needed to find symmetric embeddings. You will almost always want to use symmetries when possible because they speed up the solve time significantly.

The files ending with `*_symmetries.txt` were generated automatically for all regular maps using the following command. You do not need to run this again unless you're adding more maps.
```
./Beryl -es
```

If you just want to print the embeddable symmetries to the console for a specific regular map use `-e`.
The list of symmetries are candidates for embeddings, but you will have to test each one to determine if it can actually be embedded intersection-free.
For example:
```
./Beryl -m=R3.2 -e
```

Basic algorithmic checks are used to remove symmetries that are impossible to embed without intersections.
But if you choose to allow self-intersections, then more symmetries may be possible.
Use `-ev` to print the entire verbose list. The ones with an asterisk guarantee self-intersection.
```
./Beryl -m=R3.2 -ev
```

In this example, let's say we want to try finding an embedding with D3 geometric symmetry.
You should see this one near the end of the list:
```
 D3 : 44,4    3(3b:12) x 2(2b:8 2c:4)
```

To explain the list, the first column is the [Schoenflies notation](https://en.wikipedia.org/wiki/Schoenflies_notation) for the geometric symmetry. The next column is the list of permutation indexes in the automorphism file. So in our example, lines 44 and 4 of `R3.2_automorphisms.txt` will provide D3 symmetry. Lastly, the final column is a 'signature' of the embedding. Sometimes, there will be more than one possibility for a specific symmetry, like how the D2 symmetry actually has 5 distinct possibilities in the case of R3.2. Some may have embeddings and some won't, so you may have to check them all.

### Main Solver

To continue with the example, let's try to actually find embeddings of R3.2 with different symmetries.
Using the map name `-m=R3.2` and choosing a symmetry from `R3.2_symmetries.txt`, you can try running one of these commands.
Note that symmetries are always specified by a letter followed by a list of indices:
```
./Beryl -m=R3.2          (This will solve with no symmetry)
./Beryl -m=R3.2 -s=S2    (This will solve with S2 symmetry)
./Beryl -m=R3.2 -s=D1,42 (This will solve with D2 symmetry)
./Beryl -m=R3.2 -s=D44,4 (This will solve with D3 symmetry)
```

Sometimes, there is more than one possibility for a symmetry, so you may have to check each one to determine if it can be embedded.
And sometimes there are multiple distinct embeddings within the same symmetry group.

To find a dual just add -d. For instance:
```
./Beryl -m=R5.1 -d -s=S124
```

This may take several minutes as it tries different random seeds.
It will stop automatically once it finds a solution with 0 intersections.
Results are saved in the 'Tri' or 'Poly' folder depending on if it's a triangulation or dual.

### Refinement Solver
The refinement solver serves two purposes.
First, if a solution is found, it will adjust the shape to be more aesthetically pleasing.
This involves trying to move vertices as far from each other as possible, avoid sharp angles, opening gaps, and so on.

Second, some maps and especially larger ones, may be difficult to find intersection free solutions right away.
Beryl will by default save the best solutions found so far, those that minimize self-intersections.
These can then be fed to the refinement solver to reduce the intersection count and ideally find one that is intersection-free.

You should use the same symmetry group that was used during the main solve.
To run the refinement solver from the previous example, use `-l`.
```
./Beryl -m=R5.1 -d -s=S124 -l=path_to_your_file.obj
```

## Advanced Options
To see a full list of options, use `-h`

### Small Integer Coordinates
This can be enabled during the refinement stage to keep the vertices snapped to an integer grid.
The parameter is a scale factor of the grid. The lower this amount, the smaller the coordinates.
Keep in mind that not all symmetries can use integer coordinates.
Do not use for C3, C5 or higher, S6 or higher, D3, D5 or higher, or any duals.
```
./Beryl -m=R3.1 -s=T50,2,14 -si=2.0 -l=path_to_your_file.obj
```

### Other Exports
To export a 'wireframe' of the mesh, use -w
```
./Beryl -l=input_path.obj -w=output_path.obj
```

To export a 'cut-out' of the mesh, use -c
```
./Beryl -l=input_path.obj -c=output_path.obj
```
