import numpy as np
import textwrap
import csv


def read_csv(file):
	x = []
	z = []
	
	with open(file, "r") as csvfile:
		reader = csv.reader(csvfile)
		for row in reader:
			if len(row) >= 2:
				try:
					x.append(float(row[0]))
					z.append(float(row[1]))
				except ValueError:
					# Skip rows with invalid float conversion
					continue
	npoints = len(x)
	print(max(x))
	return x, z, npoints


x, z, npoints = read_csv("profilePoints.csv")

# Dedent the script to remove unnecessary indentation
# Generate profile C++ pointField entries
point_entries1 = "\n".join([f"        profile[{i}] = point({x[i]}, 0, {z[i]});" for i in range(len(x))])
point_entries2 = "\n".join([f"        profile[{i}] = point({x[i]}, 3.7, {z[i]});" for i in range(len(x))])
# blockmeshdict = textwrap.dedent(blockmeshdict)

input_file = "blockMeshDict"
output_file = "blockMeshDict.modified"

# The edges block you want to insert
blockmesh = f"""\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2306                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

convertToMeters 1.0;

vertices (
    (0 0 -2.5)      // 0
    (104 0 -2.5)          // 1
    (104 3.7 -2.5)     // 2
    (0 3.7 -2.5)     // 3

    (0 0 2.6)          // 4
    (104 0 2.6)       // 5
    (104 3.7 2.6)     // 6
    (0 3.7 2.6)     // 7
    
);

blocks (
    hex (0 1 2 3 4 5 6 7) (800 10 111) simpleGrading (1 1 4)

);

edges #codeStream
{{
    codeInclude
    #{{
        #include "pointField.H"
        #include "mathematicalConstants.H"
    #}};

    code
    #{{
        
        const label nPoints = {npoints};
        os  << "(" << nl << "spline 0 1" << nl;
        pointField profile(nPoints, Zero);
        {point_entries1}
        os << profile << nl;
        
        os << "spline 4 5" << nl;
        {point_entries2}
        os << profile << nl;

        os  << ");" << nl;

    #}};
}};

boundary (
    inlet
    {{
        type patch;
        faces
	(
	  (1 2 6 5)
	 );
    }}
    outlet
    {{
        type patch;
        faces
	(
	  (0 3 7 4)
	);
    }}
    bottom
    {{
        type wall;
        faces
	(
	    (0 1 2 3)
	);
    }}
    top
    {{
        type wall;
        faces (
	    (4 5 6 7)
	);
    }}
    front
    {{
        type cyclic;
	neighbourPatch  back;
        faces (
            (3 2 6 7)
        );
    }}
    back
    {{
        type cyclic;
	neighbourPatch  front;
        faces (
            (0 1 5 4)
        );
    }}
);

mergePatchPairs ();


// ************************************************************************* //
"""

# Dedent the script to remove unnecessary indentation
blockmesh_script = textwrap.dedent(blockmesh)

# Write the SLURM script to a file
with open(f"blockMeshDict.modi", "w") as file:
    file.write(blockmesh_script)
    print(f"blockmesh is written\n")