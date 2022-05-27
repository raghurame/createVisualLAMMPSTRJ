#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

typedef struct dumpEntries
{
	int id, atomType, ix, iy, iz;
	float x, y, z;
} DUMP_ENTRIES;

typedef struct coordinates
{
	float x, y, z;
} COORDINATES;

typedef struct simBoundary
{
	float xlo, xhi, ylo, yhi, zlo, zhi, xLength, yLength, zLength;
} SIM_BOUNDARY;

int main(int argc, char const *argv[])
{
	FILE *inputDump, *outputModDump;
	inputDump = fopen (argv[1], "r");
	outputModDump = fopen (argv[2], "w");

	char lineString[2000];
	int polymerStart = atoi (argv[3]), polymerEnd = atoi (argv[4]);
	int chk_beginning = 0, chk_nAtoms = 0, chk_atomEntries = 0, chk_boxBounds = 0;
	int nAtoms = 0, isDumpUninitialized = 1;

	DUMP_ENTRIES *dump;
	COORDINATES com;
	SIM_BOUNDARY boxBounds;

	while (fgets (lineString, 2000, inputDump) != NULL)
	{
		if (chk_atomEntries == 0) fprintf(outputModDump, "%s", lineString);

		if (strstr (lineString, "ITEM: TIMESTEP")) chk_beginning = 1;
		if (strstr (lineString, "ITEM: NUMBER OF ATOMS")) chk_nAtoms = 1;
		if (strstr (lineString, "ITEM: ATOMS")) chk_atomEntries = 1;
		if (strstr (lineString, "ITEM: BOX BOUNDS")) chk_boxBounds = 1;

		// Store the simulation box boundary
		if (chk_boxBounds)
		{
			fgets (lineString, 2000, inputDump);
			sscanf (lineString, "%f %f\n", &boxBounds.xlo, &boxBounds.xhi);
			fprintf(outputModDump, "%s", lineString);

			fgets (lineString, 2000, inputDump);
			sscanf (lineString, "%f %f\n", &boxBounds.ylo, &boxBounds.yhi);
			fprintf(outputModDump, "%s", lineString);

			fgets (lineString, 2000, inputDump);
			sscanf (lineString, "%f %f\n", &boxBounds.zlo, &boxBounds.zhi);
			fprintf(outputModDump, "%s", lineString);

			boxBounds.xLength = boxBounds.xhi - boxBounds.xlo;
			boxBounds.yLength = boxBounds.yhi - boxBounds.ylo;
			boxBounds.zLength = boxBounds.zhi - boxBounds.zlo;

			chk_boxBounds = 0;
		}

		// Read the number of atoms in the dump file
		// Allocate memory for 'dump' struct based on the number of atoms
		if (chk_nAtoms)
		{
			fgets (lineString, 2000, inputDump);
			fprintf(outputModDump, "%s", lineString);
			sscanf (lineString, "%d", &nAtoms);
			chk_nAtoms = 0;

			if (isDumpUninitialized)
			{
				dump = (DUMP_ENTRIES *) malloc (nAtoms * sizeof (DUMP_ENTRIES));
				isDumpUninitialized = 0;
			}
		}

		if (chk_atomEntries)
		{
			// Read the dump file and store the values
			for (int i = 0; i < nAtoms; ++i)
			{
				fgets (lineString, 2000, inputDump);
				sscanf (lineString, "%d %d %f %f %f %*f %*f %*f %d %d %d\n", 
					&dump[i].id, 
					&dump[i].atomType, 
					&dump[i].x, 
					&dump[i].y, 
					&dump[i].z, 
					&dump[i].ix, 
					&dump[i].iy, 
					&dump[i].iz);
			}

			// Resetting the center of mass
			com.x = 0; com.y = 0; com.z = 0;

			// Unwrapping all coordinates
			for (int i = 0; i < nAtoms; ++i)
			{
				dump[i].x += dump[i].ix * boxBounds.xLength;
				dump[i].y += dump[i].iy * boxBounds.yLength;
				dump[i].z += dump[i].iz * boxBounds.zLength;
			}

			// Find the center of mass of polymer chain
			for (int i = (polymerStart - 1); i < polymerEnd; ++i)
			{
				com.x += dump[i].x; com.y += dump[i].y; com.z += dump[i].z;
			}

			com.x /= (polymerEnd - polymerStart - 1); com.y /= (polymerEnd - polymerStart - 1); com.z /= (polymerEnd - polymerStart - 1);

			// Moving all atoms so the center of mass of the polymer is (0, 0, 0)
			for (int i = 0; i < nAtoms; ++i)
			{
				dump[i].x -= com.x; dump[i].y -= com.y; dump[i].z -= com.z;
			}

			// Moving the solvent molecules only
			// Iterating through atoms before polymerStart
			for (int i = 0; i < (polymerStart - 1); ++i)
			{
				while (abs (dump[i].x) > (boxBounds.xLength / 2))
				{
					if (dump[i].x > 0) dump[i].x -= boxBounds.xLength;
					else if (dump[i].x < 0) dump[i].x += boxBounds.xLength;
				}
				while (abs (dump[i].y) > (boxBounds.yLength / 2))
				{
					if (dump[i].y > 0) dump[i].y -= boxBounds.yLength;
					else if (dump[i].y < 0) dump[i].y += boxBounds.yLength;
				}
				while (abs (dump[i].z) > (boxBounds.zLength / 2))
				{
					if (dump[i].z > 0) dump[i].z -= boxBounds.zLength;
					else if (dump[i].z < 0) dump[i].z += boxBounds.zLength;
				}
			}

			// Iterating through atoms after polymerEnd
			for (int i = polymerEnd; i < nAtoms; ++i)
			{
				while (abs (dump[i].x) > (boxBounds.xLength / 2))
				{
					if (dump[i].x > 0) dump[i].x -= boxBounds.xLength;
					else if (dump[i].x < 0) dump[i].x += boxBounds.xLength;
				}
				while (abs (dump[i].y) > (boxBounds.yLength / 2))
				{
					if (dump[i].y > 0) dump[i].y -= boxBounds.yLength;
					else if (dump[i].y < 0) dump[i].y += boxBounds.yLength;
				}
				while (abs (dump[i].z) > (boxBounds.zLength / 2))
				{
					if (dump[i].z > 0) dump[i].z -= boxBounds.zLength;
					else if (dump[i].z < 0) dump[i].z += boxBounds.zLength;
				}
			}

			chk_atomEntries = 0;

			// Printing the modified atom entries
			for (int i = 0; i < nAtoms; ++i)
			{
				fprintf(outputModDump, "%d %d %.4f %.4f %.4f 0.0 0.0 0.0 0 0 0 0.0 0.0 0.0 0.0 0.0 0.0\n", dump[i].id, dump[i].atomType, dump[i].x, dump[i].y, dump[i].z);
			}
		}
	}

	fclose (inputDump);
	fclose (outputModDump);
	return 0;
}