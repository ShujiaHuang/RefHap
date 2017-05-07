SingleIndividualHaplotyper(SIH) - Efficient heuristic algorithms for the SIH problem
Version 1.0.0 (15/06/11)
===============================================================================
This package contains implementations of eight different algorithms for to solve 
the Single Individual Haplotyping (SIH) problem. By default it runs **ReFHap**, which is
a very efficient and accurate algorithm to assemble haplotypes from different kinds of
fragment data.

Building SingleIndividualHaplotyper
-----------------------------------

SingleIndividualHaplotyper has been compiled and run successfully
on the standard jdk version 1.6.0. To build the distribution
library sih.jar on a command line environment 
run the following commands in the directory
where RefHap is located:

```
cd SingleIndividualHaplotyper
make all
```

Running SingleIndividualHaplotyper
----------------------------------

The front class of SingleIndividualHaplotyper is mpg.molgen.sih.main.SIH and
on default settings calls the ReFHap algorithm for haplotyping. Basic usage just requires 
a file with the description of the fragments (the format is explained below) and the path for
the output file. The full usage is as follows: 
```
Usage: 

java -cp SIH.jar mpg.molgen.sih.main.SIH <OPTIONS> <INPUT_FILE> <OUTPUT_FILE>

OPTIONS:
		-v FILE			: Text file with genomic coordinates of variants. By default it assumes that 
							coordinates are located in the first column
		-c INT			: Column in the variants file where coordinates are located. Fields in this 
							column must be integers. 
		-a STRING		: Name of the algorithm to run. Options are Refhap, DGS, and FastHare. SHRThree, 
							Speedhap, TwoDMEC, and WMLF are also supported although they are in general less accurate
```

The input file is a text file with one line per fragment. Since fragments are allowed to have gaps, each fragment
must be divided in segments of continuous calls. The first field is the number of such segments. The second field is
the id of the fragment. Two fields are then added for each subfragment, the first is the one based relative start
of the subfragment and the second is a binary string with the calls. For HapCUT users, the format is the same 
as HapCUT. A sample file is provided with this distribution.

The output is a tab delimited text file with the list of blocks. Each block has a header with the 
one based relative position of the first variant, the size of the block and the number of variants phased.
Then, there is one line per variant with the following fields:
- Relative position of the variant
- Allele in the first haplotype
- Allele in the second haplotype

If a file is provided with genomic locations of variants through the -v option, then the relative positions are replaced 
with the positions contained in this file.

Example 1. Simple usage:
```
java -cp SIH.jar mpg.molgen.sih.main.SIH example.frags example1.phase
```

Example 2. Including genomic coordinates
```
java -cp SIH.jar mpg.molgen.sih.main.SIH -v example.allvars -c 2 example.frags example2.phase
```

Example 3. Running an algorithm different than ReFHap
```
java -cp SIH.jar mpg.molgen.sih.main.SIH -a DGS -v example.allvars -c 2 example.frags example3.phase
```

License
-------

SingleIndividualHaplotyper is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
 
SingleIndividualHaplotyper is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SingleIndividualHaplotyper.  If not, see <http://www.gnu.org/licenses/>.

