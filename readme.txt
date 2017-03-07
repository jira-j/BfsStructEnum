BfsStructEnum Manual

BfsStructEnum takes chemical formula in the form of the number of each atom types and (zero or more) mol file of biconnected cyclic substructures from users.
Then, it enumerates chemical compounds without cycle except for substructures defined by users.

Requirement:
1) boost library
2) openbabel library

How to run 
Run "Make" command
Run "./bfsenum " command follows by the options that you want
    "-c" + number of carbon atoms
    "-n" + number of nitrogen atoms
    "-o" + number of oxygen atoms
    "-h" + number of hydrogen atoms
    "-s" + (one or more) directory of mol files of biconnected cyclic substructures
    "-t" to write the enumerated compounds in SMILES format to the file named "output.smi", otherwise it returns only the number of enumerated compounds on the screen

Example 
Run "./bfsenum -c6 -h6 -m dir/benzene.mol" if the chemical formula is C6H6, desired substructure is the file ``dir/benzene.mol'', and users do not need to know the structures.

Output explanation
This work enumerates each possible combination of the number of defined substructures per round.
Output of each round consists of two main parts. 
1. the number of nodes in the enumerated compounds and 
2. the accumulated number of enumerated structure until the current round. 
Moreover, the total number of enumerated compounds is concluded in the final line.

Caution 
If a "double free allocation" error is detected during the run time, please use jemalloc to allocate the memory to avoid such error.
