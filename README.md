# BfsStructEnum Manual

## What is BfsStructEnum?
BfsStructEnum takes chemical formula in the form of the number of each atom types and (zero or more) mol file of biconnected cyclic substructures from users.
Then, it enumerates chemical compounds without cycle except for substructures defined by users.

## Requirement:
- Install boost library

for ubuntu
```
sudo apt-get install libboost-program-options-dev
```
- Install openbabel library

Please follow this [instruction](https://openbabel.org/wiki/CMake)
- Install libcurl

For ubuntu
```
sudo apt-get install libcurl4-openssl-dev
```

## How to run 
- Run `make` command
- Run "./bfsenum " command follows by the options that you want
    - "-c" + number of carbon atoms
    - "-n" + number of nitrogen atoms
    - "-o" + number of oxygen atoms
    - "-h" + number of hydrogen atoms
    - "-s" + one or more directory of mol files of biconnected cyclic substructures
    - "-t" to write the enumerated compounds in SMILES format to the file named "output.smi", otherwise it returns only the number of enumerated compounds on the screen

### Example 
```
./bfsenum -c6 -h6 -m dir/benzene.mol
```
if the chemical formula is `C6H6` and the desired substructure is in `dir/benzene.mol` without the smi file of the enumerated structures.

## Output explanation
This work enumerates each possible combination of the number of defined substructures per round.
Output of each round consists of two main parts. 
1. the number of nodes in the enumerated compounds and 
2. the accumulated number of enumerated structure until the current round. 
Moreover, the total number of enumerated compounds is concluded in the final line.

## Caution 
If a "double free allocation" error is detected during the run time, please use jemalloc to allocate the memory to avoid such error.
