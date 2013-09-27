toxicogenomic
=============

Usage of the Code

Step 1: Save "12 DBPs and other chemicals.xlsx" as "12 DBPs and other
chemicals.csv";

Step 2: Put the csv format dataset file in the same folder with Python file
"preprocessing.py";

Step 3: Use command "python preprocessing.py" in Terminal or 
"run preprocessing.py" in iPython prompt to generate a file .pkl file.
Note that it's better to use iPython for it can save variables in the 
workspace. In the following we assume the command is typed in iPython.

Step 4: Use command "run main_hierarchical.py --order 5 --linkage complete" to
run hierarchical clustering.
Note that the user can specify the order of the auto-regressive model and the
linkage method in hierarchical clustering 
