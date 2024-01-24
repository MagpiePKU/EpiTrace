Init_Matrix
------------

Input: 

- cellname: a string vector denoting the names for each column in the read count matrix. 
- peakname: a string vector denoting the names for each row in the read count matrix. 
- matrix: a read count matrix, where columns are cells, rows are peaks, and elements are "reads in a particular peak for a particular cell". 

Output: a clean-up and prepared read count matrix. 