GLOBAL DNA SEQUENCE ALIGNMENT TOOL (AFFINE GAP)

Student: Kunal Singh Gangwar
Course: Bioinformatics / Sequence Analysis Assignment

---------------------------------------------
HOW TO RUN
---------------------------------------------

1) Install Python 3
2) Install numpy:

   pip install numpy

3) Run the program:

   python align.py seq1.fasta seq2.fasta \
   --match 2 --mismatch -1 --gap_open -5 --gap_extend -1

---------------------------------------------
PROGRAM DESCRIPTION
---------------------------------------------
This program performs GLOBAL DNA sequence alignment
using the Needleman–Wunsch algorithm with affine gap penalties
(Gotoh algorithm).

INPUT:
• Two FASTA files
• User-defined scoring parameters:
  match, mismatch, gap opening, gap extension

OUTPUT:
• Best alignment
• Final alignment score

---------------------------------------------
FILES INCLUDED
---------------------------------------------
align.py      → main program
seq1.fasta    → sample input
seq2.fasta    → sample input
README.txt    → instructions
