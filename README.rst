=============================
SangerMutantLibraryAnalysis
=============================


Summary
-----------

This is a script for analyzing the distribution of mutations among Sanger-sequenced clones in a mutant library of a protein-coding gene. The distributions of mutations are analyzed and plotted.

The script for this analysis is available `on GitHub`_.

This script was written by `Jesse Bloom`_.


Requirements
---------------

This analysis simply consists of a `Python`_ script. It has been tested with `Python`_ versions 2.6 and 2.7, and probably works with other 2.* versions as well. 

The script requires `scipy`_ and `matplotlib`_. It has been tested with `scipy`_ 0.12.0 and `matplotlib`_ version 1.2.1, but will probably work with other versions as well.


Running the script
-----------------------

The analysis is performed by the script *analyze_library.py*. To run the script, simply go to the current directory and type the command::

    python analyze_library.py

The script will then ask you to enter the names of two input files: the sequence file, and the mutation list file. These are both text files that should have the following format:

    * **The sequence file:** this is simply a FASTA file that contains a single protein-coding gene. This should be the gene that you are sequencing. For example, here is an example of such a file::

        >Aichi68-NP
        atggcgtcccaaggcaccaaacggtcttatgaacagatggaaactgatggggaacgccagaatgcaactgagatcagagcatccgtcgggaagatgattgatggaattggacgattctacatccaaatgtgcactgaacttaaactcagtgattatgaggggcgactgatccagaacagcttaacaatagagagaatggtgctctctgcttttgacgaaagaaggaataaatatctggaagaacatcccagcgcggggaaggatcctaagaaaactggaggacccatatacaagagagtagatagaaagtggatgagggaactcgtcctttatgacaaagaagaaataaggcgaatctggcgccaagccaataatggtgatgatgcaacagctggtctgactcacatgatgatctggcattccaatttgaatgatacaacataccagaggacaagagctcttgttcgcaccggcatggatcccaggatgtgctctctgatgcagggttcgactctccctaggaggtctggagctgcaggcgctgcagtcaaaggagttgggacaatggtgatggagttgataaggatgatcaaacgtgggatcaatgatcggaacttctggagaggtgaaaatggacgaaaaacaaggagtgcttacgagagaatgtgcaacattctcaaaggaaaatttcaaacagctgcacaaagggcaatgatggatcaagtgagagaaagtcggaacccaggaaatgctgagatcgaagatctcatctttctggcacggtctgcactcatattgagagggtcagttgctcacaaatcttgtctgcccgcctgtgtgtatggacctgccgtagccagtggctacgacttcgaaaaagagggatactctttagtgggaatagaccctttcaaactgcttcaaaacagccaagtatacagcctaatcagaccgaacgagaatccagcacacaagagtcagctggtgtggatggcatgcaattctgctgcatttgaagatctaagagtattaagcttcatcagagggaccaaagtatccccaagggggaaactttccactagaggagtacaaattgcttcaaatgaaaacatggatgctatggaatcaagtactcttgaactgagaagcaggtactgggccataagaaccagaagtggaggaaacactaatcaacagagggcctctgcaggtcaaatcagtgtgcaacctgcattttctgtgcaaagaaacctcccatttgacaaaccaaccatcatggcagcattcactgggaatacagagggaagaacatcagacatgagggcagaaattataaggatgatggaaggtgcaaaaccagaagaaatgtccttccaggggcggggagtcttcgagctctcggacgaaagggcagcgaacccgatcgtgccctcttttgacatgagtaatgaaggatcttatttcttcggagacaatgcagaggagtacgacaattaa

    * **The mutation list file:** this is a text file that lists the mutations. The mutations should be numbered in sequential (1, 2, ...) numbering according to the sequence specified in the sequence file. Lines in this file that are empty or begin with the character # are ignored. All other lines should specify a clone and all identified mutations. The clone name should be the first entry on the line, followed by a colon. There is then a comma-delimited list of the mutations. The mutations are indicated as follows:

        * Single nucleotide substitutions are indicated as *G1378T* for mutation of site 1378 from *G* to *T*.

        * Multiple-nucleotide mutations at the same codon should be listed as *AG349GA* or *TCT1003GCC*. List mutations like this if they are sequential (or two mutations separated by a single other mutation) as these are probably mutations of the same codon.

        * Deletions should be listed like this: *delC392* for deletion of the *C* at position 392.

        * Insertions should be listed like this: *insG392* for insertion of a *G* at position 392.

      The script will check that the specified wildtype nucleotides actually match thos e indicated in the sequence file. If they do not, an error will be raised. Note that even if your sequence contains an insertion or deletion, you must ensure that subsequent sites are still numbered according to sequential numbering of the provided sequence.

      Here is an example input file::

        # Sequential 1, 2, ... numbering from start codon of Aichi68-NP coding sequence
        #
        # Sequencing from Jan-25-2013 for clones 1 to 3 for each.
        # Sequencing from Jan-26-2013 for clones 4 to 8 for each.
        WT1round3-1: delCT392, AG583CT
        WT1round3-2: ATG991TCC, AA1249CG, delA1257, TT1390GA
        WT1round3-3: AAA271CGG, TT1285AC, G1378T
        WT1round3-4: AA175GC, CTC1396GAA
        WT1round3-5: TAT232GTA, ACA388TAC, CT1007AC, AAT1492GCA
        WT1round3-6: AAC925GGG, GGA1303TCG
        WT1round3-7: None
        WT1round3-8: GGA553ATC, A696C, TG713AT, CA1208GT
        #
        WT2round3-1: GC1148AA
        WT2round3-2: AC16CT, AG269GC, TTC1390AGA
        WT2round3-3: CCC481AAG, AC1168GT
        WT2round3-4: TTT211GTA, C322G
        WT2round3-6: AAA817TTT, AGA1306CTC, delA1311
        WT2round3-7: AG92CC, ACA436GTC, ATG1441GAA
        WT2round3-8: C397T, TA854CT, CCG952ACA, A1314T, GAG1483CCC
        #
        N334H1round3-1: GT808AC
        N334H1round3-2: None
        N334H1round3-3: GAT379TAA, G662T, CG1415GC
        N334H1round3-4: A1078G
        N334H1round3-5: T825C, AT1111TA, AAT1189GCG, ACA1309CCT
        N334H1round3-6: G241A, AAA550TGG
        N334H1round3-7: TCT1237CCG, TC1286GG
        N334H1round3-8: TG130AC, AG349GA, G391C, AT487TA, AG799CA
        #
        N334H2round3-1: ATG313CGT, GC929CA, GTG985TGA, AAT1492TTG
        N334H2round3-2: A12C, A1171G, CCA1267GTG
        N334H2round3-3: AGA1024CTG
        N334H2round3-4: G6C, GG160TT, CC265TT, GGT376AAC, C512T, TCT1003GCC
        N334H2round3-6: CCC247TTT, TCT784GAG
        N334H2round3-7: ACT67GTG, AGA193CGT, T616A, TAC1153GTA, AGT1174CAC
        N334H2round3-8: AT425GC, CC470AT, AGA1165GTC, GCA1324CAG


Output of the script
----------------------

The script will print some information about the mutation statistics to standard output. It will also create some PDF plot files. For example, running the script using the example sequence file *Aichi68-NP.fasta* and the example mutation list *mutation_list.txt* provided with this script `on GitHub`_ will produce the following information printed to standard output::

   Beginning analysis.

   Enter the name of the FASTA file containing the gene sequence: Aichi68-NP.fasta
   Read a coding sequence of length 1497

   Enter the name of the file containing the list of mutations: mutation_list.txt

   Reading mutations from mutation_list.txt
   Read entries for 30 clones

   Substitutions begin at following positions: 2, 4, 6, 23, 31, 44, 54, 59, 65, 71, 78, 81, 83, 89, 90, 91, 105, 108, 117, 126, 127, 130, 131, 133, 142, 146, 157, 161, 163, 171, 184, 185, 195, 206, 221, 232, 238, 262, 267, 270, 273, 275, 285, 309, 310, 318, 329, 331, 335, 336, 342, 360, 371, 383, 385, 389, 390, 391, 392, 397, 403, 413, 417, 423, 429, 429, 435, 436, 437, 438, 442, 460, 464, 464, 466, 472, 481, 495, 498, 498

   Indels begin at following positions: 131, 419, 437

   Found a total of 80 substitutions out of 14970 codons sequenced (0.0053)

   Here are the fractions with different numbers of nucleotide mutations:
     1 nucleotide mutations: 0.00100
     2 nucleotide mutations: 0.00227
     3 nucleotide mutations: 0.00207

   Here are the fractions of mutation types
     synonymous: 0.00040
     nonsynonymous: 0.00494
     stop codon: 0.00000

    Now creating the output PDF plot files...
    The output PDF file plots have now all been created.

    Script complete.

The produced PDF files are as follows:


.. figure:: mutpositions.pdf
   :alt: mutpositions.pdf

   The file ``mutpositions.pdf``


.. figure:: mutpositions_cumulative.pdf
   :alt: mutpositions_cumulative.pdf

   The file ``mutpositions_cumulative.pdf``



.. _`on GitHub`: https://github.com/jbloom/SangerMutantLibraryAnalysis
.. _`Jesse Bloom`: http://research.fhcrc.org/bloom/en.html
.. _`Python`: http://www.python.org/
.. _`scipy`: http://www.scipy.org/
.. _`matplotlib`: http://matplotlib.org/
