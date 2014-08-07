GenePy
======

Pronounced *génépi*, like the [French alpine spirit](http://en.wikipedia.org/wiki/G%C3%A9n%C3%A9pi).

**GenePy** is a Python package that acts as an interface between BioPython, ClustalO, and PhyML, all in a neat data structure. 

Dependencies :
--------------

- [The Scipy stack](http://scipy.org) ( NumPy, Matplotlib, pandas )
- [BioPython](http://biopython.org "BioPython")
- [Clustal Omega](http://www.clustal.org/omega)
- [PhyML](https://code.google.com/p/phyml)

Optional :

- [seaborn](https://github.com/mwaskom/seaborn)
- [scikit-learn](http://scikit-learn.org)



Sequence Arrays
===============

The core data structure in **GenePy** is the *sequence array*, or `seqarray` object. It acts as a container for a list of sequences, and has member functions for sequence alignment and display, for trimming sequences, for showing sequence array statistics, and for constructing phylogenetic trees.



Example Usage
-------------

	import genepy

	# Import an array of sequences in FASTA format.
	mysequences = genepy.seqarray("rubellaE1.fasta")

	# Show the sequences.
	mysequences.show()


![.show() - visually display the sequence array](tutorial_images/figure_1.png)

They're not aligned, and one is much longer than the others. Let's dig a little deeper by printing a summary of the sequence array.

	print mysequences

This returns :

	GenePy sequence set :
	-- 31 sequences
	-- Mean length : 987.6 (min 286, max 9777)
	-- C+G content : 0.677
	-- From file : rubellaE1.fasta

We have a wide range of sequence lengths. Most are small, but based on the visualisation, we have one that's probably a full genome. Our next step will be to align them and cut away parts of sequences that we don't have enough duplicates for.


	mysequences.align()

We've aligned the sequences with the default **GenePy** arguments to Clustal Omega. A new file was written to disk : `rubellaE1_aligned_genepy.phy`. This file, in PHYLIP file format, acts as a checkpoint for your future work. Let's visualise the result.

	mysequences.show()

![Updated visualisation after sequence alignment](tutorial_images/figure_2.png)

Most of the sequences belong to the E1 glycoprotein gene. Let's trim excess from the left and right, so we're comparing like with like.

	mysequences.trimalignment(left = 8700, right = 9500)
	mysequences.show()

![Updated visualisation after sequence trimming](tutorial_images/figure_3.png)

Our sequence array is now shorter in terms of nucleotides per sequence. We're considering the same parts of the genome for each sequence, so we can start to look at the nucleotide statistics.

	mysequences.stats()

![Sequence statistics](tutorial_images/figure_4.png)

This shows us the average nucleotide content, distributions of nucleotide frequencies across difference sequences, and the frequencies of two-step nucleotide transitions. We see that there's a very high cytosine and guanine content, and that cytosines tend to be followed by more cytosines, although CpG and GpC are also common. The high C+G content is an interesting feature of the rubella virus, which is somewhat of an outlier in this respect amongst ssRNA viruses.

Let's construct a phylogenetic tree from this alignment, once again using





Function Members
----------------

	.show()

A visual representation of the sequences in the `seqarray`. 

***





	.stats()

Display :
- the nucleotide content
- the distributions of nucleotide frequencies over each sequence as estimated kernel densities ( if scikit-learn is installed ) or as a histogram ( if scikit-learn is not found )
- the two-step transition matrix between nucleotides ( frequency of having an A going to a C, etc. )

***






	.align(force = True, iter = False, full = False, full_iter = False, auto = True, threads = False)

Align the sequences in the `seqarray` by calling Clustal Omega. 

- `force` : overwrite the filename, if the alignment exists. The filename defaults to the filename of the sequence you passed on creation of the sequence array, without the extension, and with `_aligned_genepy.phy` appended. 
- `iter` :
- `full` :
- `full_iter` :
- `auto` :
- `threads` :






***

	.trimalignment(array = None, left = None, right = None)

Remove a number of nucleotides to the left and right of the sequence array. This is useful when you have aligned a number of sequences of different lengths, and want to consider only a full array, where each sequence has the same length.

- `left` : 
- `right` :





***

	.phylotree(nucleotide_frequency = "empirical", bootstrap = -4, search_algorithm = "BEST")

Construct a phylogenetic tree by calling PhyML. 

- `nucleotide_frequency` :
- `bootstrap` :
- `search_algorithm` :


	



Variable Members
------------

- `.seq` - a Python list of BioPython sequence objects ( from `Bio.Seq.Seq` )
- `.len` - the number of sequences in the sequence array
- `.seq_len` - a Python list of sequence lengths








Current State
=============

At the moment, **GenePy** is in pre-alpha. Much of the documentation is missing, and there's not yet much functionality. I'll be adding functionality, mostly as I require it, but please feel free to send a pull request my way if you'd like to contribute.


Installation
============

Soon, **GenePy** will be installable using *pip* or *easy_install*. Until then, it's still built as a full Python package, so if you're in a directory where you can see the `genepy/` folder ( but not inside it ), you can just call

	import genepy

You can add this directory to your Python path temporarily. I'm hoping to have **GenePy** on the PyPI at some point in the future for easy installation.