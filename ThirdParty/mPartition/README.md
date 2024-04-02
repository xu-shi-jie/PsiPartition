# mPartition
Model-based method for partitioning alignments

<b>I.	About mPartition</b><br>
mPartition is a model-based program to partition alignments into subsets such that sites in a subset follow the same evolutionary process. mPartition combines both site rate model and substitution model at sites to properly partition alignments. Experiments on both real and simulated datasets showed that mPartition was better than other partitioning methods tested.  Notably, mPartition overcame the pitfall of site rates-based partitioning method that groups all invariant sites into one subset leading to incorrect trees.<br>

Using the mPartition method will enhance the accuracy of maximum likelihood tree inference, especially for multiple loci or whole genome datasets. <br> 
	
<b>II.	Versions</b><br>
The mPartition program is available for Linux operating system.<br>

<b>III.	Setup</b><br>
1.	The program requires several software<br>
&emsp;&emsp;-&emsp;python 2.7<br>
&emsp;&emsp;-&emsp;IQ-TREE (http://www.iqtree.org/)<br>
&emsp;&emsp;-&emsp;TIGER (https://github.com/thulekm/mPartition/blob/master/tiger_original.tgz)<br>

2.	Run "setup.py"<br>
Note that the setup.py program will create a new config.py file.<br>

<b>IV.	Commands</b><br>
<i>python mPartition.py -f [alignment] -t [minimum_length] -mset [set_of_models] </i><br><br>
The partition scheme will be written to the file Results/par.[alignment]

<i>-f</i> &emsp;&emsp;The path to a DNA or protein alignment in Phylip format.<br>
<i>-t</i> &emsp;&emsp;The minimum length of a subset (default = 50) (optional)<br>
<i>-mset</i>&emsp;A set of substitution models separated by ‘,’ to select the best-fit model for a subset (optional)<br>
&emsp;&emsp;&emsp;The default model set for DNA: JC69,F81,HKY,GTR; and for amino acid: LG,JTT,WAG.<br>
<br><i><b>Example:</b> python mPartition.py -f example.phy</i>

 

