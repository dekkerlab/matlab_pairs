## PSB Matlab Workshop about HiC

The purpose of this hands-on workshop, is to let students get a hands-on experience exploring raw Hi-C data.

Sequencing of a Hi-C library yields a number of short paired-end reads, that are mapped to a reference genome.
Mapped locations of "left"- and "right"-side of each read are recorded in a so-called pairs file, along with the strand that each side of the read was mapped to.


In our workshop we start with a pre-processed "pairs"-file for chromosome 19 of the HFF cell line, detailed specifications of the pairs-file format can be found [here](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md).


Further details of the proposed analysis can be found in the commented sections of the matlab script file.

### todo:
during the workshop we realized that an introduction slide or two were missing.
a potential slide would demonstrate a piece of a reference genome with a pair of short reads that are "being" mapped onto it.
The slide would also display location of `pos1` and `pos2` and shed some light on the `str1/2`.



