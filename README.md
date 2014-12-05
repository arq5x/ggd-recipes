GGD recipes
===========

Work in progress. Be scared. And critical.

This repository houses YAML recipes defining how to find and transform genomics datasets to be installed by the obscenely nascent and incredibly alpha [GGD](https://github.com/arq5x/ggd) toolset.  GGD seeks to be a prototype for demonstrating how genome annotations and datasets that live on diverse websites and repositories could, in principle, be downloaded and installed via a single tool. This repository houses a "cookbook" of
"recipes" that attempt to do exactly that.  It is woefully incomplete and is a toy at this point.

The recipes follow an ontology. For example, one could use GGD to install CpG islands from UCSC for Human build 38.  The command would be:

	# source.species.genomebuild.name
	python ggd.py install ucsc.human.b38.cpg

In this case, the recipe would live in this repository under:

	ucsc/human/b38/cpg.yaml

Below is an example recipe for extracting Alu repeats from UCSC based on the RepeatMasker track. Note that the recipt command can include multiple UNIX commands in order to facilitate possibly compex transformation rules.  In other words, we can do MUCh more than simply downloading data from URLs.  We can convert it into useful formats and/or create "meta" datasets based on conversions, trsnformations, etc.  From https://github.com/arq5x/ggd-recipes/blob/master/ucsc/human/b37/alus.yaml:

	attributes:
	  name: alus
	  version: 0.1

	recipe:
	  full:
	    recipe_type: bash
	    recipe_cmds: 
	      - >
	        echo "chrom\tstart\tend\tdivergence\trep_name\trep_class\trep_family";
	        mysql --user=genome \
	              --host=genome-mysql.cse.ucsc.edu \
	              -A -N -B -D hg19 \
	              -e "SELECT genoName, genoStart, genoEnd, \
	                         milliDiv+milliIns+milliDel, repName, repClass, \
	                         repFamily \
	                  FROM rmsk WHERE repName like 'Alu%'" \
	        | sort -k1,1 -k2,2n
	    recipe_outfiles:
	      - ucsc.human.b37.alus
    


