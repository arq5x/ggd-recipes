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


    


