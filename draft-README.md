# Subspecies Classification Service

## Overview


The subspecies classification tool assigns the genotype/subtype of a virus, based on the genotype/subtype assignments maintained by the International Committee on Taxonomy of Viruses (ICTV). This tool infers the genotype/subtype for a query sequence from its position within a reference tree (using the [pplacer](https://matsen.fhcrc.org/pplacer) tool with a reference tree and reference alignment, including the query sequence as input, interpretation of the [pplacer](https://matsen.fhcrc.org/pplacer) result is handled by [Cladinator](https://github.com/cmzmasek/forester/blob/master/forester/java/src/org/forester/application/cladinator.java)).



## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

This module provides the following application specfication(s):
* [SubspeciesClassification](app_specs/SubspeciesClassification.md)


## See also



## References

1.  Pplacer
    [https://matsen.fhcrc.org/pplacer/](https://matsen.fhcrc.org/pplacer/)

2.  Cladinator
    [https://github.com/cmzmasek/forester/blob/master/forester/java/src/org/forester/application/cladinator.java](https://github.com/cmzmasek/forester/blob/master/forester/java/src/org/forester/application/cladinator.java)
