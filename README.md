# SPyC_Writer

A writer for the Thermo Galactic SPC file format, widely used in Thermo 
Scientific GRAMS and other applications, written in Python. 

The only outside library used is numpy. The remainder of libraries used all come
from the python standard library. Other outside libraries are used for testing, 
but the core module does not use them. Use Python 3.10+.

The module has been cross-validated using https://github.com/rohanisaac/spc and 
[spectragryph](https://www.effemm2.de/spectragryph/) by parsing the output file 
and verifying it matches what is expected.

Validation has been done on Y, XY, XYYY, and XYXYXY file types. All instances 
were successfully parsed.

[Wiley's Know It All](https://sciencesolutions.wiley.com/knowitall-spectroscopy-software/) 
was also used for validating. It successfully parsed all formats except XYXYXY, which it 
states is invalid.

## Tests

Tests are written using pytest and verified using https://github.com/rohanisaac/spc. This package appears to not
have been updated in a while and contains a few errors, but it provides a straight forward interface in python for verifying
parsing.

## To-Do

- Widen support for other SPC aspects that are currently listed as unsupported
	- Allow old file format
	- Provide option for not writing as IEEE float
	- etc.
- Troubleshoot why KIA considers XYXYXY formats invalid
- Expand tests
- Expand input validation and provide informative error messages for wrong input

## References

The Thermo Galactic spc-sdk code:
- https://github.com/bz-dev/spc-sdk

The newer format specification document:
- https://www.yumpu.com/en/document/read/40416248/a-brief-guide-to-spc-file-format-and-using-gspcio

The older format specification document (has greater detail in parts than newer doc):
- https://ensembles-eu.metoffice.gov.uk/met-res/aries/technical/GSPC_UDF.PDF

## Changelog

- 2024-12-16 1.0.1
    - adjusted data validation to shape instead of size
- 2023-10-02 1.0.0
    - fixed for Python 3.11
- 2022-12-05 0.2.0
    - published to PyPi
- 2022-06-06 0.1.0
    - initial version
