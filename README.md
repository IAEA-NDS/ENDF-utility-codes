## ENDF Utility Codes version

The ENDF utility codes enable the verification of files
in ENDF-5 or 6 format. For a history of updates to the
codes consult the information in the headers of the
source files.

The suite of ENDF utility codes includes:

- CHECKR - Format checking code
- FIZCON - Procedures & simple physics checking code
- PSYCHE - More complicated physics checking code
- INTER - Calculates selected cross sections and integrals
- STANEF - Creates directory, adds tape label & converts numeric fields and convert to binary format

The provided makefile enables the compilation of the source
files. By default the GNU Fortran compiler is used.
Change into the directory containing the files and run `make`.
