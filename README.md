## ENDF Utility Codes

The ENDF utility codes enable the verification of files in ENDF-5 or [ENDF-6](https://nds.iaea.org/public/endf/) format.
For the previous history of updates to the codes consult the information in the headers of the source files.

The suite of ENDF Utility Codes includes:

- CHECKR - Format checking code
- FIZCON - Procedures & simple physics checking code
- PSYCHE - More complicated physics checking code
- INTER - Calculates selected cross sections and integrals
- STANEF - Creates directory, adds tape label & converts numeric fields and convert to binary format

In addition to these utilities, the `ENDF2C` is available via the [PREPRO](https://nds.iaea.org/public/endf/prepro/) Codes.

The provided makefiles enable the compilation of the source files. By default the GNU Fortran compiler is used.
Change into the directory containing the files and run `make`.
