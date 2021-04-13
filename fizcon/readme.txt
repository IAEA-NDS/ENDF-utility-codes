   ENDF utility "fizcon"
* Copyright (C) 2004 Broohaven National Laboratory
* Copyright (C) 2005 International Atomic Energy Agency

1. Create executable
   $ gfortran -std=legacy -o fizcon fizcon.f

2. Test run
   $ fizcon <fizcon.inp >fizcon.tto
   Results:
	fizcon.out   - output file
	fizcon.tto   - terminal output
   Standard: (use to compare result)
	fizcon.out-0 - output file
	fizcon.tto-0 - terminal output
