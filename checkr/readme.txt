   ENDF utility "checkr"
* Copyright (C) 2004 Broohaven National Laboratory
* Copyright (C) 2005 International Atomic Energy Agency

1. Create executable
   $ gfortran -std=legacy -o checkr checkr.f

2. Test run
   $ checkr <checkr.inp >checkr.tto
   Results:
	checkr.out   - output file
	checkr.tto   - terminal output
   Standard: (use to compare result)
	checkr.out-0 - output file
	checkr.tto-0 - terminal output
