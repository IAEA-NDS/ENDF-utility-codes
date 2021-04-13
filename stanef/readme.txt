   ENDF utility "stanef"
* Copyright (C) 2004 Broohaven National Laboratory
* Copyright (C) 2005 International Atomic Energy Agency

1. Create executable
   $ gfortran -std=legacy -o stanef stanef.f

2. Test run
   $ stanef <stanef.inp >stanef.tto
   Results:
	stanef.out   - output file
	stanef.tto   - terminal output
   Standard: (use to compare result)
	stanef.out-0 - output file
	stanef.tto-0 - terminal output
