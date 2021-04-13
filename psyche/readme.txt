   ENDF utility "psyche"
* Copyright (C) 2004 Broohaven National Laboratory
* Copyright (C) 2005 International Atomic Energy Agency

1. Create executable
   $ gfortran -std=legacy -o psyche psyche.f

2. Test run
   $ psyche <psyche.inp >psyche.tto
   Results:
	psyche.out   - output file
	psyche.tto   - terminal output
   Standard: (use to compare result)
	psyche.out-0 - output file
	psyche.tto-0 - terminal output
