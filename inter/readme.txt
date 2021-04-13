   ENDF utility "inter"
* Copyright (C) 2004 Broohaven National Laboratory
* Copyright (C) 2005 International Atomic Energy Agency

1. Create executable
   $ gfortran -std=legacy -o inter inter.f

2. Test run
   $ inter <inter.inp >inter.tto
   Results:
	inter.lst   - report file
	inter.tto   - terminal output
   Standard: (use to compare result)
	inter.out-0 - report file
	inter.tto-0 - terminal output
