Andrea noticed (and convinced Ashi) that SLS is not catchable by a diagonal matrix alone. That is purely unique fragments. We want to find cases also where fragments are available by subtraction. So we want some kind of "triangular matrix" type check, except that triangular matrices are normally abut what is full rather than what is empty.  It may be that the only "operator" to easily do this is one that overwrites columns with zeros to approximate a molecule as having been removed, and doing this all the way done the line. This may be cheaper computationally than doing delete.  There could be a list or array of row indices remaining for mass frag rows (and also molecule columns) such that if it's been cleared it gets ignored for checks like summation etc.

I attached the edited version of the algorithm for the ExtentOfUniqueSolvable(). I realized something about it though. The algorithm we came up with does not account for solvable reference patterns that multiple reference intensities for a mass fragment. It only looks for rows containing 1 mass fragment.

I'm thinking something like this:
           Ethanol    Ethene    Crotyl Alcohol
m31    30            70            20
m32    0              40            80
m33    0              0              100


I think we need to add a deletion for the molecules that have a unique fragment, and then create a recursive function that will take the sum again with the remaining molecules to see if unique fragment occur in the other molecules.
