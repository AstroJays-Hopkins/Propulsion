------------------------------------------------------------
Creo Simulate Structure Version P-30-43:spg
Diagnostic Log
Sun Apr  8, 2018   13:47:43
------------------------------------------------------------

Begin Creating Database for Design Study
Sun Apr  8, 2018   13:47:43

Begin Integrated Mode Error Checking
Sun Apr  8, 2018   13:47:43

Begin Generating Elements
Sun Apr  8, 2018   13:47:43

Begin Integrated Mode Error Checking
Sun Apr  8, 2018   13:47:45

Begin Engine Bookkeeping
Sun Apr  8, 2018   13:47:45

Begin Analysis: "P_VAP_30C"
Sun Apr  8, 2018   13:47:45

Using Sparse Solver

Begin Mass Calculation
Sun Apr  8, 2018   13:47:46

Begin P-Loop Pass 1
Sun Apr  8, 2018   13:47:46

Begin Processing Multi-Point Constraints
Sun Apr  8, 2018   13:47:46

Begin Matrix Profile Minimization
Sun Apr  8, 2018   13:47:46

Begin Element Calculations, Pass 1
Sun Apr  8, 2018   13:47:46

Begin Global Matrix Assembly, Pass 1
Sun Apr  8, 2018   13:47:46

Begin Equation Solve, Pass 1
Sun Apr  8, 2018   13:47:46

Number of equations: 24288
Average bandwidth:   247.191
Maximum bandwidth:   1041
Size of global matrix profile (mb): 48.0302
Number of terms in global matrix profile: 6003777     
Minimum recommended solram for direct solver: 4

Size of element file (mb): 23.1019
Maximum element matrix size (kb): 14.64
Average element matrix size (kb): 14.64

NOTES:
Solver RAM allocation can be done with a single parameter,
called solram.  If the Creo Simulate Structure/Thermal
engine is the only memory-intensive application running on
your computer, performance will usually be best if you set
solram equal to half of your machine RAM.  For example,
solram 512 is a good choice for a machine with 1024 MB of RAM.

If you are running other memory-intensive applications on
your computer, decrease the solram allocation accordingly.
For example, set solram to 0.25 times machine RAM if you are
running two large applications at once.  However, you often
can run two large jobs faster one after another than if you
try to run both jobs at once.

The purpose of solram is to reduce the amount of disk I/O.
If you set solram too high, performance will usually suffer,
even on machines with very large RAM, because there will not
be enough machine RAM for other important data.  For
example, Creo Simulate allocates many large, non-solver
memory areas that will cause excessive swapping unless you
leave enough spare machine RAM.  You must also leave enough
RAM for the operating system to do disk caching.  Disk
caching improves filesystem performance by holding file data
in RAM for faster access.  Setting solram to half machine
RAM is usually the best compromise between reducing the
amount of disk I/O, and leaving enough machine RAM for disk
caching and for other data.

If you set solram too low, performance will suffer because
Creo Simulate must transfer data between machine RAM and
disk files many more times than with a larger setting.
For example, performance may degrade significantly if you
set solram to 0.1 times machine RAM or less.  A preferable
minimum is 0.25 times machine RAM.

The available swap space on your machine must be greater
than the maximum memory usage of your job.  The available
disk space on your machine must be greater than the maximum
disk usage of your job.  You can monitor the resource usage
of your job in the log (stt) file.  Your job may fail if
your machine does not have enough available disk space or
swap space, or if the maximum memory usage of your job is
greater than the memory limits set for your operating
system.

Begin Load Calculations
Sun Apr  8, 2018   13:47:47

Begin Post-Processing Calculations, Pass 1
Sun Apr  8, 2018   13:47:47

Begin Convergence Check Pass 1
Sun Apr  8, 2018   13:47:47

Begin Strain Energy Calculations
Sun Apr  8, 2018   13:47:48

Begin P-Loop Pass 2
Sun Apr  8, 2018   13:47:48

Begin Processing Multi-Point Constraints
Sun Apr  8, 2018   13:47:48

Begin Matrix Profile Minimization
Sun Apr  8, 2018   13:47:48

Begin Element Calculations, Pass 2
Sun Apr  8, 2018   13:47:48

Begin Global Matrix Assembly, Pass 2
Sun Apr  8, 2018   13:47:50

Begin Equation Solve, Pass 2
Sun Apr  8, 2018   13:47:50

Number of equations: 87792
Average bandwidth:   403.33
Maximum bandwidth:   1899
Size of global matrix profile (mb): 283.273
Number of terms in global matrix profile: 35409171    
Minimum recommended solram for direct solver: 14

Size of element file (mb): 158.516
Maximum element matrix size (kb): 519.84
Average element matrix size (kb): 100.454

Begin Load Calculations
Sun Apr  8, 2018   13:47:52

Begin Post-Processing Calculations, Pass 2
Sun Apr  8, 2018   13:47:53

Begin Convergence Check Pass 2
Sun Apr  8, 2018   13:47:54

Begin Displacement and Stress Calculation
Sun Apr  8, 2018   13:47:55

Begin Reaction Calculation
Sun Apr  8, 2018   13:47:56

Begin Strain Energy Calculations
Sun Apr  8, 2018   13:47:56

Completed P-Loop
Sun Apr  8, 2018   13:47:56

Completed Analysis: P_VAP_30C
Sun Apr  8, 2018   13:47:57
