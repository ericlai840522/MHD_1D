This is a code for one dimensional magnetohydrodynamics(MHD) simulation
1. Go to src directory, select Rusanov/HLL/HLLD solver in HLLD.cpp
2. Use "make all" to compile the code
3. Use "./main.o $test" to execute, where $test is the test case you want to execute
   Ex: "./main.o 1"
4. Use "python3 draw.py $test" to plot and save the result
   Ex: "python3 draw.py 1"
