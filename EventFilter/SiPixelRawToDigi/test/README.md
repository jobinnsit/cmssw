GPU version of RawToDigi.

1. It produces CPU output of RawToDigi for a single event -> R2D_CPU.txt

2. It produces GPU output of RawToDigi for a single event -> R2D_GPU.txt

3. validation_R2D directory contains code for validating CPU and GPU output of RawToDigi.

4. To validate copy above 2 files in validation directory and run
   g++ -std=c++11 r2d_cpu_gpu_validate_using_tuple.cpp

