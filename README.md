Radon-V-line-2D

Oh hi! This project contains code for a 2D V-line Radon transform inversion. The code is written in MATLAB with some C++ attachments. It's simply a step along the way to a full 3D radon transform inversion.

Running the code

Ensure you have a mex compiler set up for MATLAB to run the C++ code attachments used in the line integral function. See the MATLAB website for information.
Before running any script, type mex SLLineIntegralApproximator.cpp in the command window.
Run the (appropriately named) script called Main_Script.m.
Results

You'll see two plots: an image of a phantom and a reconstructed image. The reconstructed image has heavy artifacts due to the lack of a filter. In normal linear filtered backprojection algorithms, a Ram-Lak filter is used. We don't have the luxury here, so it's still being researched.

General Information

The C++ code was last compiled on the 21st September 2017 using the Apple LLVM version 9.0.0 (clang-900.0.37) compiler. The filter used for backprojection is still remaining, but there is some useful code included here, particularly the numerical line integrator.

Thank you for your attention, and I hope you have a wonderful day!
