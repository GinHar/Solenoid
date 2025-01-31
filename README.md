# Solenoid
A program in C++ that simulates the magnetic field of a solenoid. To do this a group of coils are placed one on top of the other in the z axis. The coils have the center in the z axis.

![Magnetic Field](https://github.com/GinHar/Solenoid/blob/main/Magnetic_Field.png)

## Solenoid.cpp
The program that calculates the magnetic field is solenoid.cpp where we introduce the radius, the intensity, the longitude, the number of coils and the point where we want to calculate the magnetic field. Inside solenoid.cpp there are the functions that we have to integrate for each axis that have been obtained with the Biot-Savart law. This integrals are calculate using the Gauss method.

## matrix.h
The matrix.h archive allows us to manipulate matrix.

## integ.h
The integ.h archive have some functions that allow us to calculate definite integrals. In this case we are using Int_Gauss with 100 points.

## Gauss_#n.txt
Finally, we have some text archives which we can use if we want to use more or less points in the Int_Gauss function.
