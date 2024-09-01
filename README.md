# Simplex Implementation
Implementation of the Revised Simplex Algorithm, including the Two-Phase Method.

## About
This project is a Simplex method implementation developed in **C++**. The solver is designed to process linear programming instances in `.mps` format, being able to solve problems with over 10,000 variables.
- **Input**: The program reads files in the `.mps` format. The path to the `.mps` file must be provided as a command-line argument.
- **Output**:During execution, the solver prints the current objective function value and the infeasibility value every 1,000 iterations, allowing for continuous progress monitoring. Finally, it reports the optimal solution.

## Instances
The [`instancias`](instancias) directory contains `.mps` files used to test and validate the solver. The instances used were provided by [Netlib](https://netlib.org/lp/data/readme).

## Materials Used

The following reference materials were used for the development and implementation of this solver:

- **Books**: "Introduction to linear optimization' by Dimitris Bertsimas and John N. Tsitsiklis. "Linear Programming" by Vašek Chvátal
- **Article**: Robert E. Bixby, (1992) Implementing the Simplex Method: The Initial Basis. ORSA Journal on Computing 4(3):267-284. [http://dx.doi.org/10.1287/ijoc.4.3.267](http://dx.doi.org/10.1287/ijoc.4.3.267).

## Credits
 **[Raquel Patrício](https://github.com/RaquelPM)**: For providing the `.mps` file reader ([`mpsReader.h`](src/mpsReader.h) and [`mpsReader.cpp`](src/mpsReader.cpp) files)









