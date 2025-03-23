#include <iostream>
#include <string>
#include <cmath>

double MakeNonDimensional(double value, double unutValue, const std::string& variableName = "variable") {
    double result = value / unutValue;
    
    std::cerr << "Обезразмеривание: " << variableName << " = " << value 
              << " / " << unutValue << " = " << result 
              << " (безразмерная величина)" << std::endl;
    
    return result;
}

double MakeDimensional(double nonDimValue, double unutValue) {
    double result = nonDimValue * unutValue;
    return result;
}

double VectorLength(double x, double y, double z) {
    return std::sqrt(x*x + y*y + z*z);
}
