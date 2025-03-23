#include <cmath>
#include <string>


namespace Constants {    
    double Nu = 132712440018; // км3/c2

    // 1 астрономическая единица (расстояние от Земли до Солнца)
    double UnitR = 149598261;  // (149 600 000 км)
    
    double UnitT = 84000; // c 1 день

    // Обезразмеренный гравитационный параметр
    // double UnitNu = UnitR * UnitR * UnitR / (UnitT * UnitT); 
    double UnitNu = UnitR * UnitR * UnitR / (UnitT * UnitT);
    // 
    double UnitV = UnitR / UnitT ; // 

    double UnitAcc =  UnitR / (UnitT * UnitT); // км/с2
    
    // Взята масса топлива КА
    double UnitMass = 1000; // кг

    double UnitP = UnitMass * UnitAcc; // кг * км / с2
}