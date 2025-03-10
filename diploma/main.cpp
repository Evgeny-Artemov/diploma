#include </home/evgenii/diploma/lib/space.h>
#include <string>

void printEphemeris(const SEphemeris& ephemeris, const std::string& objectName) {
    std::cout << std::fixed << std::setprecision(6); // Форматирование вывода
    std::cout << "Space Object Name: " << objectName << std::endl;
    std::cout << "Julian Date: " << ephemeris.julianDate << std::endl;

    std::cout << "Position (X, Y, Z): (" 
              << ephemeris.X << ", " 
              << ephemeris.Y << ", " 
              << ephemeris.Z << ")" << std::endl;
    std::cout << "Velocity (Vx, Vy, Vz): (" 
              << ephemeris.Vx << ", " 
              << ephemeris.Vy << ", " 
              << ephemeris.Vz << ")" << std::endl;
}

int main() {
    SpaceObject earth("/home/evgenii/diploma/earth_data.txt", "Earth");
    SpaceObject pandora("/home/evgenii/diploma/pandora.txt", "Pandora");

    std::cout << "Минимальное значение времени " << static_cast<int>(earth.GetMinMaxJulian().first) << std::endl;
    std::cout << "Максимальное значение времени " << static_cast<int>(pandora.GetMinMaxJulian().second) << std::endl;
    
    int julianDateRequest;
    std::cout << "Введите число типа double: "; 
    std::cin >> julianDateRequest;

    auto earthEphemeris = earth.GetEphemeris(julianDateRequest);
    auto pandoraEphemeris = pandora.GetEphemeris(julianDateRequest);

    printEphemeris(earthEphemeris, earth.ObjectName);
    printEphemeris(pandoraEphemeris, pandora.ObjectName);
    return 0;
}