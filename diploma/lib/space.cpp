#include "space.h"
#include <cstddef>
#include <iostream>
#include <ostream>

double LinearInterpolation(double x0, double y0, double x1, double y1, double x) {
    if (x0 == x1) {
        std::cerr << "Ошибка: x0 и x1 не должны быть равны." << std::endl;
        return 0;
    }
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}

std::pair<int, int> FindInterval(const std::vector<SEphemeris>& ephemeris, double userTimestamp) {
    if (userTimestamp < ephemeris[0].julianDate) {
        std::cerr << "x0 меньше минимального значения в массиве." << std::endl;
        return {-1, -1};
    }
    if (userTimestamp > ephemeris[ephemeris.size() - 1].julianDate) {
        std::cerr << "x0 больше максимального значения в массиве." << std::endl;
        return {-1, -1};
    }

    for (size_t i = 0; i < ephemeris.size() - 1; ++i) {
        if (ephemeris[i].julianDate <= userTimestamp && userTimestamp <= ephemeris[i + 1].julianDate) {
            return {i, i + 1}; 
        }
    }
    std::cerr << "Интервал не найден." << std::endl;
    return {-1, -1};
}

double GetItemByPatern(const std::string& line, const std::string& patern) {
    size_t pos = line.find(patern);

    if (pos == -1) {
        std::cout << "Bad parsing item " << line << std::endl;
        return 0;
    }
    return std::stod(line.substr(pos + 3));
}


double GetJulianDate(const std::string& line) {
    size_t equalPos = line.find('=');
    if (equalPos == std::string::npos) {
        std::cerr << "Некорректный формат строки: символ '=' не найден" << std::endl;
        return 0;
    }

    std::string jdPart = line.substr(0, equalPos);
    double julianDate = std::stod(jdPart);
    return julianDate;
}

SpaceObject::SpaceObject(const std::string& filePath, const std::string& objectName) {
    ObjectName = objectName;
    std::ifstream file(filePath);
    if (!file.is_open()) {
        std::cerr << "Ошибка открытия файла: " << filePath << std::endl;
        return;
    }

    bool stopFlag = false;
    while (true) {
        SEphemeris ephemeris;
        for (size_t i = 0; i < 3; i += 1) {
            std::string line;
            if (!std::getline(file, line)) {
                stopFlag = true;
                break;
            }

            if (i == 0) {
                ephemeris.julianDate = GetJulianDate(line);
            } else if (i == 1) {
                ephemeris.X = GetItemByPatern(line, "X =");
                ephemeris.Y = GetItemByPatern(line, "Y =");
                ephemeris.Z = GetItemByPatern(line, "Z =");
            } else if (i == 2) {
                ephemeris.Vx = GetItemByPatern(line, "VX=");
                ephemeris.Vy = GetItemByPatern(line, "VY=");
                ephemeris.Vz = GetItemByPatern(line, "VZ=");
            }
        };

        if (stopFlag) {
            break;
        }
        Ephemeris.push_back(ephemeris);
    }
    file.close();
};
    
SEphemeris SpaceObject::GetEphemeris(const double& date) {
    SEphemeris ephemeris;
    const auto indexPair = FindInterval(Ephemeris, date);

    const auto& leftEphemeris = Ephemeris[indexPair.first];
    const auto& rightEphemeris =  Ephemeris[indexPair.second];

    ephemeris.julianDate = date;
    ephemeris.X = LinearInterpolation(leftEphemeris.julianDate, leftEphemeris.X, rightEphemeris.julianDate, rightEphemeris.X,  date);
    ephemeris.Y = LinearInterpolation(leftEphemeris.julianDate, leftEphemeris.Y, rightEphemeris.julianDate, rightEphemeris.Y,  date);
    ephemeris.Z = LinearInterpolation(leftEphemeris.julianDate, leftEphemeris.Z, rightEphemeris.julianDate, rightEphemeris.Z,  date);
    ephemeris.Vx = LinearInterpolation(leftEphemeris.julianDate, leftEphemeris.Vx, rightEphemeris.julianDate, rightEphemeris.Vx,  date);
    ephemeris.Vy = LinearInterpolation(leftEphemeris.julianDate, leftEphemeris.Vy, rightEphemeris.julianDate, rightEphemeris.Vy,  date);
    ephemeris.Vz = LinearInterpolation(leftEphemeris.julianDate, leftEphemeris.Vz, rightEphemeris.julianDate, rightEphemeris.Vz,  date);
    return ephemeris;
}

std::pair<double, double> SpaceObject::GetMinMaxJulian() {
    return {Ephemeris[0].julianDate, Ephemeris[Ephemeris.size() - 1].julianDate};
}

