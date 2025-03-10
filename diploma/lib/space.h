#include <iostream>
#include <utility>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>


struct SEphemeris {
    double julianDate;
    double X, Y, Z;
    double Vx, Vy, Vz;
};

class SpaceObject {

    public: 
        SpaceObject(const std::string& filePath, const std::string& objectName);
        SEphemeris GetEphemeris(const double& date);
        std::pair<double, double> GetMinMaxJulian();

    public:
        std::string ObjectName;
    private:
        std::vector<SEphemeris> Ephemeris;
        
};