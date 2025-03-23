#include <iostream>
#include <ostream>
#include <random>
#include <vector>
#include <cmath>
#include <boost/numeric/odeint.hpp>

#include </home/evgenii/diploma/constants.h>
#include </home/evgenii/diploma/util/util.h>

using namespace std;
using namespace boost::numeric::odeint;

// Определяем тип для вектора состояния системы
typedef std::vector<double> state_type;

// Функция-наблюдатель для вывода результатов во время интегрирования
void observer(const state_type &x, const double t) {
    cout << "t = " << MakeDimensional(t, Constants::UnitT) / Constants::UnitT  << ", Положение = (" 
    << MakeDimensional(x[0], Constants::UnitR) << ", "
    << MakeDimensional(x[1], Constants::UnitR)  << ", "
    << MakeDimensional(x[2],  Constants::UnitR) << ")" << endl;

    cout << "t = " << MakeDimensional(t, Constants::UnitT) / Constants::UnitT << ", Скорость = ("
    << MakeDimensional(x[3], Constants::UnitV) << ", "
    << MakeDimensional(x[4], Constants::UnitV) << ", "
    << MakeDimensional(x[5], Constants::UnitV) << ")" << endl;

    cout << "t = " << MakeDimensional(t, Constants::UnitT) / Constants::UnitT << ", Масса = " << MakeDimensional(x[6], Constants::UnitMass) << endl;

    cout << "t = " << MakeDimensional(t, Constants::UnitT) / Constants::UnitT << ", Cопряженные к координатам = " << x[7] << ", " << x[8] << ", " << x[9] << ")" << endl;
    cout << "t = " << MakeDimensional(t, Constants::UnitT) / Constants::UnitT << ", Cопряженные к скорости = "<< x[10] << ", " << x[11] << ", " << x[12] << ")" << endl;
    cout << "t = " << MakeDimensional(t, Constants::UnitT) / Constants::UnitT << ", Cопряженная к массе = " << x[13] << endl;
}


double random_double() {
    static std::random_device rd;  
    static std::mt19937 gen(rd()); 
    static std::uniform_real_distribution<double> dis(-1.0, 1.0);
    return dis(gen);
}


namespace ZeroCondeitions {
    double X = -1.321562505626410E+08;
    double Y = 6.658171852356435E+07;
    double Z = -2.596377182275057E+03;
    double VX = -1.389857140256131E+01;
    double VY = -2.671292464300294E+01;
    double VZ = 1.718741410194013E-03;

    double Mass = 1000;
    double P = 0.75 * 10e-3; // кг * км / с ИД-500, разработанный Исследовательским центром им. М. В. Келдыша,
    double Weffective = 70; // км/с
}

double CalcH(double lambdaRI, double lambdaVI, double lambdaM, double lambdaRINorm, double currentMass, double rI, double r_magnitude,  double vI, double delta, double PNorm, double WNorm, double NuNorm) {
    return -(lambdaVI * NuNorm * rI / (r_magnitude * r_magnitude * r_magnitude)) + delta * PNorm * lambdaRINorm * lambdaVI / currentMass + lambdaRI * vI + lambdaM * delta * PNorm / WNorm  ;

}
int main() {
    // Инициализируем вектор состояния с 14 элементами
    state_type x(14);
    
    // Устанавливаем безразмерные начальные условия
    // Положение в ЕД (Единицах Длины)
    x[0] = MakeNonDimensional(ZeroCondeitions::X, Constants::UnitR, "Координата Х ");         // X
    x[1] = MakeNonDimensional(ZeroCondeitions::Y, Constants::UnitR, "Координата Y ");          // Y
    x[2] = MakeNonDimensional(ZeroCondeitions::Z, Constants::UnitR, "Координата Z ");      // Z
    
    // Скорость в ЕД/ЕВ (Единицах Длины/Единицах Времени)
    x[3] = MakeNonDimensional(ZeroCondeitions::VX, Constants::UnitV, "Скорость Vx ");      // VX
    x[4] = MakeNonDimensional(ZeroCondeitions::VY, Constants::UnitV, "Скорость Vy ");     // VY
    x[5] = MakeNonDimensional(ZeroCondeitions::VZ, Constants::UnitV, "Скорость Vy ");       // VZ
    
    // Начальная масса в ЕМ (Единицах Массы)
    x[6] = MakeNonDimensional(1000, Constants::UnitMass,"Масса"); //начальная масса топлива
    
    // Инициализируем сопряженные переменные
    for (int i = 7; i < 14; i++) {
        x[i] = 0.1;
    }

    auto WNorm = MakeNonDimensional(ZeroCondeitions::Weffective, Constants::UnitV, "Эфективная скорость исчечения газов ");
    auto PNorm = MakeNonDimensional(ZeroCondeitions::P, Constants::UnitP, "Тяга ");
    auto NuNorm = MakeNonDimensional(Constants::Nu, Constants::UnitNu, "К ");
    // Определяем тип интегратора (Рунге-Кутта-Фельберг 7(8))
    typedef runge_kutta_fehlberg78<state_type> stepper_type;

    std::cout << "NuNorm" << NuNorm << std::endl;
    // std::cout << "NuNorm2" << Constants::UnitNu << std::endl;
    // Определяем систему дифференциальных уравнений
    auto system = [&](const state_type &x, state_type &dxdt, double t) {
        double delta = 0.0;  // Изначально отключено (нет тяги)
  
        double curentX = x[0];
        double curentY = x[1];
        double curentZ = x[2];

        double curentVx = x[3];
        double curentVy = x[4];
        double curentVz = x[5];

        double currentMass = x[6];
        if (currentMass < 1e-10) currentMass = 1e-10;  // Избегаем деления на ноль

        double lambdaX = x[7];
        double lambdaY = x[8];
        double lambdaZ = x[9];

        double lambdaVx = x[10];
        double lambdaVy = x[11];
        double lambdaVz = x[12];

        double lambdaM = x[13];
        // Вычисляем направление тяги (единичный вектор)
        double lambdaVNorm = VectorLength(lambdaVx, lambdaVy, lambdaVz);
        
        // Единичный вектор направления тяги (если lambda_v_norm не ноль)
        // auto curent_delta = 
        double p0_x = 0.0, p0_y = 0.0, p0_z = 0.0;
        if (lambdaVNorm > 1e-10) {
            p0_x = lambdaX / lambdaVNorm;
            p0_y = lambdaVy / lambdaVNorm;
            p0_z = lambdaVz / lambdaVNorm;
            
            // Применяем принцип максимума: если функция переключения < 0, тяга включена 
            double switching_function = lambdaVNorm / currentMass + lambdaM / WNorm;
            delta = (switching_function > 0) ? 1.0 : 0.0;
        }
        // Вычисляем величину радиус-вектора
        double r = VectorLength(curentX, curentY, curentY);
        double rCubed = r * r * r;

        // Избегаем деления на ноль
        if (rCubed < 1e-10) rCubed = 1e-10;
        
        // Дифференциальные уравнения для переменных состояния
        
        // dr/dt = V
        dxdt[0] = x[3];  // dx/dt = VX
        dxdt[1] = x[4];  // dy/dt = VY
        dxdt[2] = x[5];  // dz/dt = VZ

        // dV/dt = -mu * r/|r|^3 + delta * P * p0 / current_mass
        dxdt[3] = -NuNorm * x[0] / rCubed + delta * PNorm * p0_x / x[6];
        dxdt[4] = -NuNorm * x[1] / rCubed + delta * PNorm * p0_y / x[6];
        dxdt[5] = -NuNorm * x[2] / rCubed + delta * PNorm * p0_z / x[6];
        
        // dm/dt = - delta * W / P (расход топлива)
        dxdt[6] = -delta * PNorm / WNorm;
        
        // Уравнения для сопряженных переменных (выведены из принципа максимума Понтрягина)
        
        // d(lambda_r)/dt
        dxdt[7] = (NuNorm * x[10] / rCubed) - (3 * x[0] * x[0] * NuNorm * x[10]) / (rCubed * r * r);
        dxdt[8] = (NuNorm* x[11] / rCubed) - (3 * x[1] * x[1] * NuNorm * x[11]) / (rCubed * r * r);
        dxdt[9] = (NuNorm * x[12] / rCubed) - (3 * x[2] * x[2] * NuNorm * x[12]) / (rCubed * r * r);

        // d(lambda_V)/dt = -lambda_r
        dxdt[10] = -x[7];
        dxdt[11] = -x[8];
        dxdt[12] = -x[9];
        
        dxdt[13] = -delta * PNorm / (x[6] * x[6]) * lambdaVNorm;

        auto Hx = CalcH(lambdaX, lambdaVx, lambdaM, p0_x, currentMass, curentX, r, curentVx, delta, PNorm, WNorm, NuNorm);
        auto Hy = CalcH(lambdaY, lambdaVy, lambdaM, p0_y, currentMass, curentY, r, curentVy, delta, PNorm, WNorm, NuNorm);
        auto Hz = CalcH(lambdaZ, lambdaVz, lambdaM, p0_z, currentMass, curentZ, r, curentVz, delta, PNorm, WNorm, NuNorm);

        std::cout << VectorLength(Hx, Hy, Hz) << std::endl;
    };

    // Создаем контролируемый интегратор с заданными допусками на ошибку
    auto stepper = make_controlled<stepper_type>(1.0e-10, 1.0e-10);

    // Выполняем адаптивное интегрирование с функцией-наблюдателем
    auto q = integrate_adaptive(
        stepper,
        system,
        x,
        1.0,    // Начальное время (в ЕВ)
        365.0,  // Конечное время (в ЕВ) - настроено для безразмерных единиц
        1.0,    // Начальный шаг (в ЕВ)
        observer
    );

    cout << "q " << q << endl;
    cout << "\nКонечное состояние:" << endl;
    cout << "Положение = (" << x[0] << ", " << x[1] << ", " << x[2] << ")" << endl;
    cout << "Скорость = (" << x[3] << ", " << x[4] << ", " << x[5] << ")" << endl;
    cout << "Оставшаяся масса топлива = " << x[6] << endl;
    cout << "Сопряженные переменные: " << endl;
    for (int i = 7; i < 14; i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}
