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
    auto NuNorm = MakeNonDimensional(Constants::Nu, Constants::UnitNu, "Поостоянная солнца ");
    // Определяем тип интегратора (Рунге-Кутта-Фельберг 7(8))
    typedef runge_kutta_fehlberg78<state_type> stepper_type;

    std::cout << "Время " << Constants::UnitT << std::endl;
    std::cout << "NuNorm" << NuNorm << std::endl;
    // Определяем систему дифференциальных уравнений
    double delta = 0.0;

    auto observer = [&delta, &WNorm](const state_type &x, const double t) {

        double currentMass = x[6];

        double lambdaVx = x[10];
        double lambdaVy = x[11];
        double lambdaVz = x[12];
        double lambdaM = x[13];

        double lambdaVNorm = VectorLength(lambdaVx, lambdaVy, lambdaVz);
        double p0_x = 0.0, p0_y = 0.0, p0_z = 0.0;
        if (lambdaVNorm > 1e-10) {
            p0_x = lambdaVx / lambdaVNorm;
            p0_y = lambdaVy / lambdaVNorm;
            p0_z = lambdaVz / lambdaVNorm;

            double switching_function = lambdaVNorm / currentMass - lambdaM / WNorm;
            delta = (switching_function > 0) ? 1.0 : 0.0;
        }

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
    };

    auto system = [&WNorm, &PNorm, &NuNorm, &delta](const state_type &x, state_type &dxdt, double t) {

        const double curentX = x[0];
        const double curentY = x[1];
        const double curentZ = x[2];

        const double curentVx = x[3];
        const double curentVy = x[4];
        const double curentVz = x[5];

        const double currentMass = x[6];

        const double lambdaX = x[7];
        const double lambdaY = x[8];
        const double lambdaZ = x[9];

        const double lambdaVx = x[10];
        const double lambdaVy = x[11];
        const double lambdaVz = x[12];

        const double lambdaM = x[13];

        // Вычисляем направление тяги (единичный вектор)
        const double lambdaVNorm = VectorLength(lambdaVx, lambdaVy, lambdaVz);

        // Единичный вектор направления тяги (если lambda_v_norm не ноль)

        double Lambda0_x = 0.0, Lambda0_y = 0.0, Lambda0_z = 0.0;
        if (lambdaVNorm > 1e-10) {
            Lambda0_x = lambdaVx / lambdaVNorm;
            Lambda0_y = lambdaVy / lambdaVNorm;
            Lambda0_z = lambdaVz / lambdaVNorm;
        }

        cout << "Integrator delta : " << delta << endl;

        // Вычисляем величину радиус-вектора
        const double r = VectorLength(curentX, curentY, curentZ);
        const double r3 = r * r * r;
        const double r5 = r * r * r * r * r;

        // std::cout << "H1 " << VectorLength(Hx, Hy, Hz) << std::endl;

        double dVx_dt = -NuNorm * curentX / r3 + delta * PNorm * Lambda0_x / currentMass;
        double dVy_dt = -NuNorm * curentY / r3 + delta * PNorm * Lambda0_y / currentMass;
        double dVz_dt = -NuNorm * curentZ / r3 + delta * PNorm * Lambda0_z / currentMass;

        double H = lambdaVx * dVx_dt + lambdaVy * dVy_dt + lambdaVz * dVz_dt + lambdaX * curentVx + lambdaY * curentVy + lambdaZ * curentVz + lambdaM * (delta * PNorm) / WNorm;
        std::cout << "H2 " << H << std::endl;

        // Дифференциальные уравнения для переменных состояния

        // dr/dt = V
        dxdt[0] = curentVx;  // dx/dt = VX
        dxdt[1] = curentVy;  // dy/dt = VY
        dxdt[2] = curentVz;  // dz/dt = VZ

        // dV/dt = -mu * r/|r|^3 + delta * P * p0 / current_mass
        dxdt[3] = dVx_dt;
        dxdt[4] = dVy_dt;
        dxdt[5] = dVz_dt;

        // dm/dt = - delta * W / P (расход топлива)
        dxdt[6] = - delta * PNorm / WNorm;

        // Уравнения для сопряженных переменных (выведены из принципа максимума Понтрягина)

        // d(lambda_r)/dt
        dxdt[7] = lambdaVx * NuNorm / r3 - 3 * NuNorm * curentX / r5 * (lambdaVx * curentX + lambdaVy * curentY + lambdaVz * curentZ);
        dxdt[8] = lambdaVy * NuNorm / r3 - 3 * NuNorm * curentY / r5 * (lambdaVx * curentX + lambdaVy * curentY + lambdaVz * curentZ);
        dxdt[9] = lambdaVz * NuNorm / r3 - 3 * NuNorm * curentZ / r5 * (lambdaVx * curentX + lambdaVy * curentY + lambdaVz * curentZ);

        // d(lambda_V)/dt = -lambda_r
        dxdt[10] = -lambdaX;
        dxdt[11] = -lambdaY;
        dxdt[12] = -lambdaZ;

        dxdt[13] = delta * PNorm / (currentMass * currentMass) * lambdaVNorm;
    };

    // Создаем контролируемый интегратор с заданными допусками на ошибку
    auto stepper = make_controlled<stepper_type>(1.0e-20, 1.0e-20);

    // Выполняем адаптивное интегрирование с функцией-наблюдателем
    auto q = integrate_adaptive(
        stepper,
        system,
        x,
        365.0,    // Начальное время (в ЕВ)
        3000.0,  // Конечное время (в ЕВ) - настроено для безразмерных единиц
        0.1,   // Начальный шаг (в ЕВ)
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
