#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <functional>
#include <iomanip>

using namespace std;
using namespace boost::numeric::odeint;

// Определяем тип для вектора состояния системы
typedef std::vector<double> state_type;

// Функция-наблюдатель для вывода результатов во время интегрирования
void observer(const state_type &x, const double t) {
    cout << "t = " << t << ", Положение = (" << x[0] << ", " << x[1] << ", " << x[2] << ")" << endl;
    cout << "t = " << t << ", Скорость = (" << x[3] << ", " << x[4] << ", " << x[5] << ")" << endl;
    cout << "t = " << t << ", Масса = " << x[6] << endl;
    cout << "t = " << t << ", Cопряженные к координатам = " << x[7] << ", " << x[8] << ", " << x[9] << ")" << endl;
    cout << "t = " << t << ", Cопряженные к скорости = "<< x[10] << ", " << x[11] << ", " << x[12] << ")" << endl;
    cout << "t = " << t << ", Cопряженная к массе = " << x[13] << endl;
}


double random_double() {
    static std::random_device rd;  
    static std::mt19937 gen(rd()); 
    static std::uniform_real_distribution<double> dis(-1.0, 1.0);
    return dis(gen);
}

int main() {
    // Инициализируем вектор состояния с 14 элементами
    state_type x(14);
    
    // Устанавливаем безразмерные начальные условия
    // Положение в ЕД (Единицах Длины)
    x[0] = -0.883409;         // X
    x[1] = 0.444803;          // Y
    x[2] = -1.736253e-5;      // Z
    
    // Скорость в ЕД/ЕВ (Единицах Длины/Единицах Времени)
    x[3] = -1.477834e-5;      // VX
    x[4] = -2.840275e-5;      // VY
    x[5] = 1.827423e-9;       // VZ
    
    // Начальная масса в ЕМ (Единицах Массы)
    x[6] = 0;               // Сколько истратили
    
    // Инициализируем сопряженные переменные
    for (int i = 7; i < 14; i++) {
        x[i] = random_double();
    }
    
    // Выводим описание вектора состояния
    cout << "Описание вектора состояния X:" << endl;
    cout << "x[0] - Координата X (ЕД)" << endl;
    cout << "x[1] - Координата Y (ЕД)" << endl;
    cout << "x[2] - Координата Z (ЕД)" << endl;
    cout << "x[3] - Скорость VX (ЕД/ЕВ)" << endl;
    cout << "x[4] - Скорость VY (ЕД/ЕВ)" << endl;
    cout << "x[5] - Скорость VZ (ЕД/ЕВ)" << endl;
    cout << "x[6] - Масса m (ЕМ)" << endl;
    cout << "x[7] - Сопряженная переменная lambdaX (связана с координатой X)" << endl;
    cout << "x[8] - Сопряженная переменная lambdaY (связана с координатой Y)" << endl;
    cout << "x[9] - Сопряженная переменная lambdaZ (связана с координатой Z)" << endl;
    cout << "x[10] - Сопряженная переменная lambdaVx (связана со скоростью VX)" << endl;
    cout << "x[11] - Сопряженная переменная lambdaVy (связана со скоростью VY)" << endl;
    cout << "x[12] - Сопряженная переменная lambdaVz (связана со скоростью VZ)" << endl;
    cout << "x[13] - Сопряженная переменная lambdaM (связана с массой m)" << endl;
    
    cout << "\nРешаем систему дифференциальных уравнений методом Дорманда-Принса" << endl;
    cout << "Начальные условия заданы в безразмерных единицах" << endl;
    cout << "\nРезультаты численного решения:" << endl;

    // Определяем тип интегратора (Рунге-Кутта-Фельберг 7(8))
    typedef runge_kutta_fehlberg78<state_type> stepper_type;

    // Определяем систему дифференциальных уравнений
    auto system = [](const state_type &x, state_type &dxdt, double t) {
        // Физические константы в безразмерных единицах
        const double mu = 1.0;             // Гравитационный параметр (безразмерный)
        const double P = 1.39e-11;         // Тяга (безразмерная)
        const double W = 0.0531;         // Скорость истечения газов (безразмерная)
        const double M0 = 1; 
        // Управление тягой (бинарное: 0 или 1), определяется по принципу максимума
        double delta = 0.0;  // Изначально отключено (нет тяги)
        
        // Вычисляем направление тяги (единичный вектор)
        double lambda_v_norm = sqrt(x[10]*x[10] + x[11]*x[11] + x[12]*x[12]);
        
        // Единичный вектор направления тяги (если lambda_v_norm не ноль)

        // Вычисляем текущую массу
        double current_mass = M0 - x[6];
        if (current_mass < 1e-10) current_mass = 1e-10;  // Избегаем деления на ноль
        
        double p0_x = 0.0, p0_y = 0.0, p0_z = 0.0;
        if (lambda_v_norm > 1e-10) {
            p0_x = -x[10] / lambda_v_norm;
            p0_y = -x[11] / lambda_v_norm;
            p0_z = -x[12] / lambda_v_norm;
            
            // Применяем принцип максимума: если функция переключения < 0, тяга включена
            double switching_function = lambda_v_norm / current_mass - x[13] / W;
            delta = (switching_function > 0) ? 1.0 : 0.0;
        }
        // Вычисляем величину радиус-вектора
        double r_magnitude = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        double r_magnitude_cubed = r_magnitude * r_magnitude * r_magnitude;
        
        // Избегаем деления на ноль
        if (r_magnitude_cubed < 1e-10) r_magnitude_cubed = 1e-10;
        
        // Дифференциальные уравнения для переменных состояния
        
        // dr/dt = V
        dxdt[0] = x[3];  // dx/dt = VX
        dxdt[1] = x[4];  // dy/dt = VY
        dxdt[2] = x[5];  // dz/dt = VZ
        
        // dV/dt = -mu * r/|r|^3 + delta * W * p0
        dxdt[3] = -mu * x[0] / r_magnitude_cubed + delta * W * p0_x;
        dxdt[4] = -mu * x[1] / r_magnitude_cubed + delta * W * p0_y;
        dxdt[5] = -mu * x[2] / r_magnitude_cubed + delta * W * p0_z;
        
        // dm/dt = -delta * W / P (расход топлива)
        dxdt[6] = -delta * W / P;
        
        // Уравнения для сопряженных переменных (выведены из принципа максимума Понтрягина)
        
        // d(lambda_r)/dt
        dxdt[7] = mu * x[10] * (3*x[0]*x[0]/(r_magnitude_cubed*r_magnitude*r_magnitude) - 1.0/r_magnitude_cubed);
        dxdt[8] = mu * x[11] * (3*x[1]*x[1]/(r_magnitude_cubed*r_magnitude*r_magnitude) - 1.0/r_magnitude_cubed);
        dxdt[9] = mu * x[12] * (3*x[2]*x[2]/(r_magnitude_cubed*r_magnitude*r_magnitude) - 1.0/r_magnitude_cubed);
        
        // d(lambda_V)/dt = -lambda_r
        dxdt[10] = -x[7];
        dxdt[11] = -x[8];
        dxdt[12] = -x[9];
        
        dxdt[13] = delta * P / (current_mass * current_mass) * lambda_v_norm;
    };

    // Создаем контролируемый интегратор с заданными допусками на ошибку
    auto stepper = make_controlled<stepper_type>(1.0e-10, 1.0e-10);

    // Выполняем адаптивное интегрирование с функцией-наблюдателем
    integrate_adaptive(
        stepper,
        system,
        x,
        1.0,    // Начальное время (в ЕВ)
        100.0,  // Конечное время (в ЕВ) - настроено для безразмерных единиц
        0.1,    // Начальный шаг (в ЕВ)
        observer
    );

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
