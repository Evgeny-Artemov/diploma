#include <iostream>
#include <vector>
#include <cmath> // Для функций exp() и fabs()
#include <boost/numeric/odeint.hpp>
#include <functional> // Для std::ref
#include <iomanip> // Для std::setprecision

using namespace std;
using namespace boost::numeric::odeint;

// Определяем тип для вектора состояния системы
typedef std::vector<double> state_type;

// Функция-наблюдатель для вывода результатов во время интегрирования
void observer(const state_type &x, const double t) {
    cout << "t = " << t << ", Положение = (" << x[0] << ", " << x[1] << ", " << x[2] << ")" << endl;
    cout << "t = " << t << ", Скорость = (" << x[3] << ", " << x[4] << ", " << x[5] << ")" << endl;
    cout << "t = " << t << ", Масса = " << x[6] << endl;
}

int main() {
    // Инициализируем вектор состояния с 14 элементами
    state_type x(14);
    
    // Инициализируем все переменные значением 0.0
    for (int i = 0; i < 14; i++) {
        x[i] = 0.0;
    }
    
    // Выводим описание вектора состояния
    cout << "Представление вектора состояния X:" << endl;
    cout << "x[0] - Координата x" << endl;
    cout << "x[1] - Координата y" << endl;
    cout << "x[2] - Координата z" << endl;
    cout << "x[3] - Скорость Vx" << endl;
    cout << "x[4] - Скорость Vy" << endl;
    cout << "x[5] - Скорость Vz" << endl;
    cout << "x[6] - Масса m" << endl;
    cout << "x[7] - Сопряженная переменная lambdaX (связана с координатой x)" << endl;
    cout << "x[8] - Сопряженная переменная lambdaY (связана с координатой y)" << endl;
    cout << "x[9] - Сопряженная переменная lambdaZ (связана с координатой z)" << endl;
    cout << "x[10] - Сопряженная переменная lambdaVx (связана со скоростью Vx)" << endl;
    cout << "x[11] - Сопряженная переменная lambdaVy (связана со скоростью Vy)" << endl;
    cout << "x[12] - Сопряженная переменная lambdaVz (связана со скоростью Vz)" << endl;
    cout << "x[13] - Сопряженная переменная lambdaM (связана с массой m)" << endl;
    
    cout << "\nРешаем систему дифференциальных уравнений методом Дорманда-Принса" << endl;
    cout << "Начальные условия установлены в 0.0 для всех переменных состояния" << endl;
    cout << "\nРезультаты численного решения:" << endl;

    // Определяем тип интегратора (Рунге-Кутта-Фельберг 7(8))
    typedef runge_kutta_fehlberg78<state_type> stepper_type;

    // Определяем систему дифференциальных уравнений
    auto system = [](const state_type &x, state_type &dxdt, double t) {
        // Физические константы
        const double mu = 132712440018;       // Гравитационный параметр Солнца
        const double P = ;        // Величина силы тяги
        const double m0 = 100.0;     // Начальная масса топлива
        const double W = 1.0;        // Удельный импульс двигательной установки
        
        // Управление тягой (бинарное: 0 или 1), определяется согласно принципу максимума
        double delta = 0.0;  // Изначально отключено (нет тяги)
        
        // Вычисляем направление тяги (единичный вектор)
        // Здесь определяем его на основе значений lambda (сопряженных переменных для скорости)
        double lambda_v_norm = sqrt(x[10]*x[10] + x[11]*x[11] + x[12]*x[12]);
        
        // Единичный вектор направления тяги (если lambda_v_norm не ноль)
        double p0_x = 0.0, p0_y = 0.0, p0_z = 0.0;
        if (lambda_v_norm > 1e-10) {
            p0_x = -x[10] / lambda_v_norm;
            p0_y = -x[11] / lambda_v_norm;
            p0_z = -x[12] / lambda_v_norm;
            
            // Применяем принцип максимума: если функция переключения < 0, тяга включена
            double switching_function = 1.0 - (P * lambda_v_norm) / (W * x[13]);
            delta = (switching_function < 0) ? 1.0 : 0.0;
        }
        
        // Вычисляем величину радиус-вектора (r = |r|)
        double r_magnitude = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        double r_magnitude_cubed = r_magnitude * r_magnitude * r_magnitude;
        
        // Избегаем деления на ноль
        if (r_magnitude_cubed < 1e-10) r_magnitude_cubed = 1e-10;
        
        // Вычисляем текущую массу
        double current_mass = m0 - x[6];
        if (current_mass < 1e-10) current_mass = 1e-10;  // Избегаем деления на ноль
        
        // Дифференциальные уравнения для переменных состояния
        
        // dr/dt = V
        dxdt[0] = x[3];  // dx/dt = Vx
        dxdt[1] = x[4];  // dy/dt = Vy
        dxdt[2] = x[5];  // dz/dt = Vz
        
        // dV/dt = -mu * r/|r|^3 + delta * (P/(m0-m)) * p0
        dxdt[3] = -mu * x[0] / r_magnitude_cubed + delta * (P / current_mass) * p0_x;
        dxdt[4] = -mu * x[1] / r_magnitude_cubed + delta * (P / current_mass) * p0_y;
        dxdt[5] = -mu * x[2] / r_magnitude_cubed + delta * (P / current_mass) * p0_z;
        
        // dm/dt = -P * delta / W
        dxdt[6] = -P * delta / W;
        
        // Уравнения для сопряженных переменных (выведены из принципа максимума Понтрягина)
        
        // d(lambda_r)/dt
        dxdt[7] = -x[10] * mu * (1.0/r_magnitude_cubed - 3*x[0]*x[0]/(r_magnitude_cubed*r_magnitude*r_magnitude));
        dxdt[8] = -x[11] * mu * (1.0/r_magnitude_cubed - 3*x[1]*x[1]/(r_magnitude_cubed*r_magnitude*r_magnitude));
        dxdt[9] = -x[12] * mu * (1.0/r_magnitude_cubed - 3*x[2]*x[2]/(r_magnitude_cubed*r_magnitude*r_magnitude));
        
        // d(lambda_V)/dt = -lambda_r
        dxdt[10] = -x[7];
        dxdt[11] = -x[8];
        dxdt[12] = -x[9];
        
        // d(lambda_m)/dt = delta * P * lambda_V * p0 / (m0-m)^2
        dxdt[13] = delta * P * (x[10]*p0_x + x[11]*p0_y + x[12]*p0_z) / (current_mass*current_mass);
    };

    // Создаем контролируемый интегратор с заданными допусками на ошибку
    auto stepper = make_controlled<stepper_type>(1.0e-10, 1.0e-10);

    // Выполняем адаптивное интегрирование с функцией-наблюдателем
    integrate_adaptive(
        stepper,
        system,
        x,
        0.0, // Начальное время
        10.0, // Конечное время
        0.5, // Начальный шаг
        std::function<void(const state_type&, double)>(observer)
    );

    cout << "\nКонечное состояние при t = 10:" << endl;
    cout << "Положение = (" << x[0] << ", " << x[1] << ", " << x[2] << ")" << endl;
    cout << "Скорость = (" << x[3] << ", " << x[4] << ", " << x[5] << ")" << endl;
    cout << "Израсходованная масса = " << x[6] << endl;
    cout << "Сопряженные переменные: " << endl;
    for (int i = 7; i < 14; i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    return 0;
}