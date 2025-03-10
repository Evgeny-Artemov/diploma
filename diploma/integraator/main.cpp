#include <iostream>
#include <vector>
#include <cmath> // для exp() и fabs()
#include <boost/numeric/odeint.hpp>
#include <functional> // Для std::ref

using namespace std;
using namespace boost::numeric::odeint;

// Определение типа для состояния системы
typedef std::vector<double> state_type;


// Наблюдатель для вывода результатов
void observer(const state_type &x, const double t) {
    cout << "t = " << t << ", x = " << x[0] << endl;
    cout << "t = " << t << ", x = " << x[1] << endl;
}

int main() {
    state_type x (14); // Создаем вектор состояния размерности 1
    x[0] = 1.0;      // Начальное условие: x(0) = 1.0
    x[1] = 1.0;

    cout << "Решение уравнения dx/dt = -x методом Dormand-Prince" << endl;
    cout << "Начальное условие: x(0) = " << x[0] << endl;
    cout << "Аналитическое решение: x(t) = exp(-t)" << endl;
    cout << "\nРезультаты численного решения:" << endl;
    
    typedef runge_kutta_fehlberg78<state_type> stepper_type;

    auto system = [](const state_type &x, state_type &dxdt, double t) {
        dxdt[0] = -x[1]; // Простое уравнение dx/dt = -x
        dxdt[1] = x[0]; // Простое уравнение dx/dt = x
    };

    // Создаём контроллер с заданной точностью
    auto stepper = make_controlled<stepper_type>(1.0e-10, 1.0e-10);
    
    // Интегрирование с адаптивным шагом и выводом результатов
    // Используем явное приведение функций к нужным типам
    integrate_adaptive(
        stepper,
        system,
        x,
        0.0,
        10.0,          
        0.5,
        std::function<void(const state_type&, double)>(observer)
    );
    
    // Проверка результата

    cout << "Точность 4: " << std::setprecision(10) << endl;;
    cout << "\nПроверка конечного результата:" << endl;
    cout << "Численное решение при t = 10: x = " << x[0] << endl;
    cout << "Аналитическое решение при t = 10: x = " << exp(-10.0) << endl;
    cout << "Относительная ошибка: " << fabs(x[0] - exp(-10.0))/exp(-10.0) << endl;
    
    return 0;
}