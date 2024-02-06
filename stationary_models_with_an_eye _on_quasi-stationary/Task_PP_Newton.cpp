#include "Task_PP_Newton.h"

// Конструктор класса
Task_PP_Newton::Task_PP_Newton(const Pipiline_parameters& pipiline_parameters_PP_Newton,
	const Oil_parameters& oil_parameters_PP_Newton) {
	setter(pipiline_parameters_PP_Newton, oil_parameters_PP_Newton);

};

// Методы класса
void Task_PP_Newton::setter(const Pipiline_parameters& pipiline_parameters_PP_Newton,
	const Oil_parameters& oil_parameters_PP_Newton) {
	m_pipiline_parameters_PP_Newton = pipiline_parameters_PP_Newton;
	m_oil_parameters_PP_Newton = oil_parameters_PP_Newton;
}



double Task_PP_Newton::solver_newton_rafson() {
	// Задание настроек решателя по умолчанию
	fixed_solver_parameters_t<1, 0> parameters;
	// Создание структуры для записи результатов расчета
	fixed_solver_result_t<1> result;
	// Решение системы нелинейныйх уравнений <2> с помощью решателя Ньютона - Рафсона
	// m_initial_speed_approximation - Начальное приближение
	m_initial_speed_approximation = 0;
	Task_PP_Newton task_PP_Newton(m_pipiline_parameters_PP_Newton, m_oil_parameters_PP_Newton);
	fixed_newton_raphson<1>::solve_dense(task_PP_Newton, { m_initial_speed_approximation }, parameters, &result);
	// Объявляем объект task_PP_Newton класса Bernoulli_equation
	Bernoulli_equation Bernoulli_PP_Newton(m_pipiline_parameters_PP_Newton, m_oil_parameters_PP_Newton,
		0, result.argument, m_internal_diameter);
	// Q - объемный расход [м^3/ч]
	double Q = Bernoulli_PP_Newton.volume_flow();
	return Q;
}