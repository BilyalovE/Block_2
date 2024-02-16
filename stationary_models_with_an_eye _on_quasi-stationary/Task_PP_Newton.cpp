#include "Task_PP_Newton.h"

PP_solver_Newton::PP_solver_Newton(Pipeline_parameters pipeline_parameters_PP_Newton, Oil_parameters oil_parameters_PP_Newton)
{
	m_pipeline_parameters_PP_Newton = pipeline_parameters_PP_Newton;
	m_oil_parameters_PP_Newton = oil_parameters_PP_Newton;
}

double PP_solver_Newton::solver_newton_rafson()
{
	// Задание настроек решателя по умолчанию
	fixed_solver_parameters_t<1, 0> parameters;
	// Создание структуры для записи результатов расчета
	fixed_solver_result_t<1> result;
	// Решение системы нелинейных уравнений <2> с помощью решателя Ньютона - Рафсона
	// m_initial_speed_approximation - Начальное приближение
	m_initial_speed_approximation = 0.2;
	fixed_newton_raphson<1>::solve_dense(*this, { m_initial_speed_approximation }, parameters, &result);
	// Объявляем объект task_PP_Newton класса Bernoulli_equation
	Bernoulli_equation task_PP_Newton(m_pipeline_parameters_PP_Newton, m_oil_parameters_PP_Newton, 0, result.argument, 0.7);
	// Q - объемный расход[м ^ 3 / ч]
	//cout << result.argument << endl;
	double Q = task_PP_Newton.volume_flow();
	//cout << Q*3600 << endl;
	return Q;
}

