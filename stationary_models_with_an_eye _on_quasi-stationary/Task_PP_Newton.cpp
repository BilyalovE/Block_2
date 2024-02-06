#include "Task_PP_Newton.h"

// ����������� ������
Task_PP_Newton::Task_PP_Newton(const Pipiline_parameters& pipiline_parameters_PP_Newton,
	const Oil_parameters& oil_parameters_PP_Newton) {
	setter(pipiline_parameters_PP_Newton, oil_parameters_PP_Newton);

};

// ������ ������
void Task_PP_Newton::setter(const Pipiline_parameters& pipiline_parameters_PP_Newton,
	const Oil_parameters& oil_parameters_PP_Newton) {
	m_pipiline_parameters_PP_Newton = pipiline_parameters_PP_Newton;
	m_oil_parameters_PP_Newton = oil_parameters_PP_Newton;
}



double Task_PP_Newton::solver_newton_rafson() {
	// ������� �������� �������� �� ���������
	fixed_solver_parameters_t<1, 0> parameters;
	// �������� ��������� ��� ������ ����������� �������
	fixed_solver_result_t<1> result;
	// ������� ������� ����������� ��������� <2> � ������� �������� ������� - �������
	// m_initial_speed_approximation - ��������� �����������
	m_initial_speed_approximation = 0;
	Task_PP_Newton task_PP_Newton(m_pipiline_parameters_PP_Newton, m_oil_parameters_PP_Newton);
	fixed_newton_raphson<1>::solve_dense(task_PP_Newton, { m_initial_speed_approximation }, parameters, &result);
	// ��������� ������ task_PP_Newton ������ Bernoulli_equation
	Bernoulli_equation Bernoulli_PP_Newton(m_pipiline_parameters_PP_Newton, m_oil_parameters_PP_Newton,
		0, result.argument, m_internal_diameter);
	// Q - �������� ������ [�^3/�]
	double Q = Bernoulli_PP_Newton.volume_flow();
	return Q;
}