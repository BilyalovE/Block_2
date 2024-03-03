#include "Task_PP_Newton.h"

PP_solver_Newton::PP_solver_Newton(Pipeline_parameters pipeline_parameters_PP_Newton, Oil_parameters oil_parameters_PP_Newton)
{
	m_pipeline_parameters_PP_Newton = pipeline_parameters_PP_Newton;
	m_oil_parameters_PP_Newton = oil_parameters_PP_Newton;
}

double PP_solver_Newton::solver_newton_rafson()
{
	// ������� �������� �������� �� ���������
	fixed_solver_parameters_t<1, 0> parameters;
	// �������� ��������� ��� ������ ����������� �������
	fixed_solver_result_t<1> result;
	// ������� ������� ���������� ��������� <2> � ������� �������� ������� - �������
	// m_initial_speed_approximation - ��������� �����������
	m_initial_speed_approximation = 0.2;
	fixed_newton_raphson<1>::solve_dense(*this, { m_initial_speed_approximation }, parameters, &result);
	// ��������� ������ task_PP_Newton ������ Bernoulli_equation
	Bernoulli_equation task_PP_Newton(m_pipeline_parameters_PP_Newton, m_oil_parameters_PP_Newton);
	// Q - �������� ������[� ^ 3 / �]
	//cout << result.argument << endl;
	double Q = task_PP_Newton.volume_flow(result.argument);
	//cout << Q*3600 << endl;
	return Q;
	
}

