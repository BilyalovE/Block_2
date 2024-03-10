#include "Task_PP_Newton_Eyler.h"

Task_PP_Newton_Eyler::Task_PP_Newton_Eyler(Pipeline_parameters pipeline_parameters, Oil_parameters oil_parameters, double v)
	: Task_QP_Eyler(pipeline_parameters, oil_parameters, v),
	m_pipeline_parameters(pipeline_parameters), m_oil_parameters(oil_parameters), v(v)
{

}

double Task_PP_Newton_Eyler::solver_newton_rafson()
{
	// ������� �������� �������� �� ���������
	fixed_solver_parameters_t<1, 0> parameters;
	// �������� ��������� ��� ������ ����������� �������
	fixed_solver_result_t<1> result;
	// ������� ������� ���������� ��������� <2> � ������� �������� ������� - �������
	// m_initial_speed_approximation - ��������� �����������
	fixed_newton_raphson<1>::solve_dense(*this, { v }, parameters, &result);
	// ��������� ������ task_PP_Newton ������ Bernoulli_equation
	//m_oil_parameters.p0 = result.argument;
	Bernoulli_equation task_PP_Newton_Eyler(m_pipeline_parameters, m_oil_parameters);
	// Q - �������� ������[� ^ 3 / �]
	//cout << result.argument << endl;
	//double Q = task_PP_Newton_Eyler.speed_pressure();
	//cout << Q*3600 << endl;
	return result.argument;
}