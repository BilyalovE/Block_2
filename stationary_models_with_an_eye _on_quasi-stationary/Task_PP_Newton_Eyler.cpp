#include "Task_PP_Newton_Eyler.h"

Task_PP_Newton_Eyler::Task_PP_Newton_Eyler(Pipeline_parameters pipeline_parameters, Oil_parameters oil_parameters, double v)
	: Task_QP_Eyler(pipeline_parameters, oil_parameters, v)
{

}

double Task_PP_Newton_Eyler::solver_newton_rafson(double initial_pressure_p0)
{
	// ������� �������� �������� �� ���������
	fixed_solver_parameters_t<1, 0> parameters;
	// �������� ��������� ��� ������ ����������� �������
	fixed_solver_result_t<1> result;
	// ������� ������� ���������� ��������� <2> � ������� �������� ������� - �������
	// m_initial_speed_approximation - ��������� �����������
	fixed_newton_raphson<1>::solve_dense(*this, { initial_pressure_p0 }, parameters, &result);
	// ��������� ������ task_PP_Newton ������ Bernoulli_equation
	Bernoulli_equation task_PP_Newton(pipeline_parameters, oil_parameters);
	// Q - �������� ������[� ^ 3 / �]
	//cout << result.argument << endl;
	double Q = task_PP_Newton.volume_flow(result.argument);
	//cout << Q*3600 << endl;
	return Q;
}