#include "Task_PP_Newton_QP_Eyler.h"


PP_solver_Newton_QP_Eyler::PP_solver_Newton_QP_Eyler(Pipeline_parameters pipeline_parameters_PP_Newton_QP_Eyler,
	Oil_parameters oil_parameters_PP_Newton_QP_Eyler)
{
	m_pipeline_parameters_PP_Newton_QP_Eyler = pipeline_parameters_PP_Newton_QP_Eyler;
	m_oil_parameters_PP_Newton_QP_Eyler = oil_parameters_PP_Newton_QP_Eyler;
}

double PP_solver_Newton_QP_Eyler::solver_newton_rafson()
{
	// ������� �������� �������� �� ���������
	fixed_solver_parameters_t<1, 0> parameters;
	// �������� ��������� ��� ������ ����������� �������
	fixed_solver_result_t<1> result;
	// ������� ������� ���������� ��������� <2> � ������� �������� ������� - �������
	// m_initial_speed_approximation - ��������� �����������
	m_initial_flow_approximation = 0.5;

	fixed_newton_raphson<1>::solve_dense(*this, { m_initial_flow_approximation }, parameters, &result);
	
	//cout << result.argument << endl;
	double Q = result.argument;
	return Q;

}

