#pragma once
#include <iomanip>
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>
#include <fixed/fixed_nonlinear_solver.h>

// ����������� ������ ��� ����������� ��������������� ������������� 
#include "Hydraulic_resistance_coefficient.h"
// ����������� ������ ��� ������� ����� �� ����� 2 - ���������� ������������ ������� 
// � �������� �� �������������� (��������� ��������)
#include "Bernoulli_equation.h"
#include "struct.h"
#include "const.h"
#include "Solver_Eyler.h"




/// @brief PP_solver_Newton - ����� ���������� �� �������� ������� �������-������� ��� ������ PP
// <1> - ����������� ������� ��������� - ��������� ������
class PP_solver_Newton_QP_Eyler: public fixed_system_t<1>
{

private:
	// ���������� ����� ������

	// pipeline_parameters_PP_Newton - ���� ������ (��������� � ������ Pipeline_parameters ��� ���������� pipeline_parameters_PP_Newton)
	Pipeline_parameters m_pipeline_parameters_PP_Newton_QP_Eyler;
	// m_oil_parameters_PP_Newton - ���� ������ (��������� � ������ Oil_parameters ��� ���������� oil_parameters_PP_Newton)
	Oil_parameters m_oil_parameters_PP_Newton_QP_Eyler;
	// m_hydraulic_resistance - ���� ������ - ����������� ��������������� �������������
	double m_hydraulic_resistance;
	// m_d - ���� ������ - ���������� ������� �����[�]
	double m_d;
	// m_initial_speed_approximation - ��������� ����������� �������� ������� �����, [�/�]
	double m_initial_flow_approximation;

	using fixed_system_t<1>::var_type;

public:
	/// @brief PP_solver_newton - ����������� ������ ��� ������ PP ������� �������-�������
	PP_solver_Newton_QP_Eyler(Pipeline_parameters pipeline_parameters_PP_Newton_QP_Eyler, 
		Oil_parameters oil_parameters_PP_Newton_QP_Eyler);
	/// @brief residuals - ������� �������
	/// @param - v - ������� �������� (��������, [�/�])
	var_type residuals(const var_type& m_initial_flow_approximation) {
		m_pipeline_parameters_PP_Newton_QP_Eyler.Q = m_initial_flow_approximation;
		double pressure_p0_calculated = solver_QP_Eyler2(m_pipeline_parameters_PP_Newton_QP_Eyler, m_oil_parameters_PP_Newton_QP_Eyler);
		
		double result;
		result = m_oil_parameters_PP_Newton_QP_Eyler.p0 - pressure_p0_calculated;
		return result;
	}

	double solver_newton_rafson();
};



