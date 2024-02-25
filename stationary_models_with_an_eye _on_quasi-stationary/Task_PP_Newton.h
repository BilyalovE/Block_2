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


/// @brief PP_solver_Newton - ����� ���������� �� �������� ������� �������-������� ��� ������ PP
// <1> - ����������� ������� ��������� - ��������� ������
class PP_solver_Newton : public fixed_system_t<1>
{

private:
	// ���������� ����� ������

	// pipeline_parameters_PP_Newton - ���� ������ (��������� � ������ Pipeline_parameters ��� ���������� pipeline_parameters_PP_Newton)
	Pipeline_parameters m_pipeline_parameters_PP_Newton;
	// m_oil_parameters_PP_Newton - ���� ������ (��������� � ������ Oil_parameters ��� ���������� oil_parameters_PP_Newton)
	Oil_parameters m_oil_parameters_PP_Newton;
	// m_hydraulic_resistance - ���� ������ - ����������� ��������������� �������������
	double m_hydraulic_resistance;
	// m_d - ���� ������ - ���������� ������� �����[�]
	double m_d;
	// m_initial_speed_approximation - ��������� ����������� �������� ������� �����, [�/�]
	double m_initial_speed_approximation;

	using fixed_system_t<1>::var_type;

public:
	/// @brief PP_solver_newton - ����������� ������ ��� ������ PP ������� �������-�������
	PP_solver_Newton(Pipeline_parameters pipeline_parameters_PP_Newton, Oil_parameters oil_parameters_PP_Newton);
	/// @brief residuals - ������� �������
	/// @param - v - ������� �������� (��������, [�/�])
	var_type residuals(const var_type& v) {
		// inner_diameter - ���������� ������� ����� [�]
		double inner_diameter = m_pipeline_parameters_PP_Newton.get_inner_diameter();
		// Re - ����� ����������
		Bernoulli_equation task_PP_Newton(m_pipeline_parameters_PP_Newton, m_oil_parameters_PP_Newton);
		double Re = task_PP_Newton.reynolds_number(v);
		// relative_equivalent_roughness - ������������� ������������� ������������� (e)
		double relative_equivalent_roughness = m_pipeline_parameters_PP_Newton.get_relative_roughness();
		// hydraulic_resistance - ��������������_������������� (lambda)
		Hydraulic_resistance_coefficient hydraulic_task_PP_Newton(Re, relative_equivalent_roughness);
		m_hydraulic_resistance = hydraulic_task_PP_Newton.calculation_hydraulic_resistance_coefficient();
		//m_hydraulic_resistance = hydraulic_resistance_isaev(Re, relative_equivalent_roughness);
		// result - ������� �������
		double result;
		task_PP_Newton.setter(m_pipeline_parameters_PP_Newton, m_oil_parameters_PP_Newton, m_hydraulic_resistance);
		result = v - task_PP_Newton.speed_pressure();
		return result;
	}

	double solver_newton_rafson();
};

	