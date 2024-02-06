#pragma once



#include <iostream>
#include <cmath>
#include <locale.h>
#include <vector>
#include "gtest/gtest.h"
#include <iomanip>
#include <pde_solvers/pde_solvers.h>
#include <fixed/fixed.h>
#include <fixed/fixed_nonlinear_solver.h>
class Task_PP_Newton : public fixed_system_t<1>
{
private:
	// ���������� ����� ������
	// ���� ������
	// pipiline_parameters_PP_Newton - ��������� ���������� ������������
	Pipiline_parameters m_pipiline_parameters_PP_Newton;
	// m_oil_parameters_PP_Newton - �������� ��������� �����
	Oil_parameters m_oil_parameters_PP_Newton;
	// m_hydraulic_resistance - ����������� ��������������_�������������(lambda)
	double m_hydraulic_resistance;
	// m_internal_diameter - ���������� ������� �����, [�]
	double m_internal_diameter;
	// m_initial_speed_approximation - ��������� ����������� c������� ������� �����, [�/�]
	double m_initial_speed_approximation;

	using fixed_system_t<1>::var_type;

public:
	/// @brief PP_solver_newton - ����������� ������ ��� ������ PP ������� �������-�������
	Task_PP_Newton(const Pipiline_parameters& pipiline_parameters_PP_Newton,
		const Oil_parameters& oil_parameters_PP_Newton);

	/// @brief C���� ������
	/// @param pipiline_parameters_PP_Newton - ��������� ���������� ������������
	/// @param oil_parameters_PP_Newton - �������� ��������� �����
	void setter(const Pipiline_parameters& pipiline_parameters_PP_Newton,
		const Oil_parameters& oil_parameters_PP_Newton);

	/// @brief residuals - ������� �������
	/// @param v - ������� �������� (��������, [�/�])
	var_type residuals(const var_type& v) {

		// ��������� ������ Bernoulli_PP_Newton ������ Bernoulli_equation
		Bernoulli_equation Bernoulli_PP_Newton(m_pipiline_parameters_PP_Newton, m_oil_parameters_PP_Newton);
		// m_internal_diameter - ���������� ������� ����� [�]
		m_internal_diameter = Bernoulli_PP_Newton.internal_diameter();
		double hydraulic_resistance { 0 };
		Bernoulli_PP_Newton.setter2(m_pipiline_parameters_PP_Newton, m_oil_parameters_PP_Newton,
			hydraulic_resistance, v, m_internal_diameter);
		// Re - ����� ����������
		double Re = Bernoulli_PP_Newton.reynolds_number();
		// relative_equivalent_roughness - e - ������������� ������������� �������������
		double relative_equivalent_roughness = Bernoulli_PP_Newton.relative_roughness();
		// hydraulic_resistance - ��������������_������������� (lambda)
		hydraulic_resistance = hydraulic_resistance_isaev(Re, relative_equivalent_roughness);
		m_hydraulic_resistance = hydraulic_resistance;

		// result - ������� �������
		double result;
		result = v - Bernoulli_PP_Newton.speed_pressure();
		return result;
	}

	/// @brief solver_newton - ����� �������-������� ��� ������� ������ ���������� ��������� ������������� �����������
	/// @return Q - �������� ������[�^3/�]
	double solver_newton_rafson();
};


	