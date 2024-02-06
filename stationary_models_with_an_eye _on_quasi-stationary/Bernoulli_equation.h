#pragma once
#include "struct.h"
#include "const.h"
#include <cmath>
//#include "stationary_models.cpp"

/// @brief Bernoulli_equation - ����� ��� ������� ����� �� ����� 2 - ���������� ������������ ������� 
/// � �������� �� �������������� (��������� ��������)
class Bernoulli_equation
{
	// ���� ������

	// m_pipiline_parameters - ��������� ��������� ������������
	Pipiline_parameters m_pipiline_parameters;
	// m_oil_parameters - ��������� ��������� �����
	Oil_parameters m_oil_parameters;
	// m_hydraulic_resistance - ����������� ��������������_������������� (lambda)
	double m_hydraulic_resistance;
	// m_v - c������� ������� ����� [�/�]
	double m_v;
	// m_d - ���������� ������� ����� [�]
	double m_d;
	// m_relative_roughness - o������������ ������������� �������������(e)
	double m_relative_roughness;
	// m_Re - ����� ����������
	double m_Re;
	// m_Q - �������� ������ [�^3/c]
	double m_Q;
	// m_p0 - �������� � ������ ������� ������������ [��]
	double m_p0;

public:
	/// @brief ����������� ������ �� ��������� Bernoulli_equation
	Bernoulli_equation(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters);

	/// @brief ����������� ������ Bernoulli_equation
	Bernoulli_equation(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters,
		double hydraulic_resistance, double v, double d);

	void setter1(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters);

	/// @brief setter - ������ ����������� Bernoulli_equation
	/// @param pipiline_parameters - �������� ��������� ������������
	/// @param oil_parameters - ��������� ��������� �����
	/// @param hydraulic_resistance - ����������� ��������������_������������� (lambda)
	/// @param v - c������� ������� ����� [�/�]
	/// @param d - ���������� ������� ����� [�]
	void setter2(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters,
		double& hydraulic_resistance, const double& v, double& d);



	/// @brief pressure_p0 - �����, �������������� �������� � ������ ������� ������������ [��]
	/// @return m_p0 - �������� � ������ ������� ������������ [��]
	double pressure_p0();


	double internal_diameter();

	/// @brief  relative_roughness - �����, �������������� ������������� ������������� �������������
	/// @return relative_roughness - o������������ ������������� ������������� (e)
	double relative_roughness();

	/// @brief reynolds_number - �����, �������������� ����� ����������, ��� nu ���������� � ������� ��
	/// @return m_Re - ����� ����������
	double reynolds_number();

	double reynolds_number(double m_d);

	/// @brief speed_flow - �����, �������������� �������� �� ��������� ������� �����
	/// @return m_v - c������� ������� ����� [�/�]
	double speed_flow();

	/// @brief speed_pressure - �����, �������������� �������� �� �������� � ������ PP
	/// @return v - c������� ������� ����� � ������� ��	
	double speed_pressure();

	/// @brief volume_flow - �����, �������������� �������� ������
	/// @return m_Q - �������� ������ [�^3/c]
	double volume_flow();
};