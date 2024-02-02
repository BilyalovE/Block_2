#pragma once

#include "stationary_models.cpp"

/// @brief Bernoulli_equation - ����� ��� ������� ����� �� ����� 2 - ���������� ������������ ������� 
/// � �������� �� �������������� (��������� ��������)
class Bernoulli_equation {
	
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
	/// @brief ����������� ������ Bernoulli_equation
	Bernoulli_equation(const Pipiline_parameters &pipiline_parameters, const Oil_parameters &oil_parameters, 
		const double &hydraulic_resistance, const double &v, const double &d); 

	/// @brief setter - ������ ����������� Bernoulli_equation
	/// @param pipiline_parameters - �������� ��������� ������������
	/// @param oil_parameters - ��������� ��������� �����
	/// @param hydraulic_resistance - ����������� ��������������_������������� (lambda)
	/// @param v - c������� ������� ����� [�/�]
	/// @param d - ���������� ������� ����� [�]
	void setter(Pipiline_parameters pipiline_parameters, Oil_parameters oil_parameters,
		double hydraulic_resistance, double v, double d);
	
	/// @brief pressure_p0 - �����, �������������� �������� � ������ ������� ������������ [��]
	/// @return m_p0 - �������� � ������ ������� ������������ [��]
	double pressure_p0() { return m_p0; }

	/// @brief internal_diameter - �����, �������������� ���������� ������� �����
	/// @return d - ���������� ������� ����� [�]
	double internal_diameter() { return m_d;}

	/// @brief  relative_roughness - �����, �������������� ������������� ������������� �������������
	/// @return relative_roughness - o������������ ������������� ������������� (e)
	double relative_roughness() { return m_relative_roughness; } 

	/// @brief reynolds_number - �����, �������������� ����� ����������, ��� nu ���������� � ������� ��
	/// @return m_Re - ����� ����������
	double reynolds_number() { return m_Re; }

	/// @brief speed_flow - �����, �������������� �������� �� ��������� ������� �����
	/// @return m_v - c������� ������� ����� [�/�]
	double speed_flow() { return m_v; }

	/// @brief speed_pressure - �����, �������������� �������� �� �������� � ������ PP
	/// @return v - c������� ������� ����� � ������� ��	
	double speed_pressure() { return m_v; }

	/// @brief volume_flow - �����, �������������� �������� ������
	/// @return m_Q - �������� ������ [�^3/c]
	double volume_flow() { return m_Q; }
};