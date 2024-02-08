#pragma once
#include "struct.h"
#include "const.h"
#include <cmath>
//#include "stationary_models.cpp"

/// @brief Bernoulli_equation - ����� ��� �襭�� ����� �� ����� 2 - ��������� ��樮����� ������� 
/// � ��楫�� �� �������樮��� (�ࠢ����� ���㫫�)
class Bernoulli_equation
{
	// ���� �����

	// m_pipiline_parameters - ������� ��ଥ�஢ ��㡮�஢���
	Pipiline_parameters m_pipiline_parameters;
	// m_oil_parameters - ������� ��ଥ�஢ ����
	Oil_parameters m_oil_parameters;
	// m_hydraulic_resistance - �����樥�� ���ࠢ���᪮�_ᮯ�⨢����� (lambda)
	double m_hydraulic_resistance;
	// m_v - c������ �祭�� ���� [�/�]
	double m_v;
	// m_d - ����७��� ������� ���� [�]
	double m_d;
	// m_relative_roughness - o⭮�⥫쭠� ���������⭠� ��客�����(e)
	double m_relative_roughness;
	// m_Re - �᫮ ��������
	double m_Re;
	// m_Q - ��ꥬ�� ��室 [�^3/c]
	double m_Q;
	// m_p0 - �������� � ��砫� ���⪠ ���⥯஢��� [��]
	double m_p0;

public:
	/// @brief �������⪮� ����� �� 㬮�砭�� Bernoulli_equation
	Bernoulli_equation(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters);

	/// @brief �������⪮� ����� Bernoulli_equation
	Bernoulli_equation(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters,
		double hydraulic_resistance, double v, double d);

	void setter1(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters);

	/// @brief setter - ���� ��������� Bernoulli_equation
	/// @param pipiline_parameters - ������ ��ଥ�஢ ��㡮�஢���
	/// @param oil_parameters - ������� ��ଥ�஢ ����
	/// @param hydraulic_resistance - �����樥�� ���ࠢ���᪮�_ᮯ�⨢����� (lambda)
	/// @param v - c������ �祭�� ���� [�/�]
	/// @param d - ����७��� ������� ���� [�]
	void setter2(const Pipiline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters,
		double& hydraulic_resistance, const double& v, double& d);



	/// @brief pressure_p0 - ��⮤, �����뢠�騩 �������� � ��砫� ���⪠ ���⥯஢��� [��]
	/// @return m_p0 - �������� � ��砫� ���⪠ ���⥯஢��� [��]
	double pressure_p0();


	double internal_diameter();

	/// @brief  relative_roughness - ��⮤, �����뢠�騩 �⭮�⥫��� ���������⭠� ��客�����
	/// @return relative_roughness - o⭮�⥫쭠� ���������⭠� ��客����� (e)
	double relative_roughness();

	/// @brief reynolds_number - ��⮤, �����뢠�騩 �᫮ ��������, ��� nu ��ॢ����� � ��⥬� ��
	/// @return m_Re - �᫮ ��������
	double reynolds_number();

	double reynolds_number(double m_d);

	/// @brief speed_flow - ��⮤, �����뢠�騩 ᪮���� �� ��������� ��室� ����
	/// @return m_v - c������ �祭�� ���� [�/�]
	double speed_flow();

	/// @brief speed_pressure - ��⮤, �����뢠�騩 ᪮���� �� �������� � ����� PP
	/// @return v - c������ �祭�� ���� � ��⥬� ��	
	double speed_pressure();
	double speed_pressure(double m_hydraulic_resistance);

	/// @brief volume_flow - ��⮤, �����뢠�騩 ��ꥬ�� ��室
	/// @return m_Q - ��ꥬ�� ��室 [�^3/c]
	double volume_flow();
	double volume_flow(double m_d);
};