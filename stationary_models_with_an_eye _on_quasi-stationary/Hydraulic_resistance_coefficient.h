#pragma once
#include <cmath>


/// @brief Hydraulic_resistance_coefficient - ����� ��� ����������� ��������������� ������������� 
/// � ����������� �� ����� ����������
class  Hydraulic_resistance_coefficient {

	// ���� ������
	// m_Re -  ����� ���������� (Re)
	double m_Re;
	// m_hydraulic_resistance - ����������� ��������������_������������� (lambda)
	double m_hydraulic_resistance;
	// m_relative_roughness - ������������� ������������� ������������� (e)
	double m_relative_roughness;

public:
	// ����������� ������ Hydraulic_resistance_coefficient �� ��������� (����������)
	Hydraulic_resistance_coefficient();

	// ����������� ������ Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient(double Re, double relative_roughness);

	/// @brief setter - ������ ������������ ������
	/// @param Re - ����� ����������
	/// @param relative_roughness - ������������� ������������� ������������� (e)
	void setter(double Re, double relative_roughness);

	/// @brief ������� ������
	/// @return m_hydraulic_resistance - ����������� ��������������� ������������� (lambda)
	double stokes_formula();

	/// @brief ������� ��������
	/// @return m_hydraulic_resistance - ����������� ��������������� ������������� (lambda)
	double blasius_formula();

	/// @brief ������� ��������
	/// @return m_hydraulic_resistance - ����������� ��������������� ������������� (lambda)
	double altschul_formula();

	/// @brief �������������� ������� ������������ ��������������� ������������� �� ����������� ����� ����������
	/// @return m_hydraulic_resistance - ����������� ��������������_������������� (lambda)
	double calculation_hydraulic_resistance_coefficient();



};
