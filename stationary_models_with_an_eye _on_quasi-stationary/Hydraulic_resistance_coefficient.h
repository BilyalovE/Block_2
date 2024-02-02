#pragma once



/// @brief Hydraulic_resistance_coefficient - ����� ��� ����������� ��������������� ������������� 
/// � ����������� �� ����� ����������
class  Hydraulic_resistance_coefficient {
	
	// ���� ������
	// m_Re -  ����� ����������
	double m_Re;
	// m_hydraulic_resistance - ����������� ��������������_������������� (lambda)
	double m_hydraulic_resistance;
	// m_relative_roughness - ������������� ������������� �������������
	double m_relative_roughness;

public:
	// ����������� ������ Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient(double Re, double relative_roughness);

	/// @brief ������� ������
	/// @return m_hydraulic_resistance - ����������� ��������������� ������������� (lambda)
	double stokes_formula() { return m_hydraulic_resistance; }

	/// @brief ������� ��������
	/// @return m_hydraulic_resistance - ����������� ��������������� ������������� (lambda)
	double blasius_formula() { return m_hydraulic_resistance; }

	/// @brief ������� ��������
	/// @return m_hydraulic_resistance - ����������� ��������������� ������������� (lambda)
	double altschul_formula() { return m_hydraulic_resistance; }
};
