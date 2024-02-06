
#include "Hydraulic_resistance_coefficient.h"
/// @brief Hydraulic_resistance_coefficient - ����������� ������ �� ��������� (����������)
Hydraulic_resistance_coefficient::Hydraulic_resistance_coefficient() {
}
/// @brief Hydraulic_resistance_coefficient - ����������� ������
/// @param Re - ����� ���������� (Re)
/// @param m_relative_roughness - ������������� ������������� ������������� (e)
Hydraulic_resistance_coefficient::Hydraulic_resistance_coefficient(double Re, double relative_roughness) {
	setter(Re, relative_roughness);
}

// ������ ������

void Hydraulic_resistance_coefficient::setter(double Re, double relative_roughness)
{
	m_Re = Re;
	m_relative_roughness = relative_roughness;
}

double Hydraulic_resistance_coefficient::stokes_formula()
{
	double m_hydraulic_resistance = 64 / m_Re;
	return m_hydraulic_resistance;
}

double Hydraulic_resistance_coefficient::blasius_formula()
{
	double m_hydraulic_resistance = 0.3164 / pow(m_Re, 0.25);
	return m_hydraulic_resistance;
}

double Hydraulic_resistance_coefficient::altschul_formula()
{
	double m_hydraulic_resistance = 0.11 * pow(m_relative_roughness + 68 / m_Re, 0.25);
	return m_hydraulic_resistance;
}

double Hydraulic_resistance_coefficient::calculation_hydraulic_resistance_coefficient()
{
	// ������� ������
	if (m_Re < 2000) {
		m_hydraulic_resistance = this->stokes_formula();
	}
	// ������� ��������
	else if (m_Re >= 2000 && m_Re <= 4000) {
		m_hydraulic_resistance = this->blasius_formula();
	}
	// ������� ��������
	else {
		m_hydraulic_resistance = this->altschul_formula();
	}
	return m_hydraulic_resistance;
}
