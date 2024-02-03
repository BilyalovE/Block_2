
#include "Hydraulic_resistance_coefficient.h"

/// @brief Hydraulic_resistance_coefficient - конструктор класса
/// @param Re - число Рейнольдса (Re)
/// @param m_relative_roughness - относительная эквивалентная шероховатость (e)
Hydraulic_resistance_coefficient::Hydraulic_resistance_coefficient(double Re, double relative_roughness){
	void setter(double Re, double relative_roughness);
	}

// Методы класса

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