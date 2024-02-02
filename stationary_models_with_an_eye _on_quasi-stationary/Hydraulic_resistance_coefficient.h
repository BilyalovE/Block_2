#pragma once



/// @brief Hydraulic_resistance_coefficient - класс для определения гидравлического сопротивления 
/// в зависимости от числа рейнольдса
class  Hydraulic_resistance_coefficient {
	
	// Поля класса
	// m_Re -  число рейнольдса
	double m_Re;
	// m_hydraulic_resistance - коэффициент гидравлическое_сопротивление (lambda)
	double m_hydraulic_resistance;
	// m_relative_roughness - относительная эквивалентная шероховатость
	double m_relative_roughness;

public:
	// Конструктор класса Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient(double Re, double relative_roughness);

	/// @brief Формула Стокса
	/// @return m_hydraulic_resistance - коэффициент гидравлического сопротивления (lambda)
	double stokes_formula() { return m_hydraulic_resistance; }

	/// @brief Формула Блазиуса
	/// @return m_hydraulic_resistance - коэффициент гидравлического сопротивления (lambda)
	double blasius_formula() { return m_hydraulic_resistance; }

	/// @brief Формула Альтшуля
	/// @return m_hydraulic_resistance - коэффициент гидравлического сопротивления (lambda)
	double altschul_formula() { return m_hydraulic_resistance; }
};
