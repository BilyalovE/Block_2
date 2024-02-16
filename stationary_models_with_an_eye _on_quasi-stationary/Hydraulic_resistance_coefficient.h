#pragma once
#include <cmath>


/// @brief Hydraulic_resistance_coefficient - класс для определения гидравлического сопротивления 
/// в зависимости от числа Рейнольдса
class  Hydraulic_resistance_coefficient {

	// Поля класса
	// m_Re -  число Рейнольдса (Re)
	double m_Re;
	// m_hydraulic_resistance - коэффициент гидравлическое_сопротивление (lambda)
	double m_hydraulic_resistance;
	// m_relative_roughness - относительная эквивалентная шероховатость (e)
	double m_relative_roughness;

public:
	// Конструктор класса Hydraulic_resistance_coefficient по умолчанию (перегрузка)
	Hydraulic_resistance_coefficient();

	// Конструктор класса Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient(double Re, double relative_roughness);

	/// @brief setter - сеттер конструктора класса
	/// @param Re - число Рейнольдса
	/// @param relative_roughness - относительная эквивалентная шероховатость (e)
	void setter(double Re, double relative_roughness);

	/// @brief Формула Стокса
	/// @return m_hydraulic_resistance - коэффициент гидравлического сопротивления (lambda)
	double stokes_formula();

	/// @brief Формула Блазиуса
	/// @return m_hydraulic_resistance - коэффициент гидравлического сопротивления (lambda)
	double blasius_formula();

	/// @brief Формула Альтшуля
	/// @return m_hydraulic_resistance - коэффициент гидравлического сопротивления (lambda)
	double altschul_formula();

	/// @brief Автоматический расчет коэффициента гидравлического сопротивления по переданному числу Рейнольдса
	/// @return m_hydraulic_resistance - коэффициент гидравлическое_сопротивление (lambda)
	double calculation_hydraulic_resistance_coefficient();



};
