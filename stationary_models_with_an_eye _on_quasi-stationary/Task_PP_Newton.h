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
	// Объявление полей класса
	// Поля класса
	// pipiline_parameters_PP_Newton - структура параметров трубопровода
	Pipiline_parameters m_pipiline_parameters_PP_Newton;
	// m_oil_parameters_PP_Newton - струтура парметров нефти
	Oil_parameters m_oil_parameters_PP_Newton;
	// m_hydraulic_resistance - коэффициент гидравлическое_сопротивление(lambda)
	double m_hydraulic_resistance;
	// m_internal_diameter - внутренний диаметр трубы, [м]
	double m_internal_diameter;
	// m_initial_speed_approximation - начальное приближение cкорости течения нефти, [м/с]
	double m_initial_speed_approximation;

	using fixed_system_t<1>::var_type;

public:
	/// @brief PP_solver_newton - конструктор класса для задачи PP методом Ньютона-Рафсона
	Task_PP_Newton(const Pipiline_parameters& pipiline_parameters_PP_Newton,
		const Oil_parameters& oil_parameters_PP_Newton);

	/// @brief Cетер класса
	/// @param pipiline_parameters_PP_Newton - структура параметров трубопровода
	/// @param oil_parameters_PP_Newton - струтура парметров нефти
	void setter(const Pipiline_parameters& pipiline_parameters_PP_Newton,
		const Oil_parameters& oil_parameters_PP_Newton);

	/// @brief residuals - функция невязок
	/// @param v - искомый параметр (скорость, [м/с])
	var_type residuals(const var_type& v) {

		// Объявляем объект Bernoulli_PP_Newton класса Bernoulli_equation
		Bernoulli_equation Bernoulli_PP_Newton(m_pipiline_parameters_PP_Newton, m_oil_parameters_PP_Newton);
		// m_internal_diameter - внутренний диаметр трубы [м]
		m_internal_diameter = Bernoulli_PP_Newton.internal_diameter();
		double hydraulic_resistance { 0 };
		Bernoulli_PP_Newton.setter2(m_pipiline_parameters_PP_Newton, m_oil_parameters_PP_Newton,
			hydraulic_resistance, v, m_internal_diameter);
		// Re - число Рейнольдса
		double Re = Bernoulli_PP_Newton.reynolds_number();
		// relative_equivalent_roughness - e - относительную эквивалентная шероховатость
		double relative_equivalent_roughness = Bernoulli_PP_Newton.relative_roughness();
		// hydraulic_resistance - гидравлическое_сопротивление (lambda)
		hydraulic_resistance = hydraulic_resistance_isaev(Re, relative_equivalent_roughness);
		m_hydraulic_resistance = hydraulic_resistance;

		// result - функция невязок
		double result;
		result = v - Bernoulli_PP_Newton.speed_pressure();
		return result;
	}

	/// @brief solver_newton - метод Ньютона-Рафсона для решения систем нелинейных уравнений фиксированный размерности
	/// @return Q - объемный расход[м^3/ч]
	double solver_newton_rafson();
};


	