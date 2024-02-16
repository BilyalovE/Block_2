#pragma once
#include <iomanip>
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>
#include <fixed/fixed_nonlinear_solver.h>

// Подключение класса для определения гидравлического сопротивления 
#include "Hydraulic_resistance_coefficient.h"
// Подключение класса для решения задач из блока 2 - Реализация стационарных моделей 
// с прицелом на квазистационар (Уравнение Бернулли)
#include "Bernoulli_equation.h"
#include "struct.h"
#include "const.h"


/// @brief PP_solver_Newton - класс основанный на решателе методом Ньютона-Рафсона для задачи PP
// <1> - Размерность системы уравнений - скалярный случай
class PP_solver_Newton : public fixed_system_t<1>
{

private:
	// Объявление полей класса

	// pipeline_parameters_PP_Newton - поле класса (структура с именем Pipeline_parameters для переменной pipeline_parameters_PP_Newton)
	Pipeline_parameters m_pipeline_parameters_PP_Newton;
	// m_oil_parameters_PP_Newton - поле класса (структура с именем Oil_parameters для переменной oil_parameters_PP_Newton)
	Oil_parameters m_oil_parameters_PP_Newton;
	// m_hydraulic_resistance - поле класса - коэффициент гидравлического сопротивления
	double m_hydraulic_resistance;
	// m_d - поле класса - внутренний диаметр трубы[м]
	double m_d;
	// m_initial_speed_approximation - начальное приближение скорости течения нефти, [м/с]
	double m_initial_speed_approximation;

	using fixed_system_t<1>::var_type;

public:
	/// @brief PP_solver_newton - конструктор класса для задачи PP методом Ньютона-Рафсона
	PP_solver_Newton(Pipeline_parameters pipeline_parameters_PP_Newton, Oil_parameters oil_parameters_PP_Newton);
	/// @brief residuals - функция невязок
	/// @param - v - искомый параметр (скорость, [м/с])
	var_type residuals(const var_type& v) {
		// Объявляем объект task_PP_Newton класса Bernoulli_equation
		Bernoulli_equation task_PP_Newton(m_pipeline_parameters_PP_Newton, m_oil_parameters_PP_Newton);
		// internal_diameter - внутренний диаметр трубы [м]
		double internal_diameter = task_PP_Newton.internal_diameter();// Re - число Рейнольдса
		double Re = task_PP_Newton.reynolds_number(v);
		// relative_equivalent_roughness - относительную эквивалентная шероховатость (e)
		double relative_equivalent_roughness = task_PP_Newton.relative_roughness();
		// hydraulic_resistance - гидравлическое_сопротивление (lambda)
		Hydraulic_resistance_coefficient hydraulic_task_PP_Newton(Re, relative_equivalent_roughness);
		m_hydraulic_resistance = hydraulic_task_PP_Newton.calculation_hydraulic_resistance_coefficient();
		//m_hydraulic_resistance = hydraulic_resistance_isaev(Re, relative_equivalent_roughness);
		// result - функция невязок
		double result;
		result = v - task_PP_Newton.speed_pressure(m_hydraulic_resistance);
		return result;
	}

	double solver_newton_rafson();
};

	