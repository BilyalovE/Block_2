/*!
	\brief Блок 2 - Реализация стационарных моделей с прицелом на квазистационар
	\author Bilyalov Eldar
	\version 2 - реализован многофайловый проект
	\date 05.02.2024
*/
#pragma once
// Подключаем необходимые библиотеки
#include <iostream>
#include <cmath>
#include <locale.h>
#include <vector>
#include "gtest/gtest.h"
#include <iomanip>
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>
#include <fixed/fixed_nonlinear_solver.h>

// Подключение класса для определения гидравлического сопротивления 
#include "Hydraulic_resistance_coefficient.h"
// Подключение класса для решения задач из блока 2 - Реализация стационарных моделей 
// с прицелом на квазистационар (Уравнение Бернулли)
#include "Bernoulli_equation.h"
// Подключение класса для решения задачи QP численным методом Эйлера
#include "Task_QP_Eyler.h"
#include "struct.h"
#include "const.h"

//#include "Task_PP_Newton.h"

// Используем простравнство имен std и pde_solvers
using namespace std;
using namespace pde_solvers;



/// @brief solver_QP - функция решения задачи QP
/// @param pipiline_parameters_QP - cтруктура парметров трубопровода для задачи QP
/// @param oil_parameters_QP - cтруктура парметров нефти для задачи QP
/// @return pressure_p0 - давление в начале участка нефтепровода [Па]
double solver_QP(const Pipiline_parameters& pipiline_parameters_QP, const Oil_parameters& oil_parameters_QP) {
	// Объявляем объект task_QP класса Bernoulli_equation
	Bernoulli_equation task_QP(pipiline_parameters_QP, oil_parameters_QP);
	// internal_diameter - внутренний диаметр трубы [м]
	double internal_diameter = task_QP.internal_diameter();
	// relative_roughness - относительная эквивалентная шероховатость
	double relative_roughnesse = task_QP.relative_roughness();
	// v - cкорость течения нефти [м/с]
	double v = task_QP.speed_flow();
	// Re - число Рейнольдса
	double Re = task_QP.reynolds_number();
	// pressure_p0 - давление в начале участка нефтепровода [Па]
	double pressure_p0;
	// hydraulic_resistance - коэффициент гидравлического сопротивления (lambda)
	double hydraulic_resistance;
	// Объявляем объект lambda_QP класса Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient lambda_QP(Re, relative_roughnesse);
	// Рассчет коэффициента гидравлического сопротивления (lambda)
	hydraulic_resistance = lambda_QP.calculation_hydraulic_resistance_coefficient();
	// Вызов сеттора для передачи дополнительных парметров (перегрузка конструтора) 
	task_QP.setter2(pipiline_parameters_QP, oil_parameters_QP, hydraulic_resistance, v, internal_diameter);
	//Рассчет pressure_p0 (давление в начале трубопровода)
	pressure_p0 = task_QP.pressure_p0();
	return pressure_p0;
}

/// @brief solver_QP_Eyler - функция решения задачи QP методом Эйлера
/// @param pipiline_parameters_QP_Eyler - cтруктура парметров трубопровода
/// @param oil_parameters_QP_Eyler - cтруктура парметров нефти
/// @return pressure_p0 - давление в начале участка нефтепровода [Па]
double solver_QP_Eyler(const Pipiline_parameters& pipiline_parameters_QP_Eyler,
	const Oil_parameters& oil_parameters_QP_Eyler)
{
	// Объявляем объект Bernoulli_equation_QP_Eyler класса Bernoulli_equation
	Bernoulli_equation Bernoulli_equation_QP_Eyler(pipiline_parameters_QP_Eyler, oil_parameters_QP_Eyler);
	// internal_diameter - внутренний диаметр трубы [м]
	double internal_diameter = Bernoulli_equation_QP_Eyler.internal_diameter();
	// v - cкорость течения нефти [м/с]
	double v = Bernoulli_equation_QP_Eyler.speed_flow();
	// relative_roughness - относительная эквивалентная шероховатость
	double relative_roughnesse = Bernoulli_equation_QP_Eyler.relative_roughness();
	// Re - число Рейнольдса
	double Re = Bernoulli_equation_QP_Eyler.reynolds_number();
	// pressure_p0 - давление в начале участка нефтепровода[Па]
	double pressure_p0;
	// hydraulic_resistance - коэффициент гидравлического сопротивления (lambda)
	double hydraulic_resistance;
	// Объявляем объект lambda_QP_Eyler класса Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient lambda_QP_Eyler(Re, relative_roughnesse);
	// Рассчет коэффициента гидравлического сопротивления (lambda)
	hydraulic_resistance = lambda_QP_Eyler.calculation_hydraulic_resistance_coefficient();
	// tw - касательное напряжение трения, учитывающее трение жидкости при течении по трубе
	double tw = hydraulic_resistance / 8 * oil_parameters_QP_Eyler.ro * pow(v, 2);
	// pressure_previous - давление на предыдущей итерации (граничное условие) [Па]
	double pressure_previous = oil_parameters_QP_Eyler.pl;
	// pressure_current - давление на текущей итерации(рассчитанное значение) [Па]
	double pressure_current;
	// n - кол-во точек расчётной сетки
	int n = 1000;
	// h - шаг по по координате расчетной сетки [м]
	double h = pipiline_parameters_QP_Eyler.l / n;
	Task_QP_Eyler task_QP_Eyler(pipiline_parameters_QP_Eyler, oil_parameters_QP_Eyler,
		tw, n, h, internal_diameter, pressure_previous);
	// Вызов итеративного метода Эйлера
	pressure_p0 = task_QP_Eyler.solver_eyler();
	return pressure_p0;
}

/// @brief solver_PP - функция решения задачи PP
/// @param pipiline_parameters_PP - cтруктура парметров трубопровода
/// @param oil_parameters_PP - cтруктура парметров нефти
/// @return Q -  расход [м^3/с]
double solver_PP(const Pipiline_parameters& pipiline_parameters_PP, const Oil_parameters& oil_parameters_PP) {
	// labmda_previous - гидравлическое сопротивление на предыдущем приближении
	double labmda_previous = 0.021;
	// labmda_current - гидравлическое сопротивление на текущем приближении
	double labmda_current = 0.02;
	// eps - допустимая погрешность
	double eps = 0.0005;
	// Объявляем объект task_PP класса Bernoulli_equation
	Bernoulli_equation task_PP(pipiline_parameters_PP, oil_parameters_PP);
	// internal_diameter - внутренний диаметр трубы[м]
	double internal_diameter = task_PP.internal_diameter();
	double initial_v{ 0 };
	task_PP.setter2(pipiline_parameters_PP, oil_parameters_PP, labmda_current, initial_v, internal_diameter);
	// v - cкорость течения нефти в системе СИ
	double v = task_PP.speed_pressure();
	// relative_roughness - относительная эквивалентная шероховатость
	double relative_roughnesse = task_PP.relative_roughness();
	// Объявляем объект lambda_PP класса Hydraulic_resistance_coefficient 
	Hydraulic_resistance_coefficient lambda_PP;
	while (abs(labmda_previous - labmda_current) >= eps) {
		//v = task_PP.speed_pressure();
		// Re - число Рейнольдса
		double Re = task_PP.reynolds_number();
		labmda_previous = labmda_current;
		lambda_PP.setter(Re, relative_roughnesse);
		labmda_current = lambda_PP.calculation_hydraulic_resistance_coefficient();
		labmda_current = hydraulic_resistance_isaev(Re, relative_roughnesse);
		task_PP.setter2(pipiline_parameters_PP, oil_parameters_PP, labmda_current, initial_v, internal_diameter);
		//v = task_PP.speed_pressure();
	}
	double Q = task_PP.volume_flow();
	return Q;
}


/// @brief solver_newton - метод Ньютона-Рафсона для решения систем нелинейных уравнений фиксированный размерности


/// @brief PP_solver_newton - класс основанный на решетеле методом Ньютона-Рафсона для задачи PP
// <1> - Размерность системы уравнений - cкалярный случай
class PP_solver_newton : public fixed_system_t<1>
{
private:
	// Объявление полей класса

	// pipiline_parameters_PP_Newton - поле класса (структура с именем Pipeline_parameters для переменной pipiline_parameters_PP_Newton)
	Pipiline_parameters m_pipiline_parameters_PP_Newton;
	// m_oil_parameters_PP_Newton - поле класса (структура с именем Oil_parameters для переменной oil_parameters_PP_Newton)
	Oil_parameters m_oil_parameters_PP_Newton;
	// m_hydraulic_resistance - поле класса - коэффициент гидравлического сопротивления
	double m_hydraulic_resistance;
	// m_d - поле класса - внутренний диаметр трубы[м]
	double m_d;
	// m_initial_speed_approximation - начальное приближение cкорости течения нефти, [м/с]
	double m_initial_speed_approximation;

	using fixed_system_t<1>::var_type;

public:
	/// @brief PP_solver_newton - конструктор класса для задачи PP методом Ньютона-Рафсона
	PP_solver_newton(Pipiline_parameters pipiline_parameters_PP_Newton, Oil_parameters oil_parameters_PP_Newton) {
		m_pipiline_parameters_PP_Newton = pipiline_parameters_PP_Newton;
		m_oil_parameters_PP_Newton = oil_parameters_PP_Newton;
	}
	/// @brief residuals - функция невязок
	/// @param - v - искомый параметр (скорость, [м/с])
	var_type residuals(const var_type& v) {

		// Объявляем объект task_PP_Newton класса Bernoulli_equation
		Bernoulli_equation task_PP_Newton(m_pipiline_parameters_PP_Newton, m_oil_parameters_PP_Newton);
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

	double solver_newton_rafson() {
		// Задание настроек решателя по умолчанию
		fixed_solver_parameters_t<1, 0> parameters;
		// Создание структуры для записи результатов расчета
		fixed_solver_result_t<1> result;
		// Решение системы нелинейныйх уравнений <2> с помощью решателя Ньютона - Рафсона
		// m_initial_speed_approximation - Начальное приближение
		m_initial_speed_approximation = 0.4;
		fixed_newton_raphson<1>::solve_dense(*this, { m_initial_speed_approximation }, parameters, &result);
		// Объявляем объект task_PP_Newton класса Bernoulli_equation
		Bernoulli_equation task_PP_Newton(m_pipiline_parameters_PP_Newton, m_oil_parameters_PP_Newton, 0, result.argument, 0.7);
		// Q - объемный расход[м ^ 3 / ч]
		cout << result.argument << endl;
		double Q = task_PP_Newton.volume_flow();
		cout << Q*3600 << endl;
		return Q;
	}

};

/// @brief solver_PP_Newton - функция решения задачи PP_Newton
/// @param Pipiline_parameters - cтруктура парметров трубопровода
/// @param Oil_parameters - cтруктура парметров нефти
/// @return Q - расход
double solver_PP_Newton(Pipiline_parameters pipiline_parameters_PP_Newton, Oil_parameters oil_parameters_PP_Newton) {
	// Создание экземпляра класса, который и будет решаемой системой
	PP_solver_newton solver_newton(pipiline_parameters_PP_Newton, oil_parameters_PP_Newton);
	double Q = solver_newton.solver_newton_rafson();
	return Q;
}


TEST(Block_2, Task_QP) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipiline_parameters_QP
	Pipiline_parameters pipiline_parameters_QP = { 0.720, 0.010, 100, 50, 0.00015, 80000 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_QP
	Oil_parameters oil_parameters_QP = { 870, 15e-6, 0, 0.6e6, 0.9722 };
	double pressure_p0 = solver_QP(pipiline_parameters_QP, oil_parameters_QP);
	double abs_error = 0.01e6;
	EXPECT_NEAR(6.03e6, pressure_p0, abs_error);
}

TEST(Block_2, Task_QP_Iterative_Eyler) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipiline_parameters_QP_Eyler
	Pipiline_parameters pipiline_parameters_QP_Eyler = { 0.720, 0.010, 100, 50, 0.00015, 80000 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_QP_Eyler
	Oil_parameters oil_parameters_QP_Eyler = { 870, 15e-6, 0, 0.6e6, 0.9722 };
	double pressure_p0 = solver_QP_Eyler(pipiline_parameters_QP_Eyler, oil_parameters_QP_Eyler);
	double abs_error = 0.01e6;
	EXPECT_NEAR(6.03e6, pressure_p0, abs_error);
}

TEST(Block_2, Task_PP) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipiline_parameters_PP
	Pipiline_parameters pipiline_parameters_PP = { 0.720, 0.010, 50, 100, 0.000015, 80000 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_PP
	Oil_parameters oil_parameters_PP = { 870, 15e-6, 5e6, 0.8e6 };
	double Q = solver_PP(pipiline_parameters_PP, oil_parameters_PP);
	double abs_error = 7;
	cout << Q * 3600 << endl;
	EXPECT_NEAR(2739, Q * 3600, abs_error);
}

TEST(Block_2, Task_PP_Newton) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipiline_parameters_PP_Newton
	Pipiline_parameters pipiline_parameters_PP_Newton = { 0.720, 0.010, 50, 100, 0.00015, 80000 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_PP_Newton
	Oil_parameters oil_parameters_PP_Newton = { 870, 15e-6, 5e6, 0.8e6 };
  	double Q = solver_PP_Newton(pipiline_parameters_PP_Newton, oil_parameters_PP_Newton);
	double abs_error = 1;
	EXPECT_NEAR(2739, Q * 3600, abs_error);
}