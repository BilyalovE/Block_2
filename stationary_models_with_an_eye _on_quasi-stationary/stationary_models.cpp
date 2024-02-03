/*!
	\brief Блок 2 - Реализация стационарных моделей с прицелом на квазистационар
	\author Bilyalov Eldar
	\version 1.01 - класс Hydraulic_resistance_coefficient реализован в отдельных файлах
	\date 01.02.2024
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

#include "struct.h"
#include "const.h"

// Используем простравнство имен std и pde_solvers
using namespace std;
using namespace pde_solvers;

///// @brief k_g - ускорение свободного падаения
//const double k_g { 9.81 };
///// @brief k_pi - число пи
//const double k_pi { 3.14 };

///// @brief Pipiline_parameters - Структура парметров трубопровода
///// @param D - внешний диаметр [м]
///// @param ds - толщина стенки [м]
///// @param z0 - высота в начале участка трубопровода [м]
///// @param zl - высотка в конце участка трубопровода [м]
///// @param delta - абсолютная шероховатость трубы [м]
///// @param l - длина трубопровода [м]
//struct Pipiline_parameters {
//	// D - внешний диаметр [м]
//	double D; 
//	// ds - толщина стенки [м]
//	double ds; 
//	// z0 - высота в начале участка трубопровода [м]
//	double z0; 
//	// zl - высотка в конце участка трубопровода [м]
//	double zl; 
//	// delta - абсолютная шероховатость трубы [м]
//	double delta; 
//	// l - длина трубопровода [м]
//	double l; 
//};
//
///// @brief Oil_parameters - Структура парметров нефти
///// @param ro - плотность нефти [кг/м^3]
///// @param nu - кинематическая вязкость [Ст]
///// @param p0 - давление в начале участка нефтепровода [Па]
///// @param pl - давление в конце участка нефтепровода [Па]
///// @param Q -  расход [м^3/с]
//struct Oil_parameters {
//	// ro - плотность нефти [кг/м^3]
//	double ro; 
//	// nu - кинематическая вязкость [Ст]
//	double nu; 
//	// p0 - давление в начале участка нефтепровода [Па]
//	double p0; 
//	// pl - давление в конце участка нефтепровода [Па]
//	double pl;
//	// Q - расход [м^3/с]
//	double Q; 
//};
//



///// @brief Iterative_solutions - класс для решения задач из блока 2 численными методами
//class Iterative_solutions {
//public:
//	/// @brief solver_eyler - численное интегрирование дифференциального уравнения методом Эйлера
//	/// @param pipiline_parameters_XX_Eyler - cтруктура парметров трубопровода
//	/// @param Oil_parameters - cтруктура парметров нефти
//	/// @param tw - касательное напряжение трения, учитывающее трение жидкости при течении по трубе
//	/// @param n - кол-во точек расчётной сетки
//	/// @param h - шаг по по координате расчетной сетки [м] 
//	/// @param p_prev - давление на предыдущей итерации (граничное условие) [МПа]
//	/// @return p_current - давление на текущей итерации (рассчитанное значение) [МПа]
//	double solver_eyler(Pipiline_parameters pipiline_parameters_XX, Oil_parameters oil_parameters_XX, double tw, int n, double h, double d, double p_prev) {
//		for (size_t i = 1; i <= n; i++) {
//			double p_current = (p_prev * 1000000 - h * (-4 * tw / d - oil_parameters_XX.ro * 9.81 * (pipiline_parameters_XX.zl - pipiline_parameters_XX.z0) / (n - 1) / h))/1000000;
//			return p_current;
//		}
//	}
//};


/// @brief solver_QP - функция решения задачи QP
/// @param pipiline_parameters_QP - cтруктура парметров трубопровода для задачи QP
/// @param oil_parameters_QP - cтруктура парметров нефти для задачи QP
/// @return p0 - давление в начале участка нефтепровода [Па]
double solver_QP(Pipiline_parameters &pipiline_parameters_QP, Oil_parameters &oil_parameters_QP) {
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
	// p0 - давление в начале участка нефтепровода [Па]
	double p0;
	// hydraulic_resistance - коэффициент гидравлического сопротивления (lambda)
	double hydraulic_resistance;
	// Объявляем объект lambda_QP класса Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient lambda_QP(Re, relative_roughnesse);
	// Вызов сеттора для передачи доп парметров и рассчета pressure_p0
	
	// Формула Стокса
	if (Re < 2000) {
		hydraulic_resistance = lambda_QP.stokes_formula();
		task_QP.setter2(pipiline_parameters_QP, oil_parameters_QP, hydraulic_resistance, v, internal_diameter);
		p0 = task_QP.pressure_p0();
		return p0; 
	}
	// Формула Блазиуса
	else if (Re >= 2000 && Re <= 4000) {
		hydraulic_resistance = lambda_QP.blasius_formula();
		task_QP.setter2(pipiline_parameters_QP, oil_parameters_QP, hydraulic_resistance, v, internal_diameter);
		p0 = task_QP.pressure_p0();
		return p0;
	}
	// Формула Альтшуля
	else {
		hydraulic_resistance = lambda_QP.altschul_formula();
		task_QP.setter2(pipiline_parameters_QP, oil_parameters_QP, hydraulic_resistance, v, internal_diameter);
		p0 = task_QP.pressure_p0();
		return p0;
	}
}

///// @brief solver_QP_Eyler - функция решения задачи QP методом Эйлера
///// @param pipiline_parameters_QP_Eyler - cтруктура парметров трубопровода
///// @param oil_parameters_QP_Eyler - cтруктура парметров нефти
///// @return p0 - давление в начале участка нефтепровода [МПа]
//double solver_QP_Eyler(Pipiline_parameters pipiline_parameters_QP_Eyler, Oil_parameters oil_parameters_QP_Eyler, int n, double h) {
//	// p_prev - давление на предыдущей итерации (граничное условие) [МПа]
//	double p_prev = oil_parameters_QP_Eyler.pl;
//	// p_current - давление на текущей итерации(рассчитанное значение) [МПа]
//	double p_current ;
//	// Объявляем объект task_PP класса Bernoulli_equation
//	Bernoulli_equation task_QP_Eyler;
//	// Объявляем объект lambda_QP_Eyler класса Hydraulic_resistance_coefficient
//	Hydraulic_resistance_coefficient lambda_QP_Eyler;
//	// Объявляем объект iterativetask_QP_Eyler класса Iterative_solutions
//	Iterative_solutions iterativetask_QP_Eyler;
//	// d - внутренний диаметр трубы[м]
//	double d = task_QP_Eyler.diameter(pipiline_parameters_QP_Eyler);
//	// v - cкорость течения нефти в системе СИ
//	double v = task_QP_Eyler.speed_flow(oil_parameters_QP_Eyler, d);
//	// e - относительная эквивалентная шероховатость
//	double e = task_QP_Eyler.relative_roughness(pipiline_parameters_QP_Eyler, d);
//	// Re - число Рейнольдса
//	double Re = task_QP_Eyler.reynolds_number(oil_parameters_QP_Eyler, v, d);
//	// lambda - коэффициент гидравлического сопротивления
//	double lambda;
//	// p0 - давление в начале участка нефтепровода[МПа]
//	double p0;
//	// Формула Стокса
//	if (Re < 2000) {
//		lambda = lambda_QP_Eyler.stokes_formula(Re);
//	}
//	// Формула Блазиуса
//	else if (Re >= 2000 && Re <= 4000) {
//		lambda = lambda_QP_Eyler.blasius_formula(Re);
//	}
//	// Формула Альтшуля
//	else {
//		lambda = lambda_QP_Eyler.altschul_formula(Re, e);
//	}
//	// tw - касательное напряжение трения, учитывающее трение жидкости при течении по трубе
//	double tw = lambda / 8 * oil_parameters_QP_Eyler.ro * pow(v, 2);
//	for (size_t i = 1; i <= n; i++)
//	{
//		p_current = iterativetask_QP_Eyler.solver_eyler(pipiline_parameters_QP_Eyler, oil_parameters_QP_Eyler, tw, n, h, d, p_prev);
//		p_prev = p_current;
//	}
//	p0 = p_current;
//	return p0;
//}
//
///// @brief solver_PP - функция решения задачи PP
///// @param Pipiline_parameters - cтруктура парметров трубопровода
///// @param Oil_parameters - cтруктура парметров нефти
///// @return 
//double solver_PP(Pipiline_parameters pipiline_parameters_PP, Oil_parameters oil_parameters_PP) {
//	// labmda_previous - гидравлическое сопротивление на предыдущем приближении
//	double labmda_previous = 0.021;
//	// labmda_current - гидравлическое сопротивление на текущем приближении
//	double labmda_current = 0.02;
//	// eps - допустимая погрешность
//	double eps = 0.00005;
//	// Объявляем объект task_PP класса Bernoulli_equation
//	Bernoulli_equation task_PP;
//	// d - внутренний диаметр трубы[м]
//	double d = task_PP.diameter(pipiline_parameters_PP);
//	// Объявляем объект lambda_PP класса Hydraulic_resistance_coefficient
//	Hydraulic_resistance_coefficient lambda_PP;
//	// v - cкорость течения нефти в системе СИ
//	double v = task_PP.speed_pressure(pipiline_parameters_PP, oil_parameters_PP, labmda_current, d);
//	while (abs(labmda_previous - labmda_current) >= eps) {
//
//		// Re - число Рейнольдса
//		double Re = task_PP.reynolds_number(oil_parameters_PP, v, d);
//		// e - относительную эквивалентная шероховатость
//		double e = task_PP.relative_roughness(pipiline_parameters_PP, d);
//		// Формула Стокса
//		if (Re < 2000) {
//			labmda_previous = labmda_current;
//			labmda_current = lambda_PP.stokes_formula(Re);
//		}
//		// Формула Блазиуса
//		else if (Re >= 2000 && Re <= 4000) {
//			labmda_previous = labmda_current;
//			labmda_current = lambda_PP.blasius_formula(Re);
//		}
//		// Формула Альтшуля
//		else {
//			labmda_previous = labmda_current;
//			labmda_current = lambda_PP.altschul_formula(Re, e);
//		}
//		v = task_PP.speed_pressure(pipiline_parameters_PP, oil_parameters_PP, labmda_current, d);
//	}
//	double Q = task_PP.volume_flow(v, d);
//	return Q;
//}
//
//
///// @brief solver_newton - метод Ньютона-Рафсона для решения систем нелинейных уравнений фиксированный размерности
//
//
///// @brief PP_solver_newton - класс основанный на решетеле методом Ньютона-Рафсона для задачи PP
//// <1> - Размерность системы уравнений - cкалярный случай
//class PP_solver_newton : public fixed_system_t<1>
//{
//private:
//	// Объявление полей класса
//	
//	// pipiline_parameters_PP_Newton - поле класса (структура с именем Pipeline_parameters для переменной pipiline_parameters_PP_Newton)
//	Pipiline_parameters m_pipiline_parameters_PP_Newton;
//	// m_oil_parameters_PP_Newton - поле класса (структура с именем Oil_parameters для переменной oil_parameters_PP_Newton)
//	Oil_parameters m_oil_parameters_PP_Newton;
//	// m_lambda - поле класса - коэффициент гидравлического сопротивления
//	double m_lambda;
//	// m_d - поле класса - внутренний диаметр трубы[м]
//	double m_d;
//	// m_initial_speed_approximation - начальное приближение cкорости течения нефти, [м/с]
//	double m_initial_speed_approximation;
//
//	using fixed_system_t<1>::var_type;
//
//public:
//	/// @brief PP_solver_newton - конструктор класса для задачи PP методом Ньютона-Рафсона
//	PP_solver_newton(Pipiline_parameters pipiline_parameters_PP_Newton, Oil_parameters oil_parameters_PP_Newton) {
//		m_pipiline_parameters_PP_Newton = pipiline_parameters_PP_Newton;
//		m_oil_parameters_PP_Newton = oil_parameters_PP_Newton;
//	}
//
//	/// @brief residuals - функция невязок
//	/// @param - v - искомый параметр (скорость, [м/с])
//	var_type residuals(const var_type &v) {
//		
//		// Объявляем объект task_PP_Newton класса Bernoulli_equation
//		Bernoulli_equation task_PP_Newton;
//		// d - внутренний диаметр трубы[м]
//		double d = task_PP_Newton.diameter(m_pipiline_parameters_PP_Newton);// Re - число Рейнольдса
//		double Re = task_PP_Newton.reynolds_number(m_oil_parameters_PP_Newton, v, d);
//		// relative_equivalent_roughness - e - относительную эквивалентная шероховатость
//		double relative_equivalent_roughness = task_PP_Newton.relative_roughness(m_pipiline_parameters_PP_Newton, d);
//		// hydraulic_resistance - гидравлическое_сопротивление (lambda)
//		double hydraulic_resistance = hydraulic_resistance_isaev(Re, relative_equivalent_roughness);
//		m_lambda = hydraulic_resistance;
//		m_d = d;
//		
//		// result - функция невязок
//		double result;
//		result = v - task_PP_Newton.speed_pressure(m_pipiline_parameters_PP_Newton, m_oil_parameters_PP_Newton, m_lambda, m_d);
//		return result;
//	}
//
//	double solver_newton_rafson(){
//		// Задание настроек решателя по умолчанию
//		fixed_solver_parameters_t<1, 0> parameters;
//		// Создание структуры для записи результатов расчета
//		fixed_solver_result_t<1> result;
//		// Решение системы нелинейныйх уравнений <2> с помощью решателя Ньютона - Рафсона
//		// m_initial_speed_approximation - Начальное приближение
//		m_initial_speed_approximation = 0;
//		fixed_newton_raphson<1>::solve_dense(*this, { m_initial_speed_approximation }, parameters, &result);
//		// Объявляем объект task_PP_Newton класса Bernoulli_equation
//		Bernoulli_equation task_PP_Newton;
//		// Q - объемный расход[м ^ 3 / ч]
//		double Q = task_PP_Newton.volume_flow(result.argument, m_d);
//		return Q;
//	}
//
//};
//
///// @brief solver_PP_Newton - функция решения задачи PP_Newton
///// @param Pipiline_parameters - cтруктура парметров трубопровода
///// @param Oil_parameters - cтруктура парметров нефти
///// @return Q - расход
//double solver_PP_Newton(Pipiline_parameters pipiline_parameters_PP_Newton, Oil_parameters oil_parameters_PP_Newton) {
//	// Создание экземпляра класса, который и будет решаемой системой
//	PP_solver_newton solver_newton(pipiline_parameters_PP_Newton, oil_parameters_PP_Newton);
//	double Q = solver_newton.solver_newton_rafson();
//	return Q;
//}


TEST(Block_2, Task_QP) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipiline_parameters_QP
	Pipiline_parameters pipiline_parameters_QP = { 0.720, 0.010, 100, 50, 0.015, 80000 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_QP
	Oil_parameters oil_parameters_QP = { 870, 15e-6, 0, 0.6e6, 3500/3600};
	double p0 = solver_QP(pipiline_parameters_QP, oil_parameters_QP);

	EXPECT_EQ(6.03e6, p0);
}


//TEST(Block_2, Task_QP_Eyler) {
//	/// Объявление структуры с именем Pipeline_parameters для переменной pipiline_parameters_QP
//	Pipiline_parameters pipiline_parameters_QP_Eyler = { 720, 10, 100, 50, 0.15, 80 };
//	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_QP
//	Oil_parameters oil_parameters_QP_Eyler = { 870, 15, 0, 0.6, 3500 };
//	// n - кол-во точек расчётной сетки
//	int n = 1000;
//	// h - шаг по по координате расчетной сетки [м]
//	double h = pipiline_parameters_QP_Eyler.l * 1000 / n;
//	double pk = round(100*solver_QP_Eyler(pipiline_parameters_QP_Eyler, oil_parameters_QP_Eyler, n, h))/100;
//	EXPECT_EQ(6.03, pk);
//} 
//
//TEST(Block_2, Task_PP) {
//	/// Объявление структуры с именем Pipeline_parameters для переменной pipiline_parameters_PP
//	Pipiline_parameters pipiline_parameters_PP = { 720, 10, 50, 100, 0.15, 80 };
//	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_PP
//	Oil_parameters oil_parameters_PP = { 870, 15, 5, 0.8 };
//	double Q = solver_PP(pipiline_parameters_PP, oil_parameters_PP);
//	EXPECT_EQ(2739, Q);
//}
//
//TEST(Block_2, Task_PP_Newton) {
//	/// Объявление структуры с именем Pipeline_parameters для переменной pipiline_parameters_PP_Newton
//	Pipiline_parameters pipiline_parameters_PP_Newton = { 720, 10, 50, 100, 0.15, 80 };
//	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_PP_Newton
//	Oil_parameters oil_parameters_PP_Newton = { 870, 15, 5, 0.8 };
//	double Q = solver_PP_Newton(pipiline_parameters_PP_Newton, oil_parameters_PP_Newton);
//	EXPECT_EQ(2739, Q);
//}