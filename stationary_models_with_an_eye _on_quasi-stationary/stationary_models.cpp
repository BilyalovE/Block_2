/*!
	\ brief Блок 2 - Реализация стационарных моделей с прицелом на квазистационар
	\ author Bilyalov Eldar
	\ version 4 - Решение задач, рефакторинг, многофайловый проект 
	\date 16.02.2024
*/
#pragma once
// Подключаем необходимые библиотеки
#include <iostream>
#include "gtest/gtest.h"


// Подключение класса для определения гидравлического сопротивления 
#include "Hydraulic_resistance_coefficient.h"
// Подключение класса для решения задач из блока 2 - Реализация стационарных моделей 
// с прицелом на квазистационар (Уравнение Бернулли)
#include "Bernoulli_equation.h"
// Подключение класса для решения задачи QP численным методом Эйлера
#include "Task_QP_Eyler.h"
#include "struct.h"
#include "const.h"
// Подключение класса для решения задачи PP методом Ньютона-Рафсона
#include "Task_PP_Newton.h"

// Используем пространство имен std и pde_solvers
using namespace std;
using namespace pde_solvers;



/// @brief solver_QP - функция решения задачи QP
/// @param pipeline_parameters_QP - структура параметров трубопровода для задачи QP
/// @param oil_parameters_QP - структура параметров нефти для задачи QP
/// @return pressure_p0 - давление в начале участка нефтепровода [Па]
double solver_QP(const Pipeline_parameters& pipeline_parameters_QP, const Oil_parameters& oil_parameters_QP) {
	// Объявляем объект task_QP класса Bernoulli_equation
	Bernoulli_equation task_QP(pipeline_parameters_QP, oil_parameters_QP);
	// internal_diameter - внутренний диаметр трубы [м]
	double internal_diameter = task_QP.internal_diameter();
	// relative_roughness - относительная эквивалентная шероховатость
	double relative_roughness = task_QP.relative_roughness();
	// v - скорость течения нефти [м/с]
	double v = task_QP.speed_flow();
	// Re - число Рейнольдса
	double Re = task_QP.reynolds_number();
	// pressure_p0 - давление в начале участка нефтепровода [Па]
	double pressure_p0;
	// hydraulic_resistance - коэффициент гидравлического сопротивления (lambda)
	double hydraulic_resistance;
	// Объявляем объект lambda_QP класса Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient lambda_QP(Re, relative_roughness);
	// Расчёт коэффициента гидравлического сопротивления (lambda)
	hydraulic_resistance = lambda_QP.calculation_hydraulic_resistance_coefficient();
	// Вызов сеттера для передачи дополнительных параметров (перегрузка конструктора) 
	task_QP.setter2(pipeline_parameters_QP, oil_parameters_QP, hydraulic_resistance, v, internal_diameter);
	//Расчёт pressure_p0 (давление в начале трубопровода)
	pressure_p0 = task_QP.pressure_p0();
	return pressure_p0;
}

/// @brief solver_QP_Eyler - функция решения задачи QP методом Эйлера
/// @param pipeline_parameters_QP_Eyler - структура параметров трубопровода
/// @param oil_parameters_QP_Eyler - структура параметров нефти
/// @return pressure_p0 - давление в начале участка нефтепровода [Па]
double solver_QP_Eyler(const Pipeline_parameters& pipeline_parameters_QP_Eyler,
	const Oil_parameters& oil_parameters_QP_Eyler)
{
	// Объявляем объект Bernoulli_equation_QP_Eyler класса Bernoulli_equation
	Bernoulli_equation Bernoulli_equation_QP_Eyler(pipeline_parameters_QP_Eyler, oil_parameters_QP_Eyler);
	// internal_diameter - внутренний диаметр трубы [м]
	double internal_diameter = Bernoulli_equation_QP_Eyler.internal_diameter();
	// v - скорость течения нефти [м/с]
	double v = Bernoulli_equation_QP_Eyler.speed_flow();
	// relative_roughness - относительная эквивалентная шероховатость
	double relative_roughness = Bernoulli_equation_QP_Eyler.relative_roughness();
	// Re - число Рейнольдса
	double Re = Bernoulli_equation_QP_Eyler.reynolds_number();
	// pressure_p0 - давление в начале участка нефтепровода[Па]
	double pressure_p0;
	// hydraulic_resistance - коэффициент гидравлического сопротивления (lambda)
	double hydraulic_resistance;
	// Объявляем объект lambda_QP_Eyler класса Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient lambda_QP_Eyler(Re, relative_roughness);
	// Расчет коэффициента гидравлического сопротивления (lambda)
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
	double h = pipeline_parameters_QP_Eyler.l / n;
	Task_QP_Eyler task_QP_Eyler(pipeline_parameters_QP_Eyler, oil_parameters_QP_Eyler,
		tw, n, h, internal_diameter, pressure_previous);
	// Вызов итеративного метода Эйлера
	pressure_p0 = task_QP_Eyler.solver_eyler();
	return pressure_p0;
}

/// @brief solver_PP - функция решения задачи PP
/// @param pipeline_parameters_PP - структура параметров трубопровода
/// @param oil_parameters_PP - структура параметров нефти
/// @return Q -  расход [м^3/с]
double solver_PP(const Pipeline_parameters& pipeline_parameters_PP, const Oil_parameters& oil_parameters_PP) {
	// lambda_previous - гидравлическое сопротивление на предыдущем приближении
	double lambda_previous = 0.021;
	// lambda_current - гидравлическое сопротивление на текущем приближении
	double lambda_current = 0.02;
	// eps - допустимая погрешность
	double eps = 0.0005;
	// Объявляем объект task_PP класса Bernoulli_equation
	Bernoulli_equation task_PP(pipeline_parameters_PP, oil_parameters_PP);
	// internal_diameter - внутренний диаметр трубы[м]
	double internal_diameter = task_PP.internal_diameter();
	double initial_v{ 0 };
	task_PP.setter2(pipeline_parameters_PP, oil_parameters_PP, lambda_current, initial_v, internal_diameter);
	// v - скорость течения нефти в системе СИ
	double v;
	double Q;
	// relative_roughness - относительная эквивалентная шероховатость
	double relative_roughness = task_PP.relative_roughness();
	// Объявляем объект lambda_PP класса Hydraulic_resistance_coefficient 
	Hydraulic_resistance_coefficient lambda_PP;
	while (abs(lambda_previous - lambda_current) >= eps) {
		v = task_PP.speed_pressure(lambda_current);
		// Re - число Рейнольдса
		double Re = task_PP.reynolds_number();
		lambda_previous = lambda_current;
		lambda_PP.setter(Re, relative_roughness);
		lambda_current = lambda_PP.calculation_hydraulic_resistance_coefficient();
		task_PP.setter2(pipeline_parameters_PP, oil_parameters_PP, lambda_current, v, internal_diameter);
	}
	Q = task_PP.volume_flow();
	return Q;
}


/// @brief solver_PP_Newton - функция решения задачи PP_Newton
/// @param Pipeline_parameters - структура параметров трубопровода
/// @param Oil_parameters - структура параметров нефти
/// @return Q - расход
double solver_PP_Newton(Pipeline_parameters pipeline_parameters_PP_Newton, Oil_parameters oil_parameters_PP_Newton) {
	// Создание экземпляра класса, который и будет решаемой системой
	PP_solver_Newton solver_newton(pipeline_parameters_PP_Newton, oil_parameters_PP_Newton);
	double Q = solver_newton.solver_newton_rafson();
	return Q;
}


TEST(Block_2, Task_QP) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipeline_parameters_QP
	Pipeline_parameters pipeline_parameters_QP = { 0.720, 0.010, 100, 50, 0.00015, 80000 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_QP
	Oil_parameters oil_parameters_QP = { 870, 15e-6, 0, 0.6e6, 0.9722 };
	double pressure_p0 = solver_QP(pipeline_parameters_QP, oil_parameters_QP);
	double abs_error = 0.01e6;
	EXPECT_NEAR(6.03e6, pressure_p0, abs_error);
}

TEST(Block_2, Task_QP_Iterative_Eyler) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipeline_parameters_QP_Eyler
	Pipeline_parameters pipeline_parameters_QP_Eyler = { 0.720, 0.010, 100, 50, 0.00015, 80000 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_QP_Eyler
	Oil_parameters oil_parameters_QP_Eyler = { 870, 15e-6, 0, 0.6e6, 0.9722 };
	double pressure_p0 = solver_QP_Eyler(pipeline_parameters_QP_Eyler, oil_parameters_QP_Eyler);
	double abs_error = 0.01e6;
	EXPECT_NEAR(6.03e6, pressure_p0, abs_error);
}

TEST(Block_2, Task_PP) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipeline_parameters_PP
	Pipeline_parameters pipeline_parameters_PP = { 0.720, 0.010, 50, 100, 0.00015, 80000 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_PP
	Oil_parameters oil_parameters_PP = { 870, 15e-6, 5e6, 0.8e6 };
	double Q = solver_PP(pipeline_parameters_PP, oil_parameters_PP);
	double abs_error = 4;
	cout <<"Задача PP = " << Q * 3600 << endl;
	EXPECT_NEAR(2739, Q * 3600, abs_error);
}

TEST(Block_2, Task_PP_Newton) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipeline_parameters_PP_Newton
	Pipeline_parameters pipeline_parameters_PP_Newton = { 0.720, 0.010, 50, 100, 0.00015, 80000 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_PP_Newton
	Oil_parameters oil_parameters_PP_Newton = { 870, 15e-6, 5e6, 0.8e6 };
  	double Q = solver_PP_Newton(pipeline_parameters_PP_Newton, oil_parameters_PP_Newton);
	cout << "Задача PP_Newton= " << Q * 3600 << endl;
	double abs_error = 9;
	EXPECT_NEAR(2739, Q * 3600, abs_error);
}

TEST(Block_2, Task_PP_Iterative_Eyler) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipeline_parameters_QP_Eyler
	Pipeline_parameters pipeline_parameters_PP_Eyler = { 0.720, 0.010, 100, 50, 0.00015, 80000 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_QP_Eyler
	Oil_parameters oil_parameters_PP_Eyler = { 870, 15e-6, 0, 0.6e6, 0.9722 };
	double pressure_p0 = solver_QP_Eyler(pipeline_parameters_PP_Eyler, oil_parameters_PP_Eyler);
	
	double abs_error = 4;
	//EXPECT_NEAR(2739, Q * 3600, abs_error);
}