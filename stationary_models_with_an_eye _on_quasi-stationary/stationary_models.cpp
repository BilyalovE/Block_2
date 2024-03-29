﻿/*!
	\ brief Блок 2 - Реализация стационарных моделей с прицелом на квазистационар
	\ author Bilyalov Eldar
	\ version 5 - Решение последних двух задач, рефакторинг
	\date 09.03.2024
*/
#pragma once
// Подключаем необходимые библиотеки
#include <iostream>
#include "gtest/gtest.h"
#include <iomanip>

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
#include "Task_PP_Newton_Eyler.h"
 
//#include "Task_PP_Newton_QP_Eyler.h"

// Используем пространство имен std и pde_solvers
using namespace std;
using namespace pde_solvers;



/// @brief solver_QP - функция решения задачи QP
/// @param pipeline_parameters_QP - структура параметров трубопровода для задачи QP
/// @param oil_parameters_QP - структура параметров нефти для задачи QP
/// @return pressure_p0 - давление в начале участка нефтепровода [Па]
double solver_QP(const Pipeline_parameters& pipeline_parameters_QP, const Oil_parameters& oil_parameters_QP) {
	// relative_roughness - относительная эквивалентная шероховатость
	double relative_roughness = pipeline_parameters_QP.get_relative_roughness();
	// v - скорость течения нефти [м/с]
	double v = pipeline_parameters_QP.speed_flow();
	// Re - число Рейнольдса
	double Re = pipeline_parameters_QP.reynolds_number();
	// pressure_p0 - давление в начале участка нефтепровода [Па]
	double pressure_p0;
	// hydraulic_resistance - коэффициент гидравлического сопротивления (lambda)
	double hydraulic_resistance;
	// Объявляем объект lambda_QP класса Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient lambda_QP(Re, relative_roughness);
	// Расчёт коэффициента гидравлического сопротивления (lambda)
	hydraulic_resistance = lambda_QP.calculation_hydraulic_resistance_coefficient();
	// Вызов сеттера для передачи дополнительных параметров (перегрузка конструктора) 
	// Объявляем объект task_QP класса Bernoulli_equation
	Bernoulli_equation task_QP(pipeline_parameters_QP, oil_parameters_QP, hydraulic_resistance);
	//Расчёт pressure_p0 (давление в начале трубопровода)
	pressure_p0 = task_QP.pressure_p0();
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
	// internal_diameter - внутренний диаметр трубы[м]
	double inner_diameter = pipeline_parameters_PP.get_inner_diameter();
	// v - скорость течения нефти в системе СИ
	double v{};
	double Q;
	// relative_roughness - относительная эквивалентная шероховатость
	double relative_roughness = pipeline_parameters_PP.get_relative_roughness();
	// Объявляем объект lambda_PP класса Hydraulic_resistance_coefficient 
	Hydraulic_resistance_coefficient lambda_PP;
	Bernoulli_equation task_PP(pipeline_parameters_PP, oil_parameters_PP, lambda_current);
	while (abs(lambda_previous - lambda_current) >= eps) {
		v = task_PP.speed_pressure(lambda_current);
		// Re - число Рейнольдса
		double Re = task_PP.reynolds_number(v);
		lambda_previous = lambda_current;
		lambda_PP.setter(Re, relative_roughness);
		lambda_current = lambda_PP.calculation_hydraulic_resistance_coefficient();
		task_PP.setter(pipeline_parameters_PP, oil_parameters_PP, lambda_current);
	}
	Q = task_PP.volume_flow(v);
	return Q;
}


TEST(Block_2, Task_QP) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipeline_parameters_QP
	Pipeline_parameters pipeline_parameters_QP = { 0.720, 0.010, 100, 50, 0.00015, 80000, 0.9722, 15e-6 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_QP
	Oil_parameters oil_parameters_QP = { 870, 0, 0.6e6};
	double pressure_p0 = solver_QP(pipeline_parameters_QP, oil_parameters_QP);
	double abs_error = 0.01e6;
	EXPECT_NEAR(6.03e6, pressure_p0, abs_error);
}

TEST(Block_2, Task_QP_Eyler) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipeline_parameters_QP_Eyler
	Pipeline_parameters pipeline_parameters_QP_Eyler = { 0.720, 0.010, 100, 50, 0.00015, 80000, 0.9722, 15e-6 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_QP_Eyler
	Oil_parameters oil_parameters_QP_Eyler = { 870, 0, 0.6e6 };
	// v - скорость течения нефти [м/с]
	double v = pipeline_parameters_QP_Eyler.speed_flow();
	Task_QP_Eyler task_QP_Eyler(pipeline_parameters_QP_Eyler, oil_parameters_QP_Eyler, v);
	// Вызов итеративного метода Эйлера
	double pressure_p0 = task_QP_Eyler.solver_eyler(v);
	double abs_error = 0.01e6;
	EXPECT_NEAR(6.03e6, pressure_p0, abs_error);
}

TEST(Block_2, Task_PP) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipeline_parameters_PP
	Pipeline_parameters pipeline_parameters_PP = { 0.720, 0.010, 50, 100, 0.00015, 80000, 0, 15e-6 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_PP
	Oil_parameters oil_parameters_PP = { 870, 5e6, 0.8e6 };
	double Q = solver_PP(pipeline_parameters_PP, oil_parameters_PP);
	double abs_error = 4;
	//cout <<"Задача PP = " << Q * 3600 << endl;
	EXPECT_NEAR(2739, Q * 3600, abs_error);
}

TEST(Block_2, Task_PP_Newton) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipeline_parameters_PP_Newton
	Pipeline_parameters pipeline_parameters_PP_Newton = { 0.720, 0.010, 50, 100, 0.00015, 80000, 0, 15e-6 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_PP_Newton
	Oil_parameters oil_parameters_PP_Newton = { 870, 5e6, 0.8e6 };
	// Создание экземпляра класса, который и будет решаемой системой
	PP_solver_Newton solver_newton(pipeline_parameters_PP_Newton, oil_parameters_PP_Newton);
	double Q = solver_newton.solver_newton_rafson();
	//cout << "Задача PP_Newton= " << Q * 3600 << endl;
	double abs_error = 9;
	EXPECT_NEAR(2739, Q * 3600, abs_error);
}

// Задача PP поверх Эйлера
TEST(Block_2, Task_PP_Eyler) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipeline_parameters_PP_Iterative_Eyler
	Pipeline_parameters pipeline_parameters_PP_Eyler = { 0.720, 0.010, 100, 50, 0.00015, 80000, 0.9722, 15e-6 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_PP_Iterative_Eyler
	Oil_parameters oil_parameters_PP_Eyler = { 870, 0, 0.6e6 };
	// v - скорость течения нефти [м/с]
	double v = pipeline_parameters_PP_Eyler.speed_flow();
	Task_QP_Eyler task_PP_Eyler(pipeline_parameters_PP_Eyler, oil_parameters_PP_Eyler, v);
	// Вызов итеративного метода Эйлера
	double pressure_p0 = task_PP_Eyler.solver_eyler(v);
	oil_parameters_PP_Eyler.p0 = pressure_p0;
	double abs_error = 15;
	double Q = solver_PP(pipeline_parameters_PP_Eyler, oil_parameters_PP_Eyler);
	EXPECT_NEAR(3500, Q * 3600, abs_error);
}

// Задача PP на Ньютоне поверх Эйлера
TEST(Block_2, Task_PP_Newton_Eyler) {
	Pipeline_parameters pipeline_parameters_PP_Newton_Eyler = { 0.720, 0.010, 50, 100, 0.00015, 80000, 0, 15e-6 };
	Oil_parameters oil_parameters_PP_Newton_Eyler = { 870, 5e6, 0.8e6 };
	double initial_v = 0.5;
	Task_PP_Newton_Eyler solver_PP_Newton_Eyler(pipeline_parameters_PP_Newton_Eyler, oil_parameters_PP_Newton_Eyler, initial_v);
	double Q = solver_PP_Newton_Eyler.solver_newton_rafson();
	double abs_error = 8;
	EXPECT_NEAR(2739, Q * 3600, abs_error);
}
 