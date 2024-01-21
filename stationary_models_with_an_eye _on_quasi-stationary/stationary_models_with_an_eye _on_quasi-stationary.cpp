/*!
	\brief Блок 2 - Реализация стационарных моделей с прицелом на квазистационар
	\author Bilyalov Eldar
	\version 1
	\date 21.01.2024
*/

// Подключаем необходимые библиотеки
#include <iostream>
#include <cmath>
#include <locale.h>
#include "gtest/gtest.h"
#include <iomanip>
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

// Используем простравнство имен std
using namespace std;

/// @brief Pipiline_parameters - Структура парметров трубопровода
/// @param  D - внешний диаметр [мм]
/// @param ds - толщина стенки [мм]
/// @param z0 - высота в начале участка трубопровода [м]
/// @param zl - высотка в конце участка трубопровода [м]
/// @param delta - абсолютная шероховатость трубы [мм]
/// @param l - длина трубопровода [км]
struct Pipiline_parameters {
	double D; 
	double ds; 
	double z0; 
	double zl; 
	double delta; 
	double l; 
};



/// @brief Oil_parameters - Структура парметров нефти
/// @param ro - плотность нефти [кг/м3]
/// @param nu - кинематическая вязкость [сСт]
/// @param pl - давление в конце участка нефтепровода [МПа]
/// @param p0 - давление в начале участка нефтепровода [МПа]
/// @param Q -  расход [м^3/ч]
struct Oil_parameters {
	double ro; 
	double nu; 
	double pl; 
	double p0;
	double Q; 
};

/// @brief pressure_p0 - Функция, рассчитывающая давление в начале участка нефтепровода [МПа]
/// @param pipiline_parameters_QP - структура парметров трубопровода
/// @param oil_parameters_QP - структура парметров нефти
/// @param lambda - коэффициент гидравлического сопротивления
/// @param v - скорость течения нефти в системе СИ
/// @param d - внутренний диаметр трубы в системе СИ
/// @return oil_parameters_QP.p0 - давление в начале участка нефтепровода [МПа]
double pressure_p0(Pipiline_parameters pipiline_parameters_QP, Oil_parameters oil_parameters_QP, double lambda, double v, double d) {
	oil_parameters_QP.p0 = (oil_parameters_QP.ro * 9.81) * (oil_parameters_QP.pl * 1000000 / (oil_parameters_QP.ro * 9.81) - pipiline_parameters_QP.z0 + pipiline_parameters_QP.zl + (((lambda * pipiline_parameters_QP.l * 1000) / d) * pow(v, 2)) / (2 * 9.81)) * 0.000001;
	return oil_parameters_QP.p0;
}

/// @brief Hydraulic_resistance_coefficient - класс гидравлического сопротивления
class  Hydraulic_resistance_coefficient {
public:

	/// @brief Формула Стокса
	/// @param lambda - Коэффициент гидравлического сопротивления
	/// @return lambda
	double stokes_formula(double Re) {
		double lambda = 64 / Re;
		return lambda;
	}

	/// @brief Формула Блазиуса
	/// @param lambda - Коэффициент гидравлического сопротивления
	/// @return lambda
	double blasius_formula(double Re) {
		double lambda = 0.3164 / pow(Re, 0.25);
		return lambda;
	}

	/// @brief Формула Альтшуля
	/// @param lambda - Коэффициент гидравлического сопротивления
	/// @return lambda
	double altschul_formula(double Re, double e) {
		double lambda = 0.11 * pow(e + 68 / Re, 0.25);
		return lambda;
	}
};

/// @brief Функция решения задачи QP
/// @param Pipiline_parameters - cтруктура парметров трубопровода
/// @param Oil_parameters - cтруктура парметров нефти
/// @return p0 - давление в начале участка нефтепровода
double solver(Pipiline_parameters pipiline_parameters_QP, Oil_parameters oil_parameters_QP) {
	// Внутренний диаметр трубы в системе СИ
	double d = (pipiline_parameters_QP.D - 2 * pipiline_parameters_QP.ds)/1000; 
	// Относительная эквивалентная шероховатость
	double e = (pipiline_parameters_QP.delta / 1000) / d;
	// Скорость течения нефти в системе СИ
	double v = (4 * oil_parameters_QP.Q / 3600) / (3.14 * pow(d, 2));
	// Число Рейнольдса, где nu переведено в систему СИ
	double Re = v * d / (oil_parameters_QP.nu * 0.000001);
	
	// Объявляем объект lambda_QP класса Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient lambda_QP;
	
	// Формула Стокса
	if (Re < 2000) {
		double lambda = lambda_QP.stokes_formula(Re);
		return pressure_p0(pipiline_parameters_QP,oil_parameters_QP, lambda, v, d);
	}
	// Формула Блазиуса
	else if (Re >= 2000 && Re <= 4000) {
		double lambda = lambda_QP.blasius_formula(Re);
		return pressure_p0(pipiline_parameters_QP, oil_parameters_QP, lambda, v, d);
	}
	// Формула Альтшуля
	else {
		double lambda = lambda_QP.altschul_formula(Re, e);
		return pressure_p0(pipiline_parameters_QP, oil_parameters_QP, lambda, v, d);
	}
}

TEST(Block_2,  QP) {
	
	Pipiline_parameters pipiline_parameters_QP = { 720, 10, 100, 50, 0.15, 80 };
	Oil_parameters oil_parameters_QP = { 870, 15, 0.6, 0, 3500 };
	double p0 = round(100*solver(pipiline_parameters_QP, oil_parameters_QP))/100;
	EXPECT_EQ(6.03, p0);
}
