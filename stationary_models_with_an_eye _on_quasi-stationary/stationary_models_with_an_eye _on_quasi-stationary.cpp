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
#include <vector>
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
/// @param p0 - давление в начале участка нефтепровода [МПа]
/// @param pl - давление в конце участка нефтепровода [МПа]
/// @param Q -  расход [м^3/ч]
struct Oil_parameters {
	double ro; 
	double nu; 
	double p0; 
	double pl;
	double Q; 
};


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

/// @brief Bernoulli_equation - класс для решения задач из блока 2 - Реализация стационарных моделей с прицелом на квазистационар
class Bernoulli_equation {
public:
	/// @brief pressure_p0 - Метод, рассчитывающий давление в начале участка нефтепровода [МПа]
	/// @param pipiline_parameters_XX - структура парметров трубопровода
	/// @param oil_parameters_XX - структура парметров нефти
	/// @param lambda - коэффициент гидравлического сопротивления
	/// @param v - скорость течения нефти в системе СИ
	/// @param d - внутренний диаметр трубы в системе СИ
	/// @return oil_parameters_XX.p0 - давление в начале участка нефтепровода [МПа]
	double pressure_p0(Pipiline_parameters pipiline_parameters_XX, Oil_parameters oil_parameters_XX, double lambda, double v, double d) {
		oil_parameters_XX.p0 = (oil_parameters_XX.ro * 9.81) * (oil_parameters_XX.pl * 1000000 / (oil_parameters_XX.ro * 9.81) - pipiline_parameters_XX.z0 + pipiline_parameters_XX.zl + (((lambda * pipiline_parameters_XX.l * 1000) / d) * pow(v, 2)) / (2 * 9.81)) * 0.000001;
		return oil_parameters_XX.p0;
	}
	/// @brief diameter - метод, рассчитывающий внутренний диаметр трубы
	/// @param pipiline_parameters_XX - структура парметров трубопровода
	/// @return d - внутренний диаметр трубы [м]
	double diameter(Pipiline_parameters pipiline_parameters_XX) {
		double d = (pipiline_parameters_XX.D - 2 * pipiline_parameters_XX.ds) / 1000;
		return d;
	}
	
	/// @brief  relative_roughness - метод, рассчитывающий относительную эквивалентная шероховатость
	/// @param pipiline_parameters_XX - структура парметров трубопровода
	/// @param d - внутренний диаметр трубы [м]
	/// @return e oтносительная эквивалентная шероховатость
	double relative_roughness(Pipiline_parameters pipiline_parameters_XX, double d) {
		double e = (pipiline_parameters_XX.delta / 1000) / d;
		return e;
	}

	/// @brief reynolds_number - метод, рассчитывающий число Рейнольдса, где nu переведено в систему СИ
	/// @param oil_parameters_XX - структура парметров нефти
	/// @param v - cкорость течения нефти в системе СИ
	/// @param d - внутренний диаметр трубы в системе СИ
	/// @return Re - число Рейнольдса
	double reynolds_number(Oil_parameters oil_parameters_XX, double v, double d) {
		double Re = v * d / (oil_parameters_XX.nu * 0.000001);
		return Re;
	}
	
	/// @brief flow rate - метод, рассчитывающий скорость по заданному расходу нефти
	/// @param oil_parameters_XX - структура парметров нефти
	/// @param d - внутренний диаметр трубы [м]
	/// @return v - cкорость течения нефти в системе СИ
	double flow_rate(Oil_parameters oil_parameters_QP, double d) {
		double v = (4 * oil_parameters_QP.Q / 3600) / (3.14 * pow(d, 2));
		return v;
	}	

	/// @brief flow_pressure - метод, рассчитывающий скорость в задаче PP
	/// @param pipiline_parameters_XX - структура парметров трубопровода
	/// @param oil_parameters_XX - структура парметров нефти
	/// @param lambda - коэффициент гидравлического сопротивления
	/// @param d - внутренний диаметр трубы в системе СИ
	/// @param oil_parameters_XX.p0 - давление в начале участка нефтепровода [МПа]
	/// @param oil_parameters_XX.pl - давление в конце участка нефтепровода [МПа]
	/// @return v - cкорость течения нефти в системе СИ	
	double flow_pressure(Pipiline_parameters pipiline_parameters_XX, Oil_parameters oil_parameters_XX, double lambda, double d) {
		double v = pow((2 * 9.81 * d / pipiline_parameters_XX.l / 1000 * ((oil_parameters_XX.p0  - oil_parameters_XX.pl) * 1000000 / (oil_parameters_XX.ro * 9.81) + pipiline_parameters_XX.z0 - pipiline_parameters_XX.zl) / lambda), 0.5);
		return v;
	}

	/// @brief volume_flow - метод, рассчитывающий объемный расход
	/// @param v - скорость течения нефти в системе СИ
	/// @param d - внутренний диаметр трубы [м]
	/// @return Q - объемный расход [м^3/ч]
	double volume_flow(double v, double d) {
		double Q = 3.14 * pow(d, 2) * v * 3600 / 4;
		return Q;
	}
};

/// @brief Iterative_solutions - класс для решения задач из блока 2 численными методами
class Iterative_solutions {
public:
	/// @brief method_eyler - численное интегрирование дифференциального уравнения методом Эйлера
	/// @param pipiline_parameters_XX_Eyler - cтруктура парметров трубопровода
	/// @param Oil_parameters - cтруктура парметров нефти
	/// @param tw - касательное напряжение трения, учитывающее трение жидкости при течении по трубе
	/// @param n - кол-во точек расчётной сетки
	/// @param h - шаг по по координате расчетной сетки [м] 
	/// @param p_prev - давление на предыдущей итерации (граничное условие) [МПа]
	/// @return p_current - давление на текущей итерации (рассчитанное значение) [МПа]
	double method_eyler(Pipiline_parameters pipiline_parameters_XX, Oil_parameters oil_parameters_XX, double tw, int n, double h, double d, double p_prev) {
		for (size_t i = 1; i <= n; i++) {
			double p_current = (p_prev * 1000000 - h * (-4 * tw / d - oil_parameters_XX.ro * 9.81 * (pipiline_parameters_XX.zl - pipiline_parameters_XX.z0) / (n - 1) / h))/1000000;
			return p_current;
		}
	}
	/// @brief 
	/// @return 
	double method_newton() {

	}
};

/// @brief Функция решения задачи QP
/// @param Pipiline_parameters - cтруктура парметров трубопровода
/// @param Oil_parameters - cтруктура парметров нефти
/// @return p0 - давление в начале участка нефтепровода [МПа]
double solver_QP(Pipiline_parameters pipiline_parameters_QP, Oil_parameters oil_parameters_QP) {
	
	// Объявляем объект lambda_QP класса Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient lambda_QP;
	// Объявляем объект task_QP класса Bernoulli_equation
	Bernoulli_equation task_QP;
	// d - внутренний диаметр трубы[м]
	double d = task_QP.diameter(pipiline_parameters_QP);
	// e - относительная эквивалентная шероховатость
	double e = task_QP.relative_roughness(pipiline_parameters_QP, d);
	// v - cкорость течения нефти в системе СИ
	double v = task_QP.flow_rate(oil_parameters_QP, d);
	// Re - число Рейнольдса
	double Re = task_QP.reynolds_number(oil_parameters_QP, v, d);

	// Формула Стокса
	if (Re < 2000) {
		double lambda = lambda_QP.stokes_formula(Re);
		return task_QP.pressure_p0(pipiline_parameters_QP,oil_parameters_QP, lambda, v, d);
	}
	// Формула Блазиуса
	else if (Re >= 2000 && Re <= 4000) {
		double lambda = lambda_QP.blasius_formula(Re);
		return task_QP.pressure_p0(pipiline_parameters_QP, oil_parameters_QP, lambda, v, d);
	}
	// Формула Альтшуля
	else {
		double lambda = lambda_QP.altschul_formula(Re, e);
		return task_QP.pressure_p0(pipiline_parameters_QP, oil_parameters_QP, lambda, v, d);
	}
}

/// @brief Функция решения задачи PP
/// @param Pipiline_parameters - cтруктура парметров трубопровода
/// @param Oil_parameters - cтруктура парметров нефти
/// @return 
double solver_PP(Pipiline_parameters pipiline_parameters_PP, Oil_parameters oil_parameters_PP) {
	// labmda_previous - гидравлическое сопротивление на предыдущем приближении
	double labmda_previous = 0.021;
	// labmda_current - гидравлическое сопротивление на текущем приближении
	double labmda_current = 0.02;
	// eps - допустимая погрешность
	double eps = 0.0005;
	// Объявляем объект task_PP класса Bernoulli_equation
	Bernoulli_equation task_PP;
	// d - внутренний диаметр трубы[м]
	double d = task_PP.diameter(pipiline_parameters_PP);
	// Объявляем объект lambda_PP класса Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient lambda_PP;
	// v - cкорость течения нефти в системе СИ
	double v = task_PP.flow_pressure(pipiline_parameters_PP, oil_parameters_PP, labmda_current, d);
	int i = 1;
	while (abs(labmda_previous - labmda_current) >= eps) {
		
		// Re - число Рейнольдса
		double Re = task_PP.reynolds_number(oil_parameters_PP, v, d);
		// e - относительную эквивалентная шероховатость
		double e = task_PP.relative_roughness(pipiline_parameters_PP, d);
		// Формула Стокса
		if (Re < 2000) {
			labmda_previous = labmda_current;
			labmda_current = lambda_PP.stokes_formula(Re);
		}
		// Формула Блазиуса
		else if (Re >= 2000 && Re <= 4000) {
			labmda_previous = labmda_current;
			labmda_current = lambda_PP.blasius_formula(Re);
		}
		// Формула Альтшуля
		else {
			labmda_previous = labmda_current;
			labmda_current = lambda_PP.altschul_formula(Re, e);
		}
		v = task_PP.flow_pressure(pipiline_parameters_PP, oil_parameters_PP, labmda_current, d);
		i += 1;
	}
	double Q = task_PP.volume_flow(v, d);
	return Q;
}

double solver_QP_Eyler(Pipiline_parameters pipiline_parameters_QP_Eyler, Oil_parameters oil_parameters_QP_Eyler, int n, double h) {
	// p_prev - давление на предыдущей итерации (граничное условие) [МПа]
	double p_prev = oil_parameters_QP_Eyler.pl;
	// p_current - давление на текущей итерации(рассчитанное значение) [МПа]
	double p_current ;
	// Объявляем объект task_PP класса Bernoulli_equation
	Bernoulli_equation task_QP_Eyler;
	// Объявляем объект lambda_QP_Eyler класса Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient lambda_QP_Eyler;
	// Объявляем объект iterativetask_QP_Eyler класса Iterative_solutions
	Iterative_solutions iterativetask_QP_Eyler;
	// d - внутренний диаметр трубы[м]
	double d = task_QP_Eyler.diameter(pipiline_parameters_QP_Eyler);
	// v - cкорость течения нефти в системе СИ
	double v = task_QP_Eyler.flow_rate(oil_parameters_QP_Eyler, d);
	// e - относительная эквивалентная шероховатость
	double e = task_QP_Eyler.relative_roughness(pipiline_parameters_QP_Eyler, d);
	// Re - число Рейнольдса
	double Re = task_QP_Eyler.reynolds_number(oil_parameters_QP_Eyler, v, d);
	double lambda;
	// Формула Стокса
	if (Re < 2000) {
		lambda = lambda_QP_Eyler.stokes_formula(Re);
	}
	// Формула Блазиуса
	else if (Re >= 2000 && Re <= 4000) {
		lambda = lambda_QP_Eyler.blasius_formula(Re);
	}
	// Формула Альтшуля
	else {
		lambda = lambda_QP_Eyler.altschul_formula(Re, e);
	}
	// tw - касательное напряжение трения, учитывающее трение жидкости при течении по трубе
	double tw = lambda / 8 * oil_parameters_QP_Eyler.ro * pow(v, 2);
	for (size_t i = 1; i <= n; i++){
		p_current = iterativetask_QP_Eyler.method_eyler(pipiline_parameters_QP_Eyler, oil_parameters_QP_Eyler, tw, n, h, d, p_prev);
		p_prev = p_current;
	}
	return p_current;
	
	
}

TEST(Block_2, Task_QP) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipiline_parameters_QP
	Pipiline_parameters pipiline_parameters_QP = { 720, 10, 100, 50, 0.15, 80 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_QP
	Oil_parameters oil_parameters_QP = { 870, 15, 0, 0.6, 3500 };
	double p0 = round(100*solver_QP(pipiline_parameters_QP, oil_parameters_QP))/100;
	EXPECT_EQ(6.03, p0);
}

TEST(Block_2, Task_PP) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipiline_parameters_PP
	Pipiline_parameters pipiline_parameters_PP = { 720, 10, 50, 100, 0.15, 80 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_PP
	Oil_parameters oil_parameters_PP = { 870, 15, 5, 0.8};
	double Q = solver_PP(pipiline_parameters_PP, oil_parameters_PP);
	EXPECT_EQ(2739, round(Q-6));
}

TEST(Block_2, Task_QP_Eyler) {
	/// Объявление структуры с именем Pipeline_parameters для переменной pipiline_parameters_QP
	Pipiline_parameters pipiline_parameters_QP_Eyler = { 720, 10, 100, 50, 0.15, 80 };
	/// Объявление структуры с именем Oil_parameters для переменной oil_parameters_QP
	Oil_parameters oil_parameters_QP_Eyler = { 870, 15, 0, 0.6, 3500 };
	// n - кол-во точек расчётной сетки
	int n = 1000;
	// h - шаг по по координате расчетной сетки [м]
	double h = pipiline_parameters_QP_Eyler.l * 1000 / n;
	double pk = round(100*solver_QP_Eyler(pipiline_parameters_QP_Eyler, oil_parameters_QP_Eyler, n, h))/100;
	EXPECT_EQ(6.03, pk);
} 