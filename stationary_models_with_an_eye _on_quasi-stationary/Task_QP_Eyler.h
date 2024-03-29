#pragma once
#include "const.h"
#include "struct.h"
#include "Bernoulli_equation.h"
#include "Hydraulic_resistance_coefficient.h"

/// @brief Task_QP_Eyler - Класс для решения задачи QP численным методом Эйлера
class Task_QP_Eyler
{
	// Поля класса
	// pipeline_parameters - структура параметров трубопровода
	Pipeline_parameters pipeline_parameters;
	// oil_parameters - структура параметров нефти
	Oil_parameters oil_parameters;
	// tw - касательное напряжение трения, учитывающее трение жидкости при течении по трубе
	double tw;
	// n - кол - во точек расчётной сетки
	int n;
	// h - шаг по по координате расчетной сетки [м]
	double h;
	// internal_diameter - внутренний диаметр трубы [м]
	double internal_diameter;
	// pressure_previous - давление на предыдущей итерации граничное условие) [Па]
	double pressure_previous;
	// v - скорость течения нефти [м/с]
	double v;

public:
	/// @brief Конструктор класса
	/// @param pipeline_parameters - структура параметров трубопровода
	/// @param oil_parameters - структура параметров нефти
	/// @param tw - касательное напряжение трения, учитывающее трение жидкости при течении по трубе
	/// @param n - кол - во точек расчётной сетки
	/// @param h - шаг по по координате расчетной сетки [м]
	/// @param internal_diameter - внутренний диаметр трубы [м]
	/// @param p_previous - давление на предыдущей итерации граничное условие) [Па]
	Task_QP_Eyler(const Pipeline_parameters& pipeline_parameters, const Oil_parameters& oil_parameters, double v );
	// Методы класса

	/// @brief Солвер для решения задачи QP численным методом Эйлера
	/// @return pressure_current - давление на текущей итерации(рассчитанное значение)[Па]
	double solver_eyler(double v);

};


 