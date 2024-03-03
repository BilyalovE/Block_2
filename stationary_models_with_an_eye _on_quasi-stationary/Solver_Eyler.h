#pragma once
// ����������� ������ ��� ����������� ��������������� ������������� 
#include "Hydraulic_resistance_coefficient.h"
// ����������� ������ ��� ������� ����� �� ����� 2 - ���������� ������������ ������� 
// � �������� �� �������������� (��������� ��������)
#include "Bernoulli_equation.h"
#include "struct.h"
#include "const.h"
#include "Task_QP_Eyler.h"

/// @brief solver_QP_Eyler - функция решения задачи QP методом Эйлера
/// @param pipeline_parameters_QP_Eyler - структура параметров трубопровода
/// @param oil_parameters_QP_Eyler - структура параметров нефти
/// @return pressure_p0 - давление в начале участка нефтепровода [Па]
double solver_QP_Eyler2(const Pipeline_parameters& pipeline_parameters_QP_Eyler,
	const Oil_parameters& oil_parameters_QP_Eyler)
{
	// v - скорость течения нефти [м/с]
	double v = pipeline_parameters_QP_Eyler.speed_flow();
	// relative_roughness - относительная эквивалентная шероховатость
	double relative_roughness = pipeline_parameters_QP_Eyler.get_relative_roughness();
	// Re - число Рейнольдса
	double Re = pipeline_parameters_QP_Eyler.reynolds_number();
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
		tw, n, h, pressure_previous);
	// Вызов итеративного метода Эйлера
	pressure_p0 = task_QP_Eyler.solver_eyler();
	return pressure_p0;
}