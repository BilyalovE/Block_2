#pragma once
#include "struct.h"
#include "const.h"
#include <cmath>
//#include "stationary_models.cpp"

/// @brief Bernoulli_equation - класс для решения задач из блока 2 - Реализация стационарных моделей 
/// с прицелом на квазистационар (Уравнение Бернулли)
class Bernoulli_equation
{
	// Поля класса

	// m_pipeline_parameters - структура параметров трубопровода
	Pipeline_parameters m_pipeline_parameters;
	// m_oil_parameters - структура параметров нефти
	Oil_parameters m_oil_parameters;
	// m_hydraulic_resistance - коэффициент гидравлическое_сопротивление (lambda)
	double m_hydraulic_resistance;
	// m_v - скорость течения нефти [м/с]
	double m_v;
	// m_d - внутренний диаметр трубы [м]
	double m_d;
	// m_relative_roughness - относительная эквивалентная шероховатость(e)
	double m_relative_roughness;
	// m_Re - число Рейнольдса
	double m_Re;
	// m_Q - объемный расход [м^3/c]
	double m_Q;
	// m_p0 - давление в начале участка нефтепровода [Па]
	double m_p0;

public:
	/// @brief конструктор класса по умолчанию Bernoulli_equation
	Bernoulli_equation(const Pipeline_parameters& pipeline_parameters, const Oil_parameters& oil_parameters);

	/// @brief конструктор класса Bernoulli_equation
	Bernoulli_equation(const Pipeline_parameters& pipeline_parameters, const Oil_parameters& oil_parameters,
		double hydraulic_resistance, double v, double d);

	void setter1(const Pipeline_parameters& pipiline_parameters, const Oil_parameters& oil_parameters);

	/// @brief setter - сеттер конструктора Bernoulli_equation
	/// @param pipeline_parameters - структура параметров трубопровода
	/// @param oil_parameters - структура параметров нефти
	/// @param hydraulic_resistance - коэффициент гидравлическое_сопротивление (lambda)
	/// @param v - скорость течения нефти [м/с]
	/// @param d - внутренний диаметр трубы [м]
	void setter2(const Pipeline_parameters& pipeline_parameters, const Oil_parameters& oil_parameters,
		double& hydraulic_resistance, const double& v, double& d);



	/// @brief pressure_p0 - Метод, рассчитывающий давление в начале участка нефтепровода [Па]
	/// @return m_p0 - давление в начале участка нефтепровода [Па]
	double pressure_p0();


	double internal_diameter();

	/// @brief  relative_roughness - метод, рассчитывающий относительную эквивалентная шероховатость
	/// @return relative_roughness - относительная эквивалентная шероховатость (e)
	double relative_roughness();

	/// @brief reynolds_number - метод, рассчитывающий число Рейнольдса, где nu переведено в систему СИ
	/// @return m_Re - число Рейнольдса
	double reynolds_number();

	double reynolds_number(double m_d);

	/// @brief speed_flow - метод, рассчитывающий скорость по заданному расходу нефти
	/// @return m_v - cкорость течения нефти [м/с]
	double speed_flow();

	/// @brief speed_pressure - метод, рассчитывающий скорость по давлению в задаче PP
	/// @return v - cкорость течения нефти в системе СИ	
	double speed_pressure();
	double speed_pressure(double m_hydraulic_resistance);

	/// @brief volume_flow - метод, рассчитывающий объемный расход
	/// @return m_Q - объемный расход [м^3/c]
	double volume_flow();
	double volume_flow(double m_d);
};