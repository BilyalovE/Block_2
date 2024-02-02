#pragma once

#include "stationary_models_with_an_eye _on_quasi-stationary.cpp"

/// @brief Bernoulli_equation - класс для решения задач из блока 2 - Реализация стационарных моделей 
/// с прицелом на квазистационар (Уравнение Бернулли)
class Bernoulli_equation {
	// Поля класса
	// m_pipiline_parameters - структура парметров трубопровода
	Pipiline_parameters m_pipiline_parameters;
	// m_oil_parameters - структура парметров нефти
	Oil_parameters m_oil_parameters;
	// m_hydraulic_resistance - коэффициент гидравлическое_сопротивление (lambda)
	double m_hydraulic_resistance;
	// v - cкорость течения нефти [м/с]
	double m_v;
	// m_d - внутренний диаметр трубы[м]
	double m_d;

public:
	// Конструткор класса Bernoulli_equation
	Bernoulli_equation()

	/// @brief pressure_p0 - Метод, рассчитывающий давление в начале участка нефтепровода [Па]
	/// @param pipiline_parameters_XX - структура парметров трубопровода
	/// @param oil_parameters_XX - структура парметров нефти
	/// @param lambda - коэффициент гидравлического сопротивления
	/// @param v - скорость течения нефти [м/с]
	/// @param d - внутренний диаметр трубы [м]
	/// @return oil_parameters_XX.p0 - давление в начале участка нефтепровода [Па]
	double pressure_p0() {
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

	/// @brief speed_flow - метод, рассчитывающий скорость по заданному расходу нефти
	/// @param oil_parameters_XX - структура парметров нефти
	/// @param d - внутренний диаметр трубы [м]
	/// @return v - cкорость течения нефти в системе СИ
	double speed_flow(Oil_parameters oil_parameters_QP, double d) {
		double v = (4 * oil_parameters_QP.Q / 3600) / (3.14 * pow(d, 2));
		return v;
	}

	/// @brief speed_pressure - метод, рассчитывающий скорость по давлению в задаче PP
	/// @param pipiline_parameters_XX - структура парметров трубопровода
	/// @param oil_parameters_XX - структура парметров нефти
	/// @param lambda - коэффициент гидравлического сопротивления
	/// @param d - внутренний диаметр трубы в системе СИ
	/// @param oil_parameters_XX.p0 - давление в начале участка нефтепровода [МПа]
	/// @param oil_parameters_XX.pl - давление в конце участка нефтепровода [МПа]
	/// @return v - cкорость течения нефти в системе СИ	
	double speed_pressure(Pipiline_parameters pipiline_parameters_XX, Oil_parameters oil_parameters_XX, double lambda, double d) {
		double v = pow((2 * 9.81 * d / pipiline_parameters_XX.l / 1000 * ((oil_parameters_XX.p0 - oil_parameters_XX.pl) * 1000000 / (oil_parameters_XX.ro * 9.81) + pipiline_parameters_XX.z0 - pipiline_parameters_XX.zl) / lambda), 0.5);
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