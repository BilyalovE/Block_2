#pragma once
#include <cmath>
#include "const.h"
/// @brief Pipeline_parameters - Структура парметров трубопровода
/// @param D - внешний диаметр [м]
/// @param ds - толщина стенки [м]
/// @param z0 - высота в начале участка трубопровода [м]
/// @param zl - высотка в конце участка трубопровода [м]
/// @param delta - абсолютная шероховатость трубы [м]
/// @param l - длина трубопровода [м]
/// @param Q - расход [м^3/с]
/// @param nu - кинематическая вязкость [Ст]
struct Pipeline_parameters {
	// D - внешний диаметр [м]
	double D;
	// ds - толщина стенки [м]
	double ds;
	// z0 - высота в начале участка трубопровода [м]
	double z0;
	// zl - высотка в конце участка трубопровода [м]
	double zl;
	// delta - абсолютная шероховатость трубы [м]
	double delta;
	// l - длина трубопровода [м]
	double l;
	// Q - расход [м^3/с]
	double Q;
	// nu - кинематическая вязкость [Ст]
	double nu;
	// Расчет внутреннего диамера
	double get_inner_diameter() const {
		return D - 2 * ds;
	}
	// Расчет относительной шероховатости
	double get_relative_roughness() const {
		double inner_diameter = get_inner_diameter();
		return delta / inner_diameter;
	}
	// Расчет площади сечения
	double get_inner_square() const {
		double inner_diameter = get_inner_diameter();
		return k_pi * inner_diameter * inner_diameter / 4;
	}
	// Расчет скорости по расходу
	double speed_flow() const {
		double inner_square = get_inner_square();
		return Q / inner_square;
	}
	// Расчет числа Рейнольдса
	double reynolds_number() const {
		double inner_diameter = get_inner_diameter();
		double v = speed_flow();
		return v * inner_diameter / nu;
	}



};

/// @brief Oil_parameters - Структура парметров нефти
struct Oil_parameters {
	/// @param ro - плотность нефти [кг/м^3]
	double ro;
	/// @param p0 - давление в начале участка нефтепровода [Па]
	double p0;
	/// @param pl - давление в конце участка нефтепровода [Па]
	double pl;
};
