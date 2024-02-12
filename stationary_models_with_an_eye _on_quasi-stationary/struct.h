#pragma once

/// @brief Pipiline_parameters - Структура парметров трубопровода
/// @param D - внешний диаметр [м]
/// @param ds - толщина стенки [м]
/// @param z0 - высота в начале участка трубопровода [м]
/// @param zl - высотка в конце участка трубопровода [м]
/// @param delta - абсолютная шероховатость трубы [м]
/// @param l - длина трубопровода [м]
struct Pipiline_parameters {
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
};

/// @brief Oil_parameters - Структура парметров нефти
/// @param ro - плотность нефти [кг/м^3]
/// @param nu - кинематическая вязкость [Ст]
/// @param p0 - давление в начале участка нефтепровода [Па]
/// @param pl - давление в конце участка нефтепровода [Па]
/// @param Q -  расход [м^3/с]
struct Oil_parameters {
	// ro - плотность нефти [кг/м^3]
	double ro;
	// nu - кинематическая вязкость [Ст]
	double nu;
	// p0 - давление в начале участка нефтепровода [Па]
	double p0;
	// pl - давление в конце участка нефтепровода [Па]
	double pl;
	// Q - расход [м^3/с]
	double Q;
};
