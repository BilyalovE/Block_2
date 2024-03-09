#include "Task_QP_Eyler.h"

Task_QP_Eyler::Task_QP_Eyler(const Pipeline_parameters& pipeline_parameters, const Oil_parameters& oil_parameters, double v)
{
	this->oil_parameters = oil_parameters;
	this->pipeline_parameters = pipeline_parameters;
	this->v = v;
}

// Методы класса

double Task_QP_Eyler::solver_eyler()
{
	internal_diameter = pipeline_parameters.get_inner_diameter();
	// relative_roughness - относительная эквивалентная шероховатость
	double relative_roughness = pipeline_parameters.get_relative_roughness();
	// Re - число Рейнольдса
	double Re = pipeline_parameters.reynolds_number();
	// pressure_p0 - давление в начале участка нефтепровода[Па]
	double pressure_p0;
	// hydraulic_resistance - коэффициент гидравлического сопротивления (lambda)
	double hydraulic_resistance;
	// Объявляем объект lambda_QP_Eyler класса Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient lambda_QP_Eyler(Re, relative_roughness);
	// Расчет коэффициента гидравлического сопротивления (lambda)
	hydraulic_resistance = lambda_QP_Eyler.calculation_hydraulic_resistance_coefficient();
	// tw - касательное напряжение трения, учитывающее трение жидкости при течении по трубе
	double tw = hydraulic_resistance / 8 * oil_parameters.ro * pow(v, 2);
	// pressure_previous - давление на предыдущей итерации (граничное условие) [Па]
	double pressure_previous = oil_parameters.pl;
	// pressure_current - давление на текущей итерации(рассчитанное значение) [Па]
	double pressure_current;
	// n - кол-во точек расчётной сетки
	int n = 1000;
	// h - шаг по по координате расчетной сетки [м]
	double h = pipeline_parameters.l / n;
	for (size_t i = 1; i <= n; i++) {
		pressure_current = (pressure_previous - h * (-4 * tw / internal_diameter - oil_parameters.ro * k_g * (pipeline_parameters.zl - pipeline_parameters.z0) / (n - 1) / h));
		pressure_previous = pressure_current;
	}
	return pressure_current;
}

