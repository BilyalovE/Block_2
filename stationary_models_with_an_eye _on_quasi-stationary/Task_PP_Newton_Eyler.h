#pragma once
#include "struct.h"
#include "Task_QP_Eyler.h"
#include <iomanip>
#include <fixed/fixed.h>
class Task_PP_Newton_Eyler : public Task_QP_Eyler, public fixed_system_t<1>
{
	Pipeline_parameters pipeline_parameters;
	Oil_parameters oil_parameters;

public:

	// Конструктор класса
	Task_PP_Newton_Eyler(Pipeline_parameters pipeline_parameters, Oil_parameters oil_parameters, double v);
	using fixed_system_t<1>::var_type;


	
	/// @brief residuals - ������� �������
	/// @param - v - ������� �������� (��������, [�/�])
	var_type residuals(const var_type& x) {
		double result;
		result = oil_parameters.p0 - x;
		return result;
	}

	double solver_newton_rafson(double initial_pressure_p0);
};


