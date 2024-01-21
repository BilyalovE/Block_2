/*!
	\brief ���� 2 - ���������� ������������ ������� � �������� �� ��������������
	\author Bilyalov Eldar
	\version 1
	\date 21.01.2024
*/

// ���������� ����������� ����������
#include <iostream>
#include <cmath>
#include <locale.h>
#include "gtest/gtest.h"
#include <iomanip>
#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

// ���������� ������������� ���� std
using namespace std;

/// @brief Pipiline_parameters - ��������� ��������� ������������
/// @param  D - ������� ������� [��]
/// @param ds - ������� ������ [��]
/// @param z0 - ������ � ������ ������� ������������ [�]
/// @param zl - ������� � ����� ������� ������������ [�]
/// @param delta - ���������� ������������� ����� [��]
/// @param l - ����� ������������ [��]
struct Pipiline_parameters {
	double D; 
	double ds; 
	double z0; 
	double zl; 
	double delta; 
	double l; 
};



/// @brief Oil_parameters - ��������� ��������� �����
/// @param ro - ��������� ����� [��/�3]
/// @param nu - �������������� �������� [���]
/// @param pl - �������� � ����� ������� ������������ [���]
/// @param p0 - �������� � ������ ������� ������������ [���]
/// @param Q -  ������ [�^3/�]
struct Oil_parameters {
	double ro; 
	double nu; 
	double pl; 
	double p0;
	double Q; 
};

/// @brief pressure_p0 - �������, �������������� �������� � ������ ������� ������������ [���]
/// @param pipiline_parameters_QP - ��������� ��������� ������������
/// @param oil_parameters_QP - ��������� ��������� �����
/// @param lambda - ����������� ��������������� �������������
/// @param v - �������� ������� ����� � ������� ��
/// @param d - ���������� ������� ����� � ������� ��
/// @return oil_parameters_QP.p0 - �������� � ������ ������� ������������ [���]
double pressure_p0(Pipiline_parameters pipiline_parameters_QP, Oil_parameters oil_parameters_QP, double lambda, double v, double d) {
	oil_parameters_QP.p0 = (oil_parameters_QP.ro * 9.81) * (oil_parameters_QP.pl * 1000000 / (oil_parameters_QP.ro * 9.81) - pipiline_parameters_QP.z0 + pipiline_parameters_QP.zl + (((lambda * pipiline_parameters_QP.l * 1000) / d) * pow(v, 2)) / (2 * 9.81)) * 0.000001;
	return oil_parameters_QP.p0;
}

/// @brief Hydraulic_resistance_coefficient - ����� ��������������� �������������
class  Hydraulic_resistance_coefficient {
public:

	/// @brief ������� ������
	/// @param lambda - ����������� ��������������� �������������
	/// @return lambda
	double stokes_formula(double Re) {
		double lambda = 64 / Re;
		return lambda;
	}

	/// @brief ������� ��������
	/// @param lambda - ����������� ��������������� �������������
	/// @return lambda
	double blasius_formula(double Re) {
		double lambda = 0.3164 / pow(Re, 0.25);
		return lambda;
	}

	/// @brief ������� ��������
	/// @param lambda - ����������� ��������������� �������������
	/// @return lambda
	double altschul_formula(double Re, double e) {
		double lambda = 0.11 * pow(e + 68 / Re, 0.25);
		return lambda;
	}
};

/// @brief ������� ������� ������ QP
/// @param Pipiline_parameters - c�������� ��������� ������������
/// @param Oil_parameters - c�������� ��������� �����
/// @return p0 - �������� � ������ ������� ������������
double solver(Pipiline_parameters pipiline_parameters_QP, Oil_parameters oil_parameters_QP) {
	// ���������� ������� ����� � ������� ��
	double d = (pipiline_parameters_QP.D - 2 * pipiline_parameters_QP.ds)/1000; 
	// ������������� ������������� �������������
	double e = (pipiline_parameters_QP.delta / 1000) / d;
	// �������� ������� ����� � ������� ��
	double v = (4 * oil_parameters_QP.Q / 3600) / (3.14 * pow(d, 2));
	// ����� ����������, ��� nu ���������� � ������� ��
	double Re = v * d / (oil_parameters_QP.nu * 0.000001);
	
	// ��������� ������ lambda_QP ������ Hydraulic_resistance_coefficient
	Hydraulic_resistance_coefficient lambda_QP;
	
	// ������� ������
	if (Re < 2000) {
		double lambda = lambda_QP.stokes_formula(Re);
		return pressure_p0(pipiline_parameters_QP,oil_parameters_QP, lambda, v, d);
	}
	// ������� ��������
	else if (Re >= 2000 && Re <= 4000) {
		double lambda = lambda_QP.blasius_formula(Re);
		return pressure_p0(pipiline_parameters_QP, oil_parameters_QP, lambda, v, d);
	}
	// ������� ��������
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
