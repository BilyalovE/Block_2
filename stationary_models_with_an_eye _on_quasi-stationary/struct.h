#pragma once

/// @brief Pipiline_parameters - ��������� ��������� ������������
/// @param D - ������� ������� [�]
/// @param ds - ������� ������ [�]
/// @param z0 - ������ � ������ ������� ������������ [�]
/// @param zl - ������� � ����� ������� ������������ [�]
/// @param delta - ���������� ������������� ����� [�]
/// @param l - ����� ������������ [�]
struct Pipiline_parameters {
	// D - ������� ������� [�]
	double D;
	// ds - ������� ������ [�]
	double ds;
	// z0 - ������ � ������ ������� ������������ [�]
	double z0;
	// zl - ������� � ����� ������� ������������ [�]
	double zl;
	// delta - ���������� ������������� ����� [�]
	double delta;
	// l - ����� ������������ [�]
	double l;
};

/// @brief Oil_parameters - ��������� ��������� �����
/// @param ro - ��������� ����� [��/�^3]
/// @param nu - �������������� �������� [��]
/// @param p0 - �������� � ������ ������� ������������ [��]
/// @param pl - �������� � ����� ������� ������������ [��]
/// @param Q -  ������ [�^3/�]
struct Oil_parameters {
	// ro - ��������� ����� [��/�^3]
	double ro;
	// nu - �������������� �������� [��]
	double nu;
	// p0 - �������� � ������ ������� ������������ [��]
	double p0;
	// pl - �������� � ����� ������� ������������ [��]
	double pl;
	// Q - ������ [�^3/�]
	double Q;
};
