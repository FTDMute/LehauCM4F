#include "DopFunctions.h"

double Function(double x) {
	return 1 - 2 * x * cos(x);
}

double DiffFunction(double x) {
	return -2 * (cos(x) - x * sin(x));
}

double f_1(double x, double y)
{
	return -cos(x + y) / (2. * exp(x * x + y * y));
}

double f_2(double x, double y)
{
	return -cos(x + y) / (2. * exp(x * x + y * y));
}

double F_1(double x, double y)

{
    return 2. * x * exp(x * x + y * y) + cos(x + y);
}

double F_2(double x, double y)
{
    return 2. * y * exp(x * x + y * y) + cos(x + y);
}

double dF_1_x(double x, double y)
{
    return exp(x * x + y * y) * (2 + 4. * x * x) - sin(x + y);
}

double dF_1_y(double x, double y)
{
    return 4. * x * y * exp(x * x + y * y) - sin(x + y);
}

double dF_2_x(double x, double y)
{
    return 4. * x * y * exp(x * x + y * y) - sin(x + y);
}

double dF_2_y(double x, double y)
{
    return exp(x * x + y * y) * (2 + 4. * y * y) - sin(x + y);
}