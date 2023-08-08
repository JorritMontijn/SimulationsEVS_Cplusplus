#pragma once
#include <vector>
#include <math.h>
#include <algorithm>
#include <numeric>
#include <assert.h>
//#include <experimental/execution_policy>
//#include <experimental/numeric>

// define namespace
using namespace std;

// vMax: element-wise maximum of two vectors; include overloaded functions for doubles
struct maxCmp { double operator()(const double &x, const double &y) { return max(x, y); } };
vector<double> vMax(const vector<double> &vecIn1, const vector<double> &vecIn2)
{
	vector<double> vecOut(vecIn1.size());
	transform(vecIn1.begin(), vecIn1.end(), vecIn2.begin(), vecOut.begin(), maxCmp());
	return vecOut;
}
vector<double> vMax(const vector<double> &vecIn, const double &dblIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vMax(vecIn, vecIn2); }
vector<double> vMax(const double &dblIn, const vector<double> &vecIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vMax(vecIn2, vecIn); }


// vLT: larger than; v1 > v2
struct gtOp { double operator()(const double &x, const double &y) { return (double)(x > y); } };
vector<double> vLT(const vector<double> &vecIn1, const vector<double> &vecIn2)
{
	vector<double> vecOut(vecIn1.size());
	transform(vecIn1.begin(), vecIn1.end(), vecIn2.begin(), vecOut.begin(), gtOp());
	return vecOut;
}
vector<double> vLT(const vector<double> &vecIn, const double &dblIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vLT(vecIn, vecIn2); }
vector<double> vLT(const double &dblIn, const vector<double> &vecIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vLT(vecIn2, vecIn); }

// vST: smaller than; v1 < v2
struct ltOp { double operator()(const double &x, const double &y) { return (double)(x < y); } };
vector<double> vST(const vector<double> &vecIn1, const vector<double> &vecIn2)
{
	vector<double> vecOut(vecIn1.size());
	transform(vecIn1.begin(), vecIn1.end(), vecIn2.begin(), vecOut.begin(), ltOp());
	return vecOut;
}
vector<double> vST(const vector<double> &vecIn, const double &dblIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vST(vecIn, vecIn2); }
vector<double> vST(const double &dblIn, const vector<double> &vecIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vST(vecIn2, vecIn); }

// vEq: equal to, v1 == v2
struct eqOp { double operator()(const double &x, const double &y) { return (double)(x == y); } };
vector<double> vEq(const vector<double> &vecIn1, const vector<double> &vecIn2)
{
	vector<double> vecOut(vecIn1.size());
	transform(vecIn1.begin(), vecIn1.end(), vecIn2.begin(), vecOut.begin(), eqOp());
	return vecOut;
}
vector<double> vEq(const vector<double> &vecIn, const double &dblIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vEq(vecIn, vecIn2); }
vector<double> vEq(const double &dblIn, const vector<double> &vecIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vEq(vecIn2, vecIn); }

// perform element-wise subtraction; include overloaded functions for doubles
struct minOp { double operator()(const double &x, const double &y) { return x - y; } };
vector<double> vMinus(const vector<double> &vecIn1, const vector<double> &vecIn2)
{
	vector<double> vecOut(vecIn1.size());
	transform(vecIn1.begin(), vecIn1.end(), vecIn2.begin(), vecOut.begin(), minOp());
	return vecOut;
}
vector<double> vMinus(const vector<double> &vecIn, const double &dblIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vMinus(vecIn, vecIn2); }
vector<double> vMinus(const double &dblIn, const vector<double> &vecIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vMinus(vecIn2, vecIn); }

// perform element-wise addition; include overloaded functions for doubles
struct plusOp { double operator()(const double &x, const double &y) { return x + y; } };
vector<double> vPlus(const vector<double> &vecIn1, const vector<double> &vecIn2)
{
	vector<double> vecOut(vecIn1.size());
	transform(vecIn1.begin(), vecIn1.end(), vecIn2.begin(), vecOut.begin(), plusOp());
	return vecOut;
}
vector<double> vPlus(const vector<double> &vecIn, const double &dblIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vPlus(vecIn, vecIn2); }
vector<double> vPlus(const double &dblIn, const vector<double> &vecIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vPlus(vecIn2, vecIn); }

// perform element-wise multiplication; include overloaded functions for doubles
struct multOp { double operator()(const double &x, const double &y) { return x * y; } };
vector<double> vMult(const vector<double> &vecIn1, const vector<double> &vecIn2)
{
	vector<double> vecOut(vecIn1.size());
	transform(vecIn1.begin(), vecIn1.end(), vecIn2.begin(), vecOut.begin(), multOp());
	return vecOut;
}
vector<double> vMult(const vector<double> &vecIn, const double &dblIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vMult(vecIn, vecIn2); }
vector<double> vMult(const double &dblIn, const vector<double> &vecIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vMult(vecIn2, vecIn); }


// perform element-wise division; include overloaded functions for doubles
struct divOp { double operator()(const double &x, const double &y) { return x / y; } };
vector<double> vDiv(const vector<double> &vecIn1, const vector<double> &vecIn2)
{
	vector<double> vecOut(vecIn1.size());
	transform(vecIn1.begin(), vecIn1.end(), vecIn2.begin(), vecOut.begin(), divOp());
	return vecOut;
}
vector<double> vDiv(const vector<double> &vecIn, const double &dblIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vDiv(vecIn, vecIn2); }
vector<double> vDiv(const double &dblIn, const vector<double> &vecIn) { vector<double> vecIn2(vecIn.size()); fill(vecIn2.begin(), vecIn2.end(), dblIn); return vDiv(vecIn2, vecIn); }


// exponentiate elements in vector
struct expOp { double operator()(const double &x) { return exp(x); } };
vector<double> vExp(const vector<double> &vecIn)
{
	vector<double> vecOut(vecIn.size());
	transform(vecIn.begin(), vecIn.end(), vecOut.begin(), expOp());
	return vecOut;
}

// calculate sum over elements in vector
double vSum(const vector<double> &vecIn)
{
	return accumulate(vecIn.begin(), vecIn.end(), (double)0);
	//experimental::parallel::reduce(experimental::parallel::par,vecIn.begin(), vecIn.end());
}

// select subset by index vector
vector<double> vSelect(const vector<double> &vecIn, const vector<double> &vecIdx)
{
	// pre-allocate output vector
	vector<double> vecOut;
	vecOut.reserve(vecIn.size());

	// loop
	for (int64_t intI=0; intI < vecIn.size(); intI++)
	{
		if (vecIdx[intI] == 1) 
		{
			vecOut.push_back(vecIn[intI]);
		}
	}

	//remove excess values and return
	vecOut.shrink_to_fit();
	return vecOut;
}