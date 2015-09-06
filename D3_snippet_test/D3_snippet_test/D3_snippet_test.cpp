// D3_snippet_test.cpp : メイン プロジェクト ファイルです。

#include "stdafx.h"


#ifndef _GLIBCXX_NO_ASSERT
#include <cassert>
#endif
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <ciso646>
#include <climits>
#include <clocale>
#include <cmath>
#include <csetjmp>
#include <csignal>
#include <cstdarg>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

//#ifdef __GXX_EXPERIMENTAL_CXX0X__
#include <ccomplex>
#include <cfenv>
#include <cinttypes>
#include <cstdbool>
#include <cstdint>
#include <ctgmath>
#include <cwchar>
#include <cwctype>
//#endif

#include <algorithm>
#include <bitset>
#include <complex>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <list>
#include <locale>
#include <map>
#include <memory>
#include <new>
#include <numeric>
#include <ostream>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <typeinfo>
#include <utility>
#include <valarray>
#include <vector>

#include <array>
//#include <atomic>
#include <chrono>
//#include <condition_variable>
#include <forward_list>
//#include <future>
#include <initializer_list>
//#include <mutex>
#include <random>
#include <ratio>
#include <regex>
#include <system_error>
//#include <thread>
#include <tuple>
#include <typeindex>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

//#include <bits/stdc++.h>

using namespace std;

class D3
{
public:
	long double X , Y , Z , EPS;

	D3(long double , long double , long double);
	bool operator== (D3 Partner);
	bool operator==(long double);
	D3 operator*(long double);
	D3 operator/(long double);
};

class Point :private D3
{
public:
	Point(long double , long double , long double);
	long double getX();
	long double getY();
	long double getZ();

	bool operator== (Point);
	bool operator<(Point);
	bool operator>(Point);
	long double S_point(Point , Point);
	long double V_point(Point , Point , Point);
};

class Vec
{
private:
	Point SP;//Starting Point 始点
	D3 D;//Direction 方向

	Vec(D3 Direction_ , Point SP_);
public:
	Vec(long double x_ , long double y_ , long double z_ , Point);
	Vec operator+(Vec Partner);
	Vec operator*(Vec Partner);
	Vec operator*(long double ld);
	Vec operator-(Vec Partner);
	Vec operator/(Vec Partner);
	Vec operator/(long double);
	bool operator== (Vec);
	bool operator==(long double);
	Vec Cross_product(Vec);
	double Inner_product(Vec);
	Vec(Point A , Point B);
	Point getSP();
	Point getGP();
	D3 getD();
	long double length();
	Vec Unit_vec();
	Vec Inverse_vec();
	Vec Reverse_vec();
	bool Vertical(Vec);
	bool Parallel(Vec);
	bool isIntersection(Vec);
	Point IntersectionPoint(Vec);
	bool onVecPoint(Point);
	long double S_vec(Vec);
	long double V_vec(Vec , Vec);
};


//D3

D3::D3(long double x_ = 0.0L , long double y_ = 0.0L , long double z_ = 0.0L)
{
	X = x_;
	Y = y_;
	Z = z_;
	EPS = 1e-9L;
}

bool D3::operator== (D3 Partner)
{
	return abs(X - Partner.X) < EPS&&abs(Y - Partner.Y) < EPS&&abs(Z - Partner.Z) < EPS;
}

bool D3::operator==(long double ld)
{
	return ( *this ) == D3();
}

D3 D3::operator*(long double ld)
{
	return D3(X*ld , Y*ld , Z*ld);
}

D3 D3::operator/(long double ld)
{
	return  *this*( ld*-1.L );
}


//Point

Point::Point(long double x_ = 0.0L , long double y_ = 0.0L , long double z_ = 0.0L)
{
	X = x_;
	Y = y_;
	Z = z_;
}

//x座標を返す
long double Point::getX()
{
	return X;
}

//y座標を返す
long double Point::getY()
{
	return Y;
}

//z座標を返す
long double Point::getZ()
{
	return Z;
}

bool Point::operator== (Point Partner)
{
	return D3(*this) == D3(Partner);
}
bool Point::operator<(Point Partner)
{
	if(( *this ).getX() != Partner.getX())
	{
		return ( *this ).getX() < Partner.getX();
	}
	if(( *this ).getY() != Partner.getY())
	{
		return ( *this ).getY() < Partner.getY();
	}
	return ( *this ).getZ() < Partner.getZ();

}
bool Point::operator>(Point Partner)
{
	return Partner < ( *this );
}
//三点から面積を求める
long double Point::S_point(Point B , Point C)
{
	Vec AB(*this , B) , AC(*this , C);
	return AB.S_vec(AC);
}

//四点から体積を求める
long double Point::V_point(Point B , Point C , Point D)
{
	return Vec(( *this ) , B).V_vec(Vec(( *this ) , C) , Vec(( *this ) , D));
}




//Vec

//ベクトルと始点からベクトルを作る
Vec::Vec(D3 Direction_ , Point SP_)
{
	D = Direction_;
	SP = SP_;
}

//始点と各ベクトルの長さからベクトルを作る
Vec::Vec(long double x_ = 0.0L , long double y_ = 0.0L , long double z_ = 0.0L , Point SP_ = Point())
{
	*this = Vec(D3(x_ , y_ , z_) , SP_);
}

Vec Vec::operator+(Vec Partner)
{
	return Vec(D.X + Partner.D.X , D.Y + Partner.D.Y , D.Z + Partner.D.Z , SP);
}

Vec Vec::operator*(Vec Partner)
{
	return Vec(D.X * Partner.D.X , D.Y * Partner.D.Y , D.Z * Partner.D.Z , SP);
}

Vec Vec::operator*(long double ld)
{
	return Vec(D*ld , SP);
}

Vec Vec::operator-(Vec Partner)
{
	return *this + ( Partner*( -1.L ) );
}


Vec Vec::operator/(Vec Partner)
{
	return Vec(D.X / Partner.D.X , D.Y / Partner.D.Y , D.Z / Partner.D.Z , SP);
}

Vec Vec::operator/(long double ld)
{
	return *this*( 1.L / ld );
}

bool Vec::operator== (Vec Partner)
{
	return D == Partner.D&&SP == Partner.SP;
}

bool Vec::operator==(long double)
{
	return ( *this ) == Vec();
}

//外積
Vec Vec::Cross_product(Vec Partner)
{
	return Vec(D.Y*Partner.getD().Z - D.Z*Partner.getD().Y , D.Z*Partner.getD().X - D.X*Partner.getD().Z , D.X*Partner.getD().Y - D.Y*Partner.getD().X);
}

//内積
double Vec::Inner_product(Vec Partner)
{
	return D.X*Partner.getD().X + D.Y*Partner.getD().Y + D.Z*Partner.getD().Z;
}

//２つの点からベクトルを作る
Vec::Vec(Point A , Point B)
{
	*this = Vec(B.getX() - A.getX() , B.getY() - A.getY() , B.getZ() - A.getZ() , A);
}

//始点を返す
Point Vec::getSP()
{
	return SP;
}

//終点を返す
Point Vec::getGP()
{
	return Point(SP.getX() + D.X , SP.getY() + D.Y , SP.getZ() + D.Z);
}

//方向を返す
D3 Vec::getD()
{
	return D;
}

//長さを返す
long double Vec::length()
{
	return sqrtl(D.X*D.X + D.Y*D.Y + D.Z*D.Z);
}

//単位ベクトル(長さ1)を返す
Vec Vec::Unit_vec()
{
	return Vec(*this) / ( *this ).length();
}

//逆ベクトルを返す(始点は変わらず。)
Vec Vec::Inverse_vec()
{
	return ( *this )*-1.L;
}

//逆ベクトルを返す(始点と終点の入れ替え。)
Vec Vec::Reverse_vec()
{
	return Vec(( *this ).Inverse_vec().getD() , ( *this ).getGP());
}

//垂直かどうか
bool Vec::Vertical(Vec Partner)
{
	return ( *this ).Inner_product(Partner) == 0;
}

//平行かどうか
bool Vec::Parallel(Vec Partner)
{
	return ( *this ).Unit_vec().getD() == Partner.Unit_vec().getD() || ( *this ).Unit_vec().Inverse_vec().getD() == Partner.Unit_vec().getD();
}

//２つのベクトルが交差しているかどうか
bool Vec::isIntersection(Vec Partner)
{
	int count = 0;
	Vec O = Partner.Cross_product(Vec(Partner.getSP() , ( *this ).getSP()));
	Vec P = Partner.Cross_product(Vec(Partner.getSP() , ( *this ).getGP()));
	Vec Q = ( *this ).Cross_product(Vec(( *this ).getSP() , Partner.getSP()));
	Vec R = ( *this ).Cross_product(Vec(( *this ).getSP() , Partner.getGP()));
	return ( !( *this ).Parallel(Partner) ) && (
		O.getD().X*P.getD().X < D3().EPS&& O.getD().Y*P.getD().Y < D3().EPS&& O.getD().Z*P.getD().Z < D3().EPS ) && (
		Q.getD().X*R.getD().X < D3().EPS&& Q.getD().Y*R.getD().Y < D3().EPS&& Q.getD().Z*R.getD().Z < D3().EPS );
}

//ベクトル同士の交点を返す
//TODO:未実装
Point Vec::IntersectionPoint(Vec)
{
	return Point();
}

//ベクトル上にポイントが有るかどうかを返す
bool Vec::onVecPoint(Point A)
{
	Vec Test(( *this ).getSP() , A);
	if(A == ( ( *this ).getSP() ) || Test.Parallel(*this) && ( *this ).getD().X / Test.getD().X >= 1.0L)
	{
		return true;
	}
	return false;
}

//同一始点2ベクトルから面積を求める
long double Vec::S_vec(Vec B)
{
	if(( ( *this ).getSP() == B.getSP() ))
	{
		Vec G = ( *this ).Cross_product(B);
		return sqrtl(G.Inner_product(G)) / 2.L;
	}
	if(( *this ).getGP() == B.getSP())
	{
		return ( *this ).Reverse_vec().S_vec(B);
	}
	if(( *this ).getSP() == B.getGP())
	{
		return ( *this ).S_vec(B.Reverse_vec());
	}
	if(( *this ).getGP() == B.getGP())
	{
		return ( *this ).Reverse_vec().S_vec(B.Reverse_vec());
	}

	return -1.L;
}

//同一始点かもしくは終点と始点が繋がっていなくてはいけない
//しかし今のところその判定はなし。
//TODO 始点と終点の同一判定
long double Vec::V_vec(Vec B , Vec C)
{
	Vec G = ( *this ).Cross_product(B);
	return G.Inner_product(C) / 6.L;
	return 0;
}

class Points
{
private:
	vector<Point>VP;

public:
	void push(Point p)
	{
		VP.push_back(p);
	}
	bool erase(Point p)
	{
		auto now = find(VP.begin() , VP.end() , p);
		if(now != VP.end())
		{
			VP.erase(now);
			return 0;
		} else
		{
			return 1;
		}
	}

	void print()
	{
		for(auto x : VP)
		{
			cout << x.getX() << " " << x.getY() << " " << x.getZ() << endl;
		}
	}

	//TODO:g++で怒られないようにする
	void sort()
	{
		//std::sort(VP.begin() , VP.end());
	}

};

int main()
{
	int Ax , Ay , Bx , By;
	cin >> Ax >> Ay >> Bx >> By;
	Vec AB(Point(Ax , Ay) , Point(Bx , By));
	int N;
	cin >> N;
	vector<Point>vP(N);
	vector<Vec>vV(N);
	for(size_t i = 0; i < N; i++)
	{
		int x , y;
		cin >> x >> y;
		vP[i] = Point(x , y);
	}
	for(int i = 0; i < N; i++)
	{
		vV[i] = Vec(vP[i] , vP[( i + 1 ) % N]);
	}
	long long int count = 0;
	for(size_t i = 0; i < N; i++)
	{
		count += AB.isIntersection(vV[i]);
	}
	cout << count / 2 + 1 << endl;

}

