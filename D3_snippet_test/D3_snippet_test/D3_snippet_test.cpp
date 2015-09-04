// D3_snippet_test.cpp : ���C�� �v���W�F�N�g �t�@�C���ł��B

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
	Point SP;//Starting Point �n�_
	D3 D;//Direction ����

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
	//Point IntersectionPoint(Vec);
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

//x���W��Ԃ�
long double Point::getX()
{
	return X;
}

//y���W��Ԃ�
long double Point::getY()
{
	return Y;
}

//z���W��Ԃ�
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
//�O�_����ʐς����߂�
long double Point::S_point(Point B , Point C)
{
	Vec AB(*this , B) , AC(*this , C);
	return AB.S_vec(AC);
}

//�l�_����̐ς����߂�
long double Point::V_point(Point B , Point C , Point D)
{
	return Vec(( *this ) , B).V_vec(Vec(( *this ) , C) , Vec(( *this ) , D));
}


//Vec

//�x�N�g���Ǝn�_����x�N�g�������
Vec::Vec(D3 Direction_ , Point SP_)
{
	D = Direction_;
	SP = SP_;
}

//�n�_�Ɗe�x�N�g���̒�������x�N�g�������
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

//�O��
Vec Vec::Cross_product(Vec Partner)
{
	return Vec(D.Y*Partner.D.Z - D.Z*Partner.D.Y , D.Z*Partner.D.X - D.X*Partner.D.Z , D.X*Partner.D.Y - D.Y*Partner.D.X);
}

//����
double Vec::Inner_product(Vec Partner)
{
	return D.X*Partner.D.X + D.Y*Partner.D.Y + D.Z*Partner.D.Z;
}

//�Q�̓_����x�N�g�������
Vec::Vec(Point A , Point B)
{
	*this = Vec(B.getX() - A.getX() , B.getY() - A.getY() , B.getZ() - A.getZ() , A);
}

//�n�_��Ԃ�
Point Vec::getSP()
{
	return SP;
}

//�I�_��Ԃ�
Point Vec::getGP()
{
	return Point(SP.getX() + D.X , SP.getZ() + D.Z , SP.getZ() + D.Z);
}

//������Ԃ�
D3 Vec::getD()
{
	return D;
}

//������Ԃ�
long double Vec::length()
{
	return sqrtl(D.X*D.X + D.Y*D.Y + D.Z*D.Z);
}

//�P�ʃx�N�g��(����1)��Ԃ�
Vec Vec::Unit_vec()
{
	return Vec(*this) / ( *this ).length();
}

//�t�x�N�g����Ԃ�(�n�_�͕ς�炸�B)
Vec Vec::Inverse_vec()
{
	return ( *this )*-1.L;
}

//�t�x�N�g����Ԃ�(�n�_�ƏI�_�̓���ւ��B)
Vec Vec::Reverse_vec()
{
	return Vec(( *this ).Inverse_vec().getD() , ( *this ).getGP());
}

//�������ǂ���
bool Vec::Vertical(Vec Partner)
{
	return ( *this ).Inner_product(Partner) == 0;
}

//���s���ǂ���
bool Vec::Parallel(Vec Partner)
{
	return ( *this ).Unit_vec().getD() == Partner.Unit_vec().getD() || ( *this ).Unit_vec().Inverse_vec().getD() == Partner.Unit_vec().getD();
}

//�Q�̃x�N�g�����������Ă��邩�ǂ���
//TODO:������
bool Vec::isIntersection(Vec Partner)
{
	int count = 0;
	Vec O = Vec(( *this ).getSP() , Partner.getSP()).Cross_product(Partner);
	Vec P = Vec(( *this ).getGP() , Partner.getSP()).Cross_product(Partner);
	Vec Q = Vec(Partner.getSP() , ( *this ).getSP()).Cross_product(*this);
	Vec R = Vec(Partner.getGP() , ( *this ).getSP()).Cross_product(*this);
	if(O.getD().X < 0 != P.getD().X < 0)
	{
		return false;
	}
	if(O.getD().Y < 0 != P.getD().Y < 0)
	{
		return false;
	}
	if(O.getD().Z < 0 != P.getD().Z < 0)
	{
		return false;
	}
	if(Q.getD().X < 0 != R.getD().X < 0)
	{
		return false;
	}
	if(Q.getD().Y < 0 != R.getD().Y < 0)
	{
		return false;
	}
	if(Q.getD().Z < 0 != R.getD().Z < 0)
	{
		return false;
	}
	return true;
}

//����n�_2�x�N�g������ʐς����߂�
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

//����n�_���������͏I�_�Ǝn�_���q�����Ă��Ȃ��Ă͂����Ȃ�
//���������̂Ƃ��낻�̔���͂Ȃ��B
//TODO �n�_�ƏI�_�̓��ꔻ��
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

	//TODO:g++�œ{���Ȃ��悤�ɂ���
	void sort()
	{
		//std::sort(VP.begin() , VP.end());
	}

};

int main()
{
	long double xa , ya , xb , yb , xc , yc , xd , yd;
	char z;
	while(cin >> xa >> z >> ya >> z >> xb >> z >> yb >> z >> xc >> z >> yc >> z >> xd >> z >> yd)
	{
		Point A(xa , ya) , B(xb , yb) , C(xc , yc) , D(xd , yd);
		if(A.S_point(B , D) + C.S_point(B , D) == B.S_point(A , C) + D.S_point(A , C))
		{
			cout << "YES" << endl;
		} else
		{
			cout << "NO" << endl;
		}
	}
}

