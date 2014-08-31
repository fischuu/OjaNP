/* $Id: misc.h,v 1.1 2008/01/25 11:47:50 ruthe Exp $ */

#ifndef MISC_H
#define MISC_H

#include <list>
#include <vector>
#include <iostream>

using namespace std;

int fact(int k);
int choices(int n,int k);
double round(double x,double prec);
bool is_file(const char* s);

void get_partitions(list<vector<int> >& P,int n);
#endif
