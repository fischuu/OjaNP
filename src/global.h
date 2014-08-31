/* $Id: global.h,v 1.1 2008/01/25 11:47:49 ruthe Exp $ */

#ifndef GLOBAL_H
#define GLOBAL_H

//using namespace std; //df

extern bool debug;
extern bool trace;
extern bool verbose;
extern bool quiet;
extern bool adaptive;

bool parse_arguments(int argc,char** argv);

const char* program_name();

int options();
const char* option(int idx);

int parameters();
const char* parameter(int idx);
	 
bool has_option(const char* opt);
int int_option(const char* opt,int def=0);
double double_option(const char* opt,double def=0);
const char* string_option(const char* opt,const char* def="");

#endif
