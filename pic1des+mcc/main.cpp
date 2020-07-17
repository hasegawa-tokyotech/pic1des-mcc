//
//  main.cpp
//  es1d2v
//
//  Created by Jun Hasegawa on 2013/03/12.
//	Last Updated by Jun Hasegawa on 2019/06/13.
//
//  Copyright (c) 2019 Jun Hasegawa, Tokyo Tech. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <ctime>
#include <cstdio>
#include <cmath>

#include "World_1DES.h"

int main()
{
	using namespace std;
	
	//char fn[255];
	//cout << "input file name = ";
	//cin >> fn;
	
    try {
		//World_1DES pic_world( fn );
		World_1DES pic_world( "input.txt" );

		time_t begin, end;
		time( &begin );
    
		pic_world.run();
    
		time( &end );
		std::cout << "\ncalculation time = " << end - begin << " s\n";
		
        int n = pic_world.tstep()/pic_world.tout(); // total number of data output
        
		// output script files for gif animation
		ofstream fout( "ps_ani.txt", ios::out );
		fout << "set terminal gif animate delay 5 optimize size 640,480\n";
		fout << "set output 'ps.gif'\n";
		fout << "set xrange [ 0 : " << pic_world.lx() << " ]\n";
		fout << "set yrange [ -0.02 : 0.02 ]\n";
		fout << "do for [i=0:" << n << "] {\n";
		fout << "\tplot ";
		for ( int i = 0; i < pic_world.ns(); i++ ) {
			fout << "sprintf(\"pa" << i << "_%03d.txt\",i) u 4:5 notitle w points pt 7 ps 0.2, ";
		}
		fout << "\n}\nunset output\n";
        fout.close();
		
        fout.open( "fi_ani.txt", ios::out );
		fout << "set terminal gif animate delay 5 optimize size 640,480\n";
		fout << "set output 'fi.gif'\n";
		fout << "set xrange [ 0 : " << pic_world.lx() << " ]\n";
        char phin[10];
        sprintf( phin, "%.0e", 2*fabs(pic_world.phin()) );
		fout << "set yrange [ -" << phin << " : " << phin << " ]\n";
		fout << "do for [i=0:" << n << "] {\n";
		fout << "\tplot sprintf(\"fi0_%03d.txt\",i) u 1:3 notitle w l\n}\n";
		fout << "unset output\n";
		fout.close();
    }
	catch( exception& err )
	{
		cout << "*** Exception: " << err.what() << " ***" << endl;
		return 1;
	}
	
	return 0;
}

