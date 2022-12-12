//////////////////////////////////////////
// SEE BOTTOM FOR CHRISTMAS SURPRISE <3 //
//////////////////////////////////////////
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include <map>

#ifndef __funcs_hpp__
#define __funcs_hpp__

class funcs
{

	private: 
	// Nothing to see here.

	public:

	arma::cx_vec vector_filler(int M, arma::sp_cx_mat V);

	void index_translator(int M, int k, int & i_p, int & j);

	void matrix_filler(int M, arma::cx_double r_val, int L, arma::sp_cx_mat & A, arma::sp_cx_mat & B);

	void diagonal_fill_AB(int M, arma::cx_double h, arma::cx_double dt, int L, arma::sp_cx_mat V,
				arma::sp_cx_mat & A, arma::sp_cx_mat & B);

	arma::cx_vec Bu_b(int M, int L, arma::sp_cx_mat V, arma::sp_cx_mat B);

	arma::cx_vec Au_b(arma::sp_cx_mat A, arma::cx_vec b);

	void initial_u(int M, double h, int L, arma::cx_vec U_0, arma::cx_double x_c, arma::cx_double y_c,
                arma::cx_double sigma_x, arma::cx_double sigma_y, arma::cx_double p_x, arma::cx_double p_y);

	void initial_V();
}; // end of class Header

#endif

/*
     _  _   ___  ___   ___  _  _
    |,\/,| |[_' |[_]) |[_]) \\//
    ||\/|| |[_, ||'\, ||'\,  ||

            ___ __ __ ____  __  __  ____  _  _    __    __
           // ' |[_]| |[_]) || ((_' '||' |,\/,|  //\\  ((_'
           \\_, |[']| ||'\, || ,_))  ||  ||\/|| //``\\ ,_))
                                                               

                                         ,;7,
                                       _ ||:|,
                     _,---,_           )\'  '|
                   .'_.-.,_ '.         ',')  j
                  /,'   ___}  \        _/   /
      .,         ,1  .''  =\ _.''.   ,`';_ |
    .'  \        (.'T ~, (' ) ',.'  /     ';',
    \   .\(\O/)_. \ (    _Z-'`>--, .'',      ;
     \  |   I  _|._>;--'`,-j-'    ;    ',  .'
    __\_|   _.'.-7 ) `'-' "       (      ;'
  .'.'_.'|.' .'   \ ',_           .'\   /
  | |  |.'  /      \   \          l  \ /
  | _.-'   /        '. ('._   _ ,.'   \i
,--' ---' / k  _.-,.-|__L, '-' ' ()    ;
 '._     (   ';   (    _-}             |
  / '     \   ;    ',.__;         ()   /
 /         |   ;    ; ___._._____.: :-j
|           \,__',-' ____: :_____.: :-\
|               F :   .  ' '        ,  L
',             J  |   ;             j  |
  \            |  |    L            |  J
   ;         .-F  |    ;           J    L
    \___,---' J'--:    j,---,___   |_   |
              |   |'--' L       '--| '-'|
               '.,L     |----.__   j.__.'
                | '----'   |,   '-'  }
                j         / ('-----';
               { "---'--;'  }       |
               |        |   '.----,.'
               ',.__.__.'    |=, _/
                |     /      |    '.
                |'= -x       L___   '--,
          snd   L   __\          '-----'
                 '.____)
*/
