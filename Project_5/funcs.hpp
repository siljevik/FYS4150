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

	arma::cx_vec vector_filler(int M, arma::cx_mat V);

	void index_translator(int M, int k, int & i_p, int & j);

	void matrix_filler(int M, arma::cx_double r_val, int L, arma::cx_mat & A, arma::cx_mat & B);

	void diagonal_fill_AB(int M, arma::cx_double h, arma::cx_double dt, int L, arma::cx_mat V,arma::cx_mat & A, arma::cx_mat & B);

	arma::cx_vec Bu_b(int M, int L, arma::cx_mat V, arma::cx_mat B);

	arma::cx_vec Au_b(arma::cx_mat A, arma::cx_vec b);

	//void initial_u(int M, arma::cx_double h, int L, arma::cx_vec u_0);

  void double_slit(arma::cx_double h, arma::cx_double v0, arma::sp_cx_mat V)

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