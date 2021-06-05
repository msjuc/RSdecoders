#include "header.h"
#include <random>
#include <time.h>
/*BCE7qA npoBep9u Ha BIG ENDIAN!!!!*/

int main() {
	//for checking!
	random_device rd{};
	mt19937 gen{ rd() };
	double Q;
	srand(time(NULL));
	cout << "input the \'noise level\': "; cin >> Q;
	normal_distribution<> nd{ 0,Q };
	///R = { 0,0,2,0,6 };
	///goto decoding;
	int N = 0;
	int good = 0, noerr=0;
	vector<int> OM = return_OM(); cout << "OM: "; print_vector(OM); cout << endl;
begin:
	printf("\n---%d---\n", N);
	vector<int>  U, Uv, V, R;
	vector<bool>V_bin;
create_word: {
	U = {}, V = {}, R = {}, Uv = {}, V_bin = {};
	for (int i = 0; i < k; i++) {
		U.push_back( rand() % n + 1 );
	}
	cout << "U: "; print_vector(U); cout << endl;
}
			/* she makes the sound the sound sea make
			 to calm me down*/
make_code_word: {
	 Uv = shift_n_zeros_to(U, n - k);
	 V = Uv + Uv % OM;  /*code message*/
	 cout << "V: "; print_vector(V); cout << endl;
	 U.~vector(); Uv.~vector();
}
	/*sending through the channel*/
transmit_word: {
	V_bin = {}; R = {};
	//oTo6pa#eHue B qBou4Hb1u vector
	for (int i = 0; i < n; i++) {
		vector<bool> Vi = A(V[i]);
		for (int j = 0; j < m; j++)
		{
			V_bin.push_back(Vi[j]);
		}
	}
	///cout << "Vb: "; print_bector(V_bin); cout << endl;
	cout << "...\n";
	//distorting
	double		distortion;
	for (int i = 0; i < n*m; i++) {
		//bpsk
		distortion = V_bin[i] ? 1 : -1 + nd(gen);
		V_bin[i] = distortion >= 0 ? 1 : 0;
	}
	///cout << "Vb: "; print_bector(V_bin); cout << endl;
	for (int i = 0; i < n; i++) {
		vector<bool>temp{};
		for (int j = 0; j < m; j++) {
			temp.push_back(V_bin[i*m + j]);
		}
		R.push_back(L(temp));
		temp.~vector();
	}
	cout << "R: "; print_vector(R); cout << endl;
	V_bin.~vector();
}
detecting:
		 int		isOK;
vector<int>		Syns = return_syndroms(R);
cout << "S:"; print_vector(Syns); cout << endl;
if (!Syns[dmin]) { cout << "no errors!\n"; noerr++;  goto end0; }
{
	vector<int>		Cix, locators, Zx, E{};
	Cix = BMA(Syns);	//Алгоритм Берлекэмпа-Мэсси для РС
								//TODO: Реализовать Евклида.
//	Cix = EA	(Syns);
	//Cix = PGZ	(Syns);		//Алгоритм ПГЦ для РС
	cout << "Cx:"; print_vector(Cix); cout << endl;
	locators = CHEN(Cix); //Чень - поиск позиций ошибок
		if (!locators.empty()) { cout << "j: "; print_vector(locators); cout << endl; }
		cout << "CHECK: ";
		if (locators.empty() || locators.size() < Cix.size() - 1) { 
				cout << "too much errors, can't correct\n";
			goto end0;
		}
		cout << "errors can be correcred\n";
		
		// Форни - значения ошибок:
		Zx = return_Z(Syns, Cix);
		E = FORNEY(locators, Zx);
		
		Zx.~vector();
		Cix.~vector();
		Syns.~vector();
		locators.~vector();
		//вывод:
		R = R + E;
		cout << "R': "; print_vector(R); cout << endl;
		R = R + V;
		isOK = 0;
		for (int i = 0; i < n; i++)isOK += R[i];
		if (!isOK) { cout << "all correct\n"; good++; }
		else cout << "errors remained\n";
		E.~vector();
	}
end0:
	V.~vector();
	R.~vector();
	if (++N < 20)goto begin;
	printf("%d of %d words have been corrected successfully",good, N-noerr);
}