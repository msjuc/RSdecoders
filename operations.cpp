#include "header.h"

vector<bool> *makegf (int n) {
	vector<bool>*gf = new vector<bool>[n];
	int gfi = 0;
	while (gfi < n)
	{
		while (gfi<m) {
			gf[gfi].assign(m, 0);
			gf[gfi][gfi] = true;
			gfi++;
		}

		gf[gfi] = gf[gfi - m] ^ gf[gfi - m + 1];
		gfi++;
	}
	return gf;
}
void print_vector(vector<int>v) {		//печатает вектор
	int *p = &v[0];
	while (p != &v.back() + 1)cout << *p++ << ",";
}
void print_bector(vector<bool>v) {		//печатает двоичный вектор
	for (int i = 0; i < v.size(); i++)cout << v[i] ? 1 : 0;
}
vector<bool> operator ^ (vector<bool> a, vector<bool> b) {	//сложение двоичных векторов по модулю 2
	int min;
	vector<bool> result{};
	if (a.size() <= b.size()) { 
		min = a.size();
		result = b;
	}
	else {
		min = b.size();
		result = a;
	}
	int i = 0;
	while (i < min) {
		result[i] = a[i] ^ b[i];
		i++;
	}
	return result;
}; 
int in_GF_range(int power) {		//вход - степень альфа (элемента GF), не число! (прим: число = степень + 1)
	while (power > n-1) power -= n;
	while (power < 0)	 power += n;
	return power;
}
vector<bool> operator / (vector<bool> a, vector<bool> b) {	//деление двоичных векторов (элементов GF)
	return A(in_GF_range(L(a) - L(b))+1);
}
vector<bool> operator * (vector<bool> a, vector<bool> b) {
	if (L(a) == 0 || L(b) == 0)
		return A(0);
	return A(in_GF_range(L(a) + L(b) - 2)+1);
}
vector<int> cut(vector<int>a) {					//удаление старших нулей (нулевых коэффициентов при Х)
	if (!a.empty()) 
		while (!a.empty() && a.back() == 0)a.pop_back();
	return a;
}
vector<int> operator + (vector<int> a, vector<int> b) {		//сложение полиномов (где коэффициенты - элементы GF)
	vector<int>result{};
	int min;
	if (a.size() <= b.size()) {				//количество итераций = длине минимального вектора (т.е. минимальной старшей степени)
		min = a.size();
		result = b;
	}
	else {
		min = b.size();
		result = a;
	}
	for (int i = 0; i < min; i++) result[i] = (L(A(a[i]) ^ A(b[i])));
	return result;
}
vector<int> operator * (vector<int> a, vector<int> b) {		//перемножение полиномов
	vector<int>result{};
	result.assign(a.size() + b.size() - 1, 0);
	for (int i = 0; i < a.size(); i++)
		for (int j = 0; j < b.size(); j++)
			result[i + j] = L((A(a[i]) * A(b[j])) ^ A(result[i + j]));
	return result;
}
vector<int> operator / (vector<int> a, vector<int> b) {		//деление полиномов (целая часть)
	vector<int> res_polynom{};
	res_polynom.assign(a.size() - b.size() + 1, 0);
	while (a.size() >= b.size()) {
		a = cut(a); //detele old zeros!
		if (a.empty())break;	//если а был равен нулю - выход
		int old = a.size() - b.size(); //степень монома результирующего полинома
		if (old < 0)break;		//степень монома меньше нуля? выход!
		int mon = L(A(a.back()) / A(b.back()));	//коэффициент монома
		vector<int> monom{};
		monom.assign(old, 0);
		monom.push_back(mon);
		res_polynom[old] = mon;
		vector<int> btw = b * monom;
		a = a + btw;
	}
	return res_polynom;
}
vector<int> operator % (vector<int> a, vector<int> b) {		//деление полиномов (остаток)
	vector<int> remain = a;
	while (remain.size() >= b.size()) {
		remain = cut(remain); //detele old zeros!
		int old = remain.size() - b.size(); //degree of monom
		if (old < 0)break;
		int mon = L(A(remain.back()) / A(b.back())); //coefficient of monom
		vector<int> monom{};
		monom.assign(old, 0);
		monom.push_back(mon);//monom[old] = mon;
		vector<int> btw = b * monom;
		remain = remain + btw;
	}
	return remain;
}
vector <int> shift_n_zeros_to(vector<int>a, int n) {	//умножить полином на х^n
	vector<int> result{};
	result.assign(n, 0);
	for (int i = 0; i < a.size(); i++)result.push_back(a[i]);
	return result;
}
vector<int> return_OM() {				//вычисление ОМ
	vector<int> OM{ b+1, 1 };
	for (int i = b+1; i < dmin - 1; i++)OM = OM * vector<int>{i + 1, 1};
	return OM;
}
vector<bool> A(int num) {				//сопоставляет числу двоичный вектор из GF
	if (num)					//если число не нуль - возврат альфа в степени (число-1)
		return gf[num - 1];
	else {						//иначе - возврат нулевого вектора
		vector<bool>zero{};
		zero.assign(m, 0);
		return zero;
	}
}
int L(vector<bool> v) { 				//сопоставляет двоичному вектору из gf число (функция, обратная A())
	for (int i = 0; i < n; i++)if (v == gf[i]) return i + 1;
	return 0; //NULL
}
int Fx(vector<int> f, int x) { //где X - число, а не степень альфа. подставляет элемент GF в полином f
	int siz = f.size();
	vector<bool> result = A(f[0]);
	for (int i = 1; i < siz; i++) {
		result = result ^ (A(f[i]) * pow_gf(x,i));
	}
	return L(result);
}
vector<bool> pow_gf (int a, int power) {		//возводит число-элемент GF a в степень power
	return A(in_GF_range((a - 1)*power) + 1);
}
vector<int> return_syndroms(vector<int> code) {		//возвращает синдромный полином на основе принятого полиномиального кода
	vector<int> result{ 1 };
	int checkbit = 0;
	int temp;
	for (int i = b; i < dmin - 1; i++) {
		temp = Fx(code, L(gf[i]));
		result.push_back(temp);
		checkbit += temp;
	}
	result.push_back(checkbit);
	return result;
}
vector<int> return_Z(vector<int>S, vector<int> C) {	//полином Z для альтернативного алг. Форни (вход - синдромы и локаторы ошибок)
	vector<int> result{ 1 };
	vector<bool> temp;
	int size = C.size();
	for (int j = 1; j < size; j++) {	//C.size == v+1
		temp = A(S[j]) ^ A(C[j]);
		for (int i = 2; i < j + 1; i++) {
			temp = temp ^ (A(S[j - i + 1]) * A(C[i-1])); //S[ j - (i-1) ] * C[ i-1 ]
		}
		result.push_back(L(temp));
	}
	return result;
}
vector<int> BMA(vector<int> S) {		//алгоритм БМ
	vector<int>	Cix{ 1 },		//sigma
			Px{ 0, 1 };	//corrector
	vector<bool>	sum;
	vector<int>	dv, Cix_new;
	int		d, l = 0, i = 1;
	{
		while (1) {
			sum = A(0); //двоичный вектор нуля (вернет 0000 для m=4)
			for (int j = 1; j <= l; j++)sum = sum ^ (A(Cix[j]) * A(S[i - j]));	// вычисление d
			d = L(A(S[i]) ^ sum);							//
			if (d) {
				dv = { d };
				Cix_new = Cix + dv * Px;
				if (2 * l < i) {
					l = i - l;
					Px = Cix / dv;
				}
				Cix = Cix_new;
			}
			//вывод:
			/*printf("----\ni=%d: d=%d,l=%d\nsigma: ", i, d, l);
			print_vector(Cix); cout << "\np(old): ";
			print_vector(Px); cout << "\np(new): ";*/
			Px = vector<int>{ 0,1 }*Px;
			/*print_vector(Px); cout << "\n---\n";*/

			i++;
			if (i >= dmin)break;
		}
		cout << "sigma: "; print_vector(Cix); cout << endl;
		sum.~vector();
		dv.~vector();
		Cix_new.~vector();
		Px.~vector();
	}
	return Cix;
}
vector<int> CHEN(vector<int> C) {				//метод Ченя
	vector<int>	locators{};
	int		loc_at;
	for (int j = 1; j < n + 1; j++) {
		loc_at = Fx(C, j);
		if (!loc_at) {
			loc_at = in_GF_range(-(j - 1));
			locators.push_back(loc_at);
		}
	}
	return locators;
}
vector<int> FORNEY(vector<int>j, vector<int> Zx) {	//альтернативный алг. Форни
	int		jl, alpha, size = j.size();
	vector<bool>	e1, e2, e3;
	vector<int>	E; E.assign(n, 0);
	for (int l = 1; l < size + 1; l++) {
		jl = j[l - 1];
		alpha = in_GF_range(-jl) + 1;
		e1 = gf[jl * (1 - b)];
		e2 = A(Fx(Zx, alpha));
		e3 = A(1);
		for (int i = 1; i < size + 1; i++) {
			if (i == l)continue;
			e3 = e3 * (A(1) ^ A(in_GF_range(j[i - 1] - jl) + 1));
		}
		E[jl] = L((e1*e2) / e3);
	}
	e1.~vector();
	e2.~vector();
	e3.~vector();
	return E;
}
vector<bool> det(vector<int> *mtx) { //для матриц 1х1 в 1 итерацию. иначе - рекурсия!
	int size = mtx->size();
	int offset;
	vector<bool> result = A(0);
	vector<bool> det_alcomp=A(1);
	vector<vector<bool>>rangesum;
	for (int eri = 0; eri < size; eri++) { //er-th row, i-th column
		if (size > 1) {
			vector<int> *mtz = new vector<int>[size - 1]; //алгебраическое дополнение
			for (int i = 0; i < size-1; i++) {
				mtz[i].assign(size - 1, 0);
				offset = 0;
				for (int j = 0; j < size-1; j++) {
					if (j >= eri)offset = 1;
					mtz[i][j] = mtx[i+1][j+ offset];
				}
			}
			det_alcomp = A(mtx[0][eri]) * det(mtz);
		}
		else return A(mtx[0][0]);
		rangesum.push_back(det_alcomp);
	}
	for (int i = 0; i < size; i++)result = result ^ rangesum[i];
	return result;
}
vector<int> PGZ(vector<int> S) {
	int i = t;		//изначально, i равно максимальному числу ошибок, которые код может исправить
	vector<bool> dt;	//здесь будет значение определителя в двоичной форме (формирование конечного поля с p(x)=x^m+x+1 см. функция makegf
	vector<int> Cix;	//сигма - полином локаторов ошибок
	while (i > 0) {
		vector<int> *mtx = new vector<int>[i];	//матрица синдромов
		for (int j = 0; j < i; j++) {		//заполнение...
			mtx[j].assign(i, 0);
			for (int k = 0; k < i; k++)
				mtx[j][k] = S[j + k + 1];
		}
		dt = det(mtx);	//высилили определитель
		if (L(dt)) {	//если логарифм вектора определителя != 0...
			Cix.assign(i+1, 0); Cix[0] = 1;	//строим сигму
			vector<int> *mtx1 = mtx;	//используем для каждого монома при x^j используется dt и определитель матрицы,
			for (int j = i; j > 0; j--) {	//(i-j)-ый столбец которого равен вектору-столбцу из правой части ключевого уравнения. 
				for (int k = i + 1; k < 2 * i + 1; k++)mtx1[k-i-1][i - j] = S[k];
				Cix[j] = L( pow_gf(L(dt),-1) * det(mtx1) );
				mtx1 = mtx;
			}
			return Cix;
		}
		i--;
	}
	return Cix = { 1 };
}
