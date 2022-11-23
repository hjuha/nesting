#include <bits/stdc++.h>
#include "Breal.hpp"

#define EPSILON 0.000000001

long double time_since(chrono::steady_clock::time_point start) {
	return std::chrono::duration_cast<std::chrono::milliseconds>(chrono::steady_clock::now() - start).count() / 1000.;
}

mt19937 gen;

using namespace std;

inline long double log_sum_exp(const long double x1, const long double x2) {
	if (isnan(x1)) return x2;
	if (isnan(x2)) return x1;
	long double x = max(x1, x2);
	return x + log(exp(x1 - x) + exp(x2 - x));
}

template <typename T>
inline T sum(const vector<T> v) {
	T s = 0;
	for (T t : v) s += t;
	return s;
}

int weighted_random(vector<long double> weights, bool gumbel = true) {
	int n = weights.size();

	if (gumbel) { // gumbel max trick, log weights
		uniform_real_distribution<> uniform_dist(0, 1);
		int initial = 0;
		while (initial < n && isnan(weights[initial])) initial++;
		if (initial == n) return -1;
		int max_i = initial;
		long double max_weight = -log(-log(uniform_dist(gen))) + weights[initial];
		for (int i = initial + 1; i < n; i++) {
			if (isnan(weights[i])) continue;
			long double w = -log(-log(uniform_dist(gen))) + weights[i];
			if (w > max_weight) {
				max_i = i;
				max_weight = w;
			}
		}
		return max_i;
	} else {
		long double s = 0;
		for (long double d : weights) s += d;
		uniform_real_distribution<> uniform_dist(0, s);
		long double p = uniform_dist(gen);
		for (int i = 0; i < n; i++) {
			if (p <= weights[i]) return i;
			p -= weights[i];
		}
		cout<<"Error in weighted_random"<<endl;
		exit(0);
	}
}

template <typename T>
long double sis1(vector<vector<T>> v) {
	int n = v.size();
	T c[n];
	bool used[n];
	for (int j = 0; j < n; j++) {
		c[j] = 0;
		used[j] = 0;
		for (int i = 0; i < n; i++) {
			c[j] += v[i][j];
		}
	}
	long double weight = 0;
	for (int i = 0; i < n - 1; i++) {
		vector<long double> p(n);
		long double ps = NAN;
		for (int j = 0; j < n; j++) {
			if (used[j] || v[i][j] < EPSILON) {
				p[j] = NAN;
			} else {
				p[j] = (long double)(n - i - c[j]) / (long double)(n - i - v[i][j]);
				c[j] -= v[i][j];
				ps = log_sum_exp(ps, p[j]);
			}
		}
		int k = weighted_random(p);
		used[k] = 1;
		weight += log(v[i][k]) - p[k] + ps;
		if (isnan(ps)) return ps;
	}
	for (int j = 0; j < n; j++) {
		if (!used[j] && v[n - 1][j] < EPSILON) return NAN;
		else if (!used[j]) return weight + log(v[n - 1][j]);
	}
}

template <typename T>
long double sis2(vector<vector<T>> v) {
	int n = v.size();
	T c[n];
	T r[n];
	bool used[n];
	bool usedr[n];
	for (int i = 0; i < n; i++) {
		r[i] = 0;
	}
	for (int j = 0; j < n; j++) {
		c[j] = 0;
		used[j] = usedr[j] = 0;
		for (int i = 0; i < n; i++) {
			c[j] += v[i][j];
			r[i] += v[i][j];
		}
	}
	long double weight = 0;
	for (int ii = 0; ii < n - 1; ii++) {
		int i = 0;
		for (int iii = 0; iii < n; iii++) {
			if (usedr[i] || (!usedr[iii] && r[i] > r[iii])) {
				i = iii;
			}
		}
		i = ii;
		vector<long double> p(n);
		long double ps = NAN;
		for (int j = 0; j < n; j++) {
			if (used[j] || v[i][j] < EPSILON) {
				p[j] = NAN;
			} else {
				p[j] = /* log(v[i][j]) */ -log(c[j] - v[i][j]);
				c[j] -= v[i][j];
				ps = log_sum_exp(ps, p[j]);
			}
		}
		int k = weighted_random(p);
		for (int j = 0; j < n; j++) {
			r[j] -= v[j][k];
		}
		used[k] = 1;
		usedr[i] = 1;
		weight += log(v[i][k]) - p[k] + ps;
		if (isnan(ps)) return ps;
	}
	int i = 0;
	while (usedr[i]) i++;
	for (int j = 0; j < n; j++) {
		if (!used[j] && v[i][j] < EPSILON) return NAN;
		else if (!used[j]) return weight + log(v[i][j]);
	}
}

vector<int> random_weighted_order(vector<long double> r) {
	int n = r.size();
	vector<int> o(n);
	int id[n];
	for (int i = 0; i < n; i++) id[i] = i;
	for (int i = 0; i < n; i++) {
		int j = weighted_random(vector<long double>(r.begin() + i, r.end()), false);
		j += i;
		o[i] = id[j];
		swap(r[i], r[j]);
		swap(id[i], id[j]);
	}
	return o;
}

template <typename T>
long double pps(vector<vector<T>> v) {
	int n = v.size();
	vector<long double> r(n, 0);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			r[i] += v[i][j];
		}
	}
	vector<int> o = random_weighted_order(r);
	long double p = 0;
	for (int i = 0; i < n; i++) {
		if (sum(v[o[i]]) < EPSILON) return NAN;
		int j = weighted_random(v[o[i]], false);
		p += log(sum(v[o[i]]));
		for (int k = 0; k < n; k++) v[k][j] = 0;
	}
	return p;
}

// Monte Carlo (3.2), simple IS (3.2.1) and 2-dimensional (3.2.4) approaches are almost uniformly inferior (higher standard deviations) for all test matrices. – Smith & Dawkins (2001)

/*
Currently produces (natural logarithm of) estimate for the first 1000 samples and then every 0.1 seconds
Otherwise we could get millions or even billions of data points on longer runs 

INPUT (note that J can be anything as it's ignored: It was used for rejection samplers)
n J time_limit
a_11 a_12 ... a_1n
a_21 a_22 ... a_2n
...
a_n1 a_n2 ... a_nn
*/

int main() {
	long double time_limit;
	int n;
	int J;
	cin>>n>>J>>time_limit;
	vector<vector<long double>> v;

	for (int i = 0; i < n; i++) {
		v.push_back(vector<long double>(n));
		for (int j = 0; j < n; j++) {
			cin>>v[i][j];
		}
	}
	chrono::steady_clock::time_point start_time = chrono::steady_clock::now();
	long double co = 0;
	long double su = NAN;
	cout<<setprecision(8)<<fixed;
	int pout = -1;
	long double ts = 0;
	while (time_since(start_time) < time_limit) {
		su = log_sum_exp(su, pps(v));
		co++;
		ts = time_since(start_time);
		if (isfinite(su) &&  (co < 1000 || pout != (int)(10 * ts))) {
			cout<<(su - log(co))<<" "<<ts<<endl;
			pout = (int)(10 * ts);
		}
	}
}