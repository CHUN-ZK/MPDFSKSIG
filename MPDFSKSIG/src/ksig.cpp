#include <mutex>
#include "header.h"

/*建立帕斯卡三角形表*/
void single_pascal_table(TABLE& pt, int N) {
	pt[0].distx[0] = 1;
	for (int i = 1; i < N; i++) {
		pt[i].distx[0] = 1;
		for (int t = 1; t < i + 1; t++) {
			pt[i].distx[t] = pt[i - 1].distx[t] + pt[i - 1].distx[t - 1];
		}
	}
}
/*建立帕斯卡三角形表(平行版)*/
void parallel_pascal_table(TABLE& pt, int N) {
	pt[0].distx[0] = 1;

	for (int i = 1; i < N; i++) {
		pt[i].distx[0] = 1;
		tbb::parallel_for(1, i + 1, [&](int t) {
			pt[i].distx[t] = pt[i - 1].distx[t] + pt[i - 1].distx[t - 1];
			});
	}
}

/*利用建立單一積項的K特徵值表*/
void single_ksig_table(TABLE& NCkt, TABLE pt, int N) {
	NCkt[0].distx[0] = 0;

	for (int i = 1; i < N; i++) {
		for (int t = 0; t < i; t++)
			NCkt[i].distx[t] = NCkt[i - 1].distx[t] * 2 + pt[i - 1].distx[t] * mp::pow(cpp_int(2), i - 1);
	}
}

/*利用建立單一積項的K特徵值表(平行版)*/
void parallel_ksig_table(TABLE& NCkt, TABLE pt, int N) {
	NCkt[0].distx[0] = 0;

	for (int i = 1; i < N; i++) {
		tbb::parallel_for(0, i, [&](int t) {
			NCkt[i].distx[t] = NCkt[i - 1].distx[t] * 2 + pt[i - 1].distx[t] * mp::pow(cpp_int(2), i - 1);
			});
	}
}

/*利用KSO分群*/
void kso_group(TABLE kso) {
	vector<int> check(NUMOUTPUTS, 0);
	vector<int> count(NUMOUTPUTS + 1, 0);
	int m = 0;
	for (int i = 0; i < kso.size(); i++) {
		m = 1;
		if (check[i] == 0) {
			for (int j = i + 1; j < kso.size(); j++) {
				if (kso[i] == kso[j]) {
					m++;
					check[j] = 1;
				}
			}
			check[i] = 1;
			count[m]++;
		}
	}

	for (int i = 1; i < count.size(); i++) {
		if (count[i] != 0) {
			std::cout << count[i] << "(" << i << ")";
		}
	}
	cout << " ";
}

/*利用KSI分群*/
void ksi_group(vector<TABLE> ksi) {
	vector<int> check(NUMINPUTS, 0);
	vector<int> count(NUMINPUTS + 1, 0);

	int m = 0;
	int match = 1;
	for (int v = 0; v < NUMINPUTS; v++) {
		if (check[v] == 0) {
			m = 1;
			check[v] = 1;
			for (int j = v + 1; j < NUMINPUTS; j++) {
				match = 1;
				for (int i = 0; i < NUMOUTPUTS && match; i++) {
					if (ksi[i][v] != ksi[i][j]) {
						match = 0;
					}
				}
				if (match) {
					check[j] = 1;
					m++;
				}
			}
			count[m]++;
		}
	}

	for (int i = 1; i < count.size(); i++) {
		if (count[i] != 0) {
			std::cout << count[i] << "(" << i << ")";
		}
	}
}

/*KSO與KSI計算*/
/*單執行緒循序(Naive)計算KSO*/
void sigle_kso(pset_family pf, ksig_value& kso, TABLE pt, TABLE kt) {
	pset one, two, r = set_new(cube.size);
	int dist, cup_num, cap_num, dc_num;
	int p, q, i, j;


	for (p = 0, one = pf->data; p < pf->count; p++, one += pf->wsize) {
		dc_num = dc_count(one);

		/*KSO*/
		kso += kt[dc_num];

		for (q = p + 1, two = pf->data + pf->wsize * q; q < pf->count; q++, two += pf->wsize) {
			dist = cdistN(one, two);
			if (dist == 0) continue;
			cup_num = set_ordc(one, two);
			cap_num = set_anddc(one, two);

			/*KSO*/
			for (i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
				kso.distx[i] += pt[cup_num].distx[j] * mp::pow(cpp_int(2), cap_num);
			}
		}
	}
	FREE(r);
}
/*單執行緒循序(Naive)計算KSI*/
void sigle_ksi(pset_family pf, TABLE& ksi, TABLE pt, TABLE kt, int N, pset sup) {
	pset one, two, r = set_new(cube.size);
	int dist, cup_num, cap_num, dc_num;
	int p, q, i, j, t, k;
	char c1, c2;

	for (int j = 0; j < NUMINPUTS; j++) {
		if (getvar(sup, j) == '?') {
			for (int i = 0; i < ksi[j].distx.size(); i++) {
				ksi[j].distx[i] = -1;
			}
		}
	}

	vector<cpp_int> dc_ksig, dist_ksig;
	cpp_int powdc2;

	dc_ksig.resize(N);
	dist_ksig.resize(N);

	for (p = 0, one = pf->data; p < pf->count; p++, one += pf->wsize) {
		dc_num = dc_count(r, one);

		/*KSI*/
		if (dc_num > 0) {
			for (i = 0; i < dc_num && i < N; i++) {
				dc_ksig[i] = pt[dc_num - 1].distx[i] * powdc2;
			}
			for (i = 1, t = 0; i <= cube.inword && t != dc_num; i++) {
				if (r[i]) {
					for (j = 0; j < varnum; j++) {
						if (is_in_set(r, (i - 1) * BPI + j * 2)) {
							for (k = 0; k < dc_num && k < N; k++) {
								ksi[(i - 1) * varnum + j].distx[k] += dc_ksig[k];
							}
							t++;
						}
					}
				}
			}
		}

		for (q = p + 1, two = pf->data + pf->wsize * q; q < pf->count; q++, two += pf->wsize) {
			dist = cdistN(one, two);
			if (dist == 0) continue;
			cup_num = set_ordc(one, two);
			cap_num = set_anddc(one, two);

			/*KSI*/
			if (dist <= N) {
				powdc2 = mp::pow(mp::cpp_int(2), cap_num);
				for (i = 0; i < N; i++) {
					if ((i + 1) > dist + cup_num || (i + 1) < dist) {

						dist_ksig[i] = 0;
						dc_ksig[i] = 0;
					}
					else if ((i + 1) - dist == 0) {
						dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
						dc_ksig[i] = 0;
					}
					else {
						dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
						dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * powdc2;
					}
				}

				for (i = 0; i < NUMINPUTS; i++) {
					c1 = getvar(one, i);
					c2 = getvar(two, i);

					if (c1 == '2' || c2 == '2') {
						for (j = dist; j <= N; j++) {
							ksi[i].distx[j - 1] += dc_ksig[j - 1];
						}
					}
					else if (c1 != c2) {
						for (j = dist; j <= N; j++) {
							ksi[i].distx[j - 1] += dist_ksig[j - 1];
						}
					}
				}
			}
		}
	}
	FREE(r);
}

/*Intel TBB \版本計算KSO*/
void tbb_kso(pset_family pf, ksig_value& kso, TABLE pt, TABLE kt) {
	spin_mutex sm;

	tbb::parallel_for(blocked_range<ullg>(0, pf->count, 1), [&](blocked_range<ullg> B) {
		ksig_value thread_kso(NUMINPUTS);
		pset one, two;
		int dist, cup_num, cap_num, dc_num;
		int i, j;
		cpp_int powdc2;
		for (ullg pid1 = B.begin(), p_end1 = B.end(); pid1 < p_end1; pid1++) {
			dc_num = dc_count(pf->data + pf->wsize * pid1);
			thread_kso += kt[dc_num];
			for (ullg pid2 = pid1; pid2 < pf->count; pid2++) {
				one = pf->data + pf->wsize * pid1;
				two = pf->data + pf->wsize * pid2;
				dist = cdistN(one, two);
				if (dist == 0) continue;
				cup_num = set_ordc(one, two);
				cap_num = set_anddc(one, two);

				powdc2 = mp::pow(mp::cpp_int(2), cap_num);

				for (i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
					thread_kso.distx[i] += pt[cup_num].distx[j] * powdc2;
				}
			}
		}
		{
			tbb::spin_mutex::scoped_lock lock(sm);
			kso += thread_kso;
		}
		});
}

/*Intel TBB版本計算KSI*/
void tbb_ksi(pset_family pf, TABLE& ksi, TABLE pt, TABLE kt, int N, pset sup) {
	for (int j = 0; j < NUMINPUTS; j++) {
		if (getvar(sup, j) == '?') {
			for (int i = 0; i < ksi[j].distx.size(); i++) {
				ksi[j].distx[i] = -1;
			}
		}
	}
	spin_mutex sm;
	tbb::parallel_for(blocked_range<ullg>(0, pf->count, 1), [&](blocked_range<ullg> B) {
		TABLE thread_ksi(NUMINPUTS, ksig_value(N));
		pset one, two, r = set_new(cube.size);
		int dist, cup_num, cap_num, dc_num;
		int i, j, t, k;
		char c1, c2;
		vector<cpp_int> dc_ksig(N), dist_ksig(N);
		cpp_int powdc2;
		for (ullg pid1 = B.begin(), p_end1 = B.end(); pid1 < p_end1; pid1++) {
			dc_num = dc_count(r, pf->data + pf->wsize * pid1);

			if (dc_num > 0) {
				powdc2 = mp::pow(mp::cpp_int(2), dc_num - 1);

				for (i = 0; i < dc_num && i < N; i++) {
					dc_ksig[i] = pt[dc_num - 1].distx[i] * powdc2;
				}

				for (i = 1, t = 0; i <= cube.inword && t != dc_num; i++) {
					if (r[i]) {
						for (j = 0; j < varnum; j++) {
							if (is_in_set(r, (i - 1) * BPI + j * 2)) {
								for (k = 0; k < N && k < dc_num; k++) {
									thread_ksi[(i - 1) * varnum + j].distx[k] += dc_ksig[k];
								}
								t++;
							}
						}
					}
				}
			}

			for (ullg pid2 = pid1; pid2 < pf->count; pid2++) {
				one = pf->data + pf->wsize * pid1;
				two = pf->data + pf->wsize * pid2;
				dist = cdistN(one, two);
				if (dist == 0) continue;
				cup_num = set_ordc(one, two);
				cap_num = set_anddc(one, two);

				if (dist <= N) {
					powdc2 = mp::pow(mp::cpp_int(2), cap_num);

					for (i = 0; i < N; i++) {
						if ((i + 1) > dist + cup_num || (i + 1) < dist) {

							dist_ksig[i] = 0;
							dc_ksig[i] = 0;
						}
						else if ((i + 1) - dist == 0) {
							dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * mp::pow(cpp_int(2), cap_num);
							dc_ksig[i] = 0;
						}
						else {
							dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * mp::pow(cpp_int(2), cap_num);
							dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * mp::pow(cpp_int(2), cap_num);
						}
					}

					for (i = 0; i < NUMINPUTS; i++) {
						c1 = getvar(one, i);
						c2 = getvar(two, i);
						for (j = dist; j <= N; j++) {
							if (c1 == '2' || c2 == '2') {
								thread_ksi[i].distx[j - 1] += dc_ksig[j - 1];
							}
							else if (c1 != c2) {
								thread_ksi[i].distx[j - 1] += dist_ksig[j - 1];
							}
						}
					}
				}
			}
		}
		FREE(r);
		{
			tbb::spin_mutex::scoped_lock lock(sm);
			for (int i = 0; i < N; i++) {
				ksi[i] += thread_ksi[i];
			}
		}
		});
}
/*OpenMP版本計算KSO*/
void omp_kso(pset_family pf, ksig_value& kso, TABLE pt, TABLE kt) {
	omp_lock_t olt;
	omp_init_lock(&olt);

#pragma omp parallel for schedule (dynamic)
	for (int p1 = 0; p1 < pf->count; p1++) {
		cpp_int powdc2;
		ksig_value thread_kso(NUMINPUTS);

		int dc_num = dc_count(pf->data + pf->wsize * p1);

		/*kso*/
		thread_kso += kt[dc_num];

		for (int p2 = p1 + 1; p2 < pf->count; p2++) {
			pset one, two;
			int dist, cup_num, cap_num;

			one = pf->data + pf->wsize * p1;
			two = pf->data + pf->wsize * p2;
			dist = cdistN(one, two);
			if (dist == 0) continue;
			cup_num = set_ordc(one, two);
			cap_num = set_anddc(one, two);
			powdc2 = mp::pow(mp::cpp_int(2), cap_num);

			/*kso*/
			for (int i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
				thread_kso.distx[i] += pt[cup_num].distx[j] * powdc2;
			}
		}
		{
			omp_set_lock(&olt);
			kso += thread_kso;
			omp_unset_lock(&olt);
		}
	}
	omp_destroy_lock(&olt);
}
/*OpenMP版本計算KSI*/
void omp_ksi(pset_family pf, TABLE& ksi, TABLE pt, TABLE kt, int N, pset sup) {
	for (int j = 0; j < NUMINPUTS; j++) {
		if (getvar(sup, j) == '?') {
			for (int i = 0; i < ksi[j].distx.size(); i++) {
				ksi[j].distx[i] = -1;
			}
		}
	}

	omp_lock_t olt;
	omp_init_lock(&olt);

#pragma omp parallel for schedule (dynamic)
	for (int p1 = 0; p1 < pf->count; p1++) {
		cpp_int powdc2;
		TABLE thread_ksi(NUMINPUTS, ksig_value(N));
		vector<cpp_int> dc_ksig(N), dist_ksig(N);
		int i, j, t, k;
		pset r = set_new(cube.size);
		int dc_num = dc_count(r, pf->data + pf->wsize * p1);

		/*ksi*/
		if (dc_num > 0) {
			powdc2 = mp::pow(mp::cpp_int(2), dc_num - 1);
			for (i = 0; i < dc_num && i < N; i++) {
				dc_ksig[i] = pt[dc_num - 1].distx[i] * powdc2;
			}
			for (i = 1, t = 0; i <= cube.inword && t != dc_num; i++) {
				if (r[i]) {
					for (j = 0; j < varnum; j++) {
						if (is_in_set(r, (i - 1) * BPI + j * 2)) {
							for (k = 0; k < N && k < dc_num; k++) {
								thread_ksi[(i - 1) * varnum + j].distx[k] += dc_ksig[k];
							}
							t++;
						}
					}
				}
			}
		}
		FREE(r);
		for (int p2 = p1 + 1; p2 < pf->count; p2++) {
			pset one, two;
			int dist, cup_num, cap_num;
			char c1, c2;

			one = pf->data + pf->wsize * p1;
			two = pf->data + pf->wsize * p2;
			dist = cdistN(one, two);
			if (dist == 0) continue;
			cup_num = set_ordc(one, two);
			cap_num = set_anddc(one, two);
			powdc2 = mp::pow(mp::cpp_int(2), cap_num);

			/*ksi*/
			if (dist <= N) {
				for (int i = 0; i < N; i++) {
					if ((i + 1) > dist + cup_num || (i + 1) < dist) {
						dist_ksig[i] = 0;
						dc_ksig[i] = 0;
					}
					else if ((i + 1) - dist == 0) {
						dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
						dc_ksig[i] = 0;
					}
					else {
						dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
						dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * powdc2;
					}
				}
				for (i = 0; i < NUMINPUTS; i++) {
					c1 = getvar(one, i);
					c2 = getvar(two, i);

					if (c1 == '2' || c2 == '2') {
						for (j = dist; j <= N; j++) {
							thread_ksi[i].distx[j - 1] += dc_ksig[j - 1];
						}
					}
					else if (c1 != c2) {
						for (j = dist; j <= N; j++) {
							thread_ksi[i].distx[j - 1] += dist_ksig[j - 1];
						}
					}
				}
			}
		}
		{
			omp_set_lock(&olt);
			for (int i = 0; i < N; i++) {
				ksi[i] += thread_ksi[i];
			}
			omp_unset_lock(&olt);
		}
	}
	omp_destroy_lock(&olt);
}

/*C++並行函式庫版本+DCP計算KSO*/
void chunk_kso(pset_family pf, ksig_value& kso, TABLE pt, TABLE kt, int tsize) {
	int chunksize = 10; /*資料塊大小*/
	int chunknum = ceil((double)pf->count / chunksize); /*計算on-set的數量為chunksize的幾倍，也為生成執行緒的數目*/
	int cs = ceil((double)pf->count / chunknum); /*計算平均一塊資料塊的積項數量*/

	int TN = chunknum > tsize ? tsize : chunknum;

	TABLE thread_kso(TN, ksig_value(NUMINPUTS));
	atomic<int> Index1 = 0, Index2 = 0;
	int job2_size = (chunknum - 1) * chunknum / 2;
	vector<pair<int, int>> Job2(job2_size, make_pair(-1, -1));

	for (int i = 0, index = 0; i < chunknum; i++) {
		for (int j = i + 1; j < chunknum; j++) {
			Job2[index] = make_pair(i, j);
			index++;
		}
	}

	vector<thread> thread_vec;
	for (int id = 0; id < TN; id++) {
		thread_vec.emplace_back(std::thread([&](int tid) {
			pset one, two;
			int dist, cup_num, cap_num, dc_num;
			int i, j, e1, e2, se1, se2;

			/*KSI*/
			cpp_int powdc2;

			int pid;

			while ((pid = Index1.fetch_add(1)) < chunknum) {
				for (e1 = pid * cs, e2 = 0, one = pf->data + pf->wsize * e1; e1 < pf->count && e2 < cs; e1++, e2++, one += pf->wsize) {
					dc_num = dc_count(one);

					/*KSO*/
					thread_kso[tid] += kt[dc_num];


					for (se1 = e1 + 1, se2 = e2 + 1, two = pf->data + pf->wsize * se1; se1 < pf->count && se2 < cs; se1++, se2++, two += pf->wsize) {
						dist = cdistN(one, two);
						cup_num = set_ordc(one, two);
						cap_num = set_anddc(one, two);
						powdc2 = mp::pow(mp::cpp_int(2), cap_num);

						if (dist > 0) {
							/*KSO*/
							for (i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
								thread_kso[tid].distx[i] += pt[cup_num].distx[j] * powdc2;
							}
						}
					}
				}
			}
			while ((pid = Index2.fetch_add(1)) < job2_size) {
				for (e1 = Job2[pid].first * cs, e2 = 0, one = pf->data + pf->wsize * e1; e1 < pf->count && e2 < cs; e1++, e2++, one += pf->wsize) {
					for (se1 = Job2[pid].second * cs, se2 = 0, two = pf->data + pf->wsize * se1; se1 < pf->count && se2 < cs; se1++, se2++, two += pf->wsize) {
						dist = cdistN(one, two);
						cup_num = set_ordc(one, two);
						cap_num = set_anddc(one, two);
						powdc2 = mp::pow(mp::cpp_int(2), cap_num);

						/*kso*/
						if (dist > 0) {
							/*KSO*/
							for (i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
								thread_kso[tid].distx[i] += pt[cup_num].distx[j] * powdc2;
							}
						}
					}
				}
			}
			}, id));
	}
	for (auto& thread : thread_vec) {
		thread.join();
	}
	for (int i = 0; i < TN; i++) {
		kso += thread_kso[i];
	}
}

/*C++並行函式庫版本+DCP計算KSI*/
void chunk_ksi(pset_family pf, TABLE& ksi, TABLE pt, TABLE kt, int N, pset sup, int tsize) {

	for (int j = 0; j < NUMINPUTS; j++) {
		if (getvar(sup, j) == '?') {
			for (int i = 0; i < ksi[j].distx.size(); i++) {
				ksi[j].distx[i] = -1;
			}
		}
	}

	int chunksize = 10; /*資料塊大小*/
	int chunknum = ceil((double)pf->count / chunksize); /*計算on-set的數量為chunksize的幾倍，也為生成執行緒的數目*/
	int cs = ceil((double)pf->count / chunknum); /*計算平均一塊資料塊的積項數量*/

	int TN = chunknum > tsize ? tsize : chunknum;

	atomic<int> Index1 = 0, Index2 = 0;
	int job2_size = (chunknum - 1) * chunknum / 2;
	vector<pair<int, int>> Job2(job2_size, make_pair(-1, -1));

	for (int i = 0, index = 0; i < chunknum; i++) {
		for (int j = i + 1; j < chunknum; j++) {
			Job2[index] = make_pair(i, j);
			index++;
		}
	}

	vector<TABLE> thread_ksi(TN, TABLE(NUMINPUTS, ksig_value(N)));
	vector<thread> thread_vec;
	for (int id = 0; id < TN; id++) {
		thread_vec.emplace_back(std::thread([&](int tid) {
			pset one, two, r = set_new(cube.size);
			int dist, cup_num, cap_num, dc_num;
			int i, j, t, k, e1, e2, se1, se2;
			char c1, c2;

			/*KSI*/
			vector<cpp_int> dc_ksig(N), dist_ksig(N);
			cpp_int powdc2;

			int pid;

			while ((pid = Index1.fetch_add(1)) < chunknum) {
				for (e1 = pid * cs, e2 = 0, one = pf->data + pf->wsize * e1; e1 < pf->count && e2 < cs; e1++, e2++, one += pf->wsize) {
					dc_num = dc_count(r, one);
					if (dc_num > 0) {
						/*KSO*/
						powdc2 = mp::pow(mp::cpp_int(2), dc_num - 1);
						for (i = 0; i < dc_num && i < N; i++) {
							dc_ksig[i] = pt[dc_num - 1].distx[i] * powdc2;
						}
						for (i = 1, t = 0; i <= cube.inword && t != dc_num; i++) {
							if (r[i]) {
								for (j = 0; j < varnum; j++) {
									if (is_in_set(r, (i - 1) * BPI + j * 2)) {
										for (k = 0; k < dc_num && k < N; k++) {
											thread_ksi[tid][(i - 1) * varnum + j].distx[k] += dc_ksig[k];
										}
										t++;
									}
								}
							}
						}
					}
					for (se1 = e1 + 1, se2 = e2 + 1, two = pf->data + pf->wsize * se1; se1 < pf->count && se2 < cs; se1++, se2++, two += pf->wsize) {
						dist = cdistN(one, two);
						cup_num = set_ordc(one, two);
						cap_num = set_anddc(one, two);
						powdc2 = mp::pow(mp::cpp_int(2), cap_num);

						if (dist > 0) {
							/*KSI*/
							if (dist <= N) {
								for (i = 0; i < N; i++) {
									if ((i + 1) > dist + cup_num || (i + 1) < dist) {
										dist_ksig[i] = 0;
										dc_ksig[i] = 0;
									}
									else if ((i + 1) - dist == 0) {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = 0;
									}
									else {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * powdc2;
									}
								}
								for (i = 0; i < NUMINPUTS; i++) {
									c1 = getvar(one, i);
									c2 = getvar(two, i);

									if (c1 == '2' || c2 == '2') {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dc_ksig[j - 1];
										}
									}
									else if (c1 != c2) {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dist_ksig[j - 1];
										}
									}
								}
							}
						}
					}
				}
			}
			FREE(r);
			while ((pid = Index2.fetch_add(1)) < job2_size) {
				for (e1 = Job2[pid].first * cs, e2 = 0, one = pf->data + pf->wsize * e1; e1 < pf->count && e2 < cs; e1++, e2++, one += pf->wsize) {
					for (se1 = Job2[pid].second * cs, se2 = 0, two = pf->data + pf->wsize * se1; se1 < pf->count && se2 < cs; se1++, se2++, two += pf->wsize) {
						dist = cdistN(one, two);
						cup_num = set_ordc(one, two);
						cap_num = set_anddc(one, two);
						powdc2 = mp::pow(mp::cpp_int(2), cap_num);

						if (dist > 0) {
							/*KSI*/
							if (dist <= N) {
								for (i = 0; i < N; i++) {
									if ((i + 1) > dist + cup_num || (i + 1) < dist) {
										dist_ksig[i] = 0;
										dc_ksig[i] = 0;
									}
									else if ((i + 1) - dist == 0) {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = 0;
									}
									else {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * powdc2;
									}
								}
								for (i = 0; i < NUMINPUTS; i++) {
									c1 = getvar(one, i);
									c2 = getvar(two, i);

									if (c1 == '2' || c2 == '2') {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dc_ksig[j - 1];
										}
									}
									else if (c1 != c2) {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dist_ksig[j - 1];
										}
									}
								}
							}
						}
					}
				}
			}
			}, id));
	}
	for (auto& thread : thread_vec) {
		thread.join();
	}
	for (int i = 0; i < TN; i++) {
		for (int j = 0; j < NUMINPUTS; j++) {
			ksi[j] += thread_ksi[i][j];
		}
	}
}

/*KSO與KSI同時計算*/
/*單執行緒循序(Naive)計算KSO與KSI*/
void single_cal_ksignature(pset_family pf, ksig_value& kso, TABLE& ksi, TABLE pt, TABLE kt, int N, pset sup) {

	pset one, two, r = set_new(cube.size);
	int dist, cup_num, cap_num, dc_num, dist_count, not_dist_count;
	int p, q, i, j, t, k;
	char c1, c2;

	for (int j = 0; j < NUMINPUTS; j++) {
		if (getvar(sup, j) == '?') {
			for (int i = 0; i < ksi[j].distx.size(); i++) {
				ksi[j].distx[i] = -1;
			}
		}
	}

	vector<cpp_int> dc_ksig, dist_ksig;
	cpp_int powdc2;

	dc_ksig.resize(N);
	dist_ksig.resize(N);

	for (p = 0, one = pf->data; p < pf->count; p++, one += pf->wsize) {

		/*單一積項*/
		dc_num = dc_count(r, one);
		if (dc_num > 0) {
			/*KSO*/
			kso += kt[dc_num];
			/*KSI*/
			powdc2 = mp::pow(mp::cpp_int(2), dc_num - 1);
			for (i = 0; i < dc_num && i < N; i++) {
				dc_ksig[i] = pt[dc_num - 1].distx[i] * powdc2;
			}
			for (i = 1, t = 0; i <= cube.inword && t != dc_num; i++) {
				if (r[i]) {
					for (j = 0; j < varnum; j++) {
						if (is_in_set(r, (i - 1) * BPI + j * 2)) {
							for (k = 0; k < dc_num && k < N; k++) {
								ksi[(i - 1) * varnum + j].distx[k] += dc_ksig[k];
							}
							t++;
						}
					}
				}
			}
		}

		for (q = p + 1, two = pf->data + pf->wsize * q; q < pf->count; q++, two += pf->wsize) {
			dist = cdistN(one, two);
			if (dist == 0) continue;
			cup_num = set_ordc(one, two);
			cap_num = set_anddc(one, two);

			/*KSO*/
			for (i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
				kso.distx[i] += pt[cup_num].distx[j] * mp::pow(cpp_int(2), cap_num);
			}

			/*KSI*/
			if (dist <= N) {
				powdc2 = mp::pow(mp::cpp_int(2), cap_num);
				for (i = 0; i < N; i++) {
					if ((i + 1) > dist + cup_num || (i + 1) < dist) {

						dist_ksig[i] = 0;
						dc_ksig[i] = 0;
					}
					else if ((i + 1) - dist == 0) {
						dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * mp::pow(cpp_int(2), cap_num);
						dc_ksig[i] = 0;
					}
					else {
						dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * mp::pow(cpp_int(2), cap_num);
						dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * mp::pow(cpp_int(2), cap_num);
					}
				}

				for (i = 0; i < NUMINPUTS; i++) {
					c1 = getvar(one, i);
					c2 = getvar(two, i);

					if (c1 == '2' || c2 == '2') {
						for (j = dist; j <= N; j++) {
							ksi[i].distx[j - 1] += dc_ksig[j - 1];
						}
					}
					else if (c1 != c2) {
						for (j = dist; j <= N; j++) {
							ksi[i].distx[j - 1] += dist_ksig[j - 1];
						}
					}
				}
			}
		}
	}
	FREE(r);
}

/*Intel TBB版本計算KSO與KSI*/
void tbb_cal_ksignature(pset_family pf, ksig_value& kso, TABLE& ksi, TABLE pt, TABLE kt, int N, pset sup) {
	for (int j = 0; j < NUMINPUTS; j++) {
		if (getvar(sup, j) == '?') {
			for (int i = 0; i < ksi[j].distx.size(); i++) {
				ksi[j].distx[i] = -1;
			}
		}
	}

	spin_mutex sm;
	tbb::parallel_for(blocked_range<ullg>(0, pf->count, 1), [&](blocked_range<ullg> B) {
		for (ullg pid1 = B.begin(), p_end1 = B.end(); pid1 < p_end1; pid1++) {
			ksig_value thread_kso(NUMINPUTS);
			TABLE thread_ksi(NUMINPUTS, ksig_value(N));

			pset one, two, r = set_new(cube.size);
			int dist, cup_num, cap_num, dc_num;
			int i, j, t, k;
			char c1, c2;
			vector<cpp_int> dc_ksig(N), dist_ksig(N);
			cpp_int powdc2;

			for (ullg pid2 = pid1; pid2 < pf->count; pid2++) {
				if (pid1 == pid2) {
					dc_num = dc_count(r, pf->data + pf->wsize * pid2);
					thread_kso += kt[dc_num];

					if (dc_num > 0) {
						powdc2 = mp::pow(mp::cpp_int(2), dc_num - 1);

						for (i = 0; i < N; i++) {
							if ((i + 1) > dc_num) {
								dc_ksig[i] = 0;
							}
							else {
								dc_ksig[i] = kt[dc_num].distx[i] * (i + 1) / dc_num;
							}
						}

						for (i = 1, t = 0; i <= cube.inword && t != dc_num; i++) {
							if (r[i]) {
								for (j = 0; j < varnum; j++) {
									if (is_in_set(r, (i - 1) * BPI + j * 2)) {
										for (k = 0; k < N; k++) {
											thread_ksi[(i - 1) * varnum + j].distx[k] += dc_ksig[k];
										}
										t++;
									}
								}
							}
						}
					}
				}
				else {
					one = pf->data + pf->wsize * pid1;
					two = pf->data + pf->wsize * pid2;
					dist = cdistN(one, two);
					if (dist == 0) continue;
					cup_num = set_ordc(one, two);
					cap_num = set_anddc(one, two);

					for (i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
						thread_kso.distx[i] += pt[cup_num].distx[j] * mp::pow(cpp_int(2), cap_num);
					}

					if (dist <= N) {
						powdc2 = mp::pow(mp::cpp_int(2), cap_num);

						for (i = 0; i < N; i++) {
							if ((i + 1) > dist + cup_num || (i + 1) < dist) {

								dist_ksig[i] = 0;
								dc_ksig[i] = 0;
							}
							else if ((i + 1) - dist == 0) {
								dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * mp::pow(cpp_int(2), cap_num);
								dc_ksig[i] = 0;
							}
							else {
								dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * mp::pow(cpp_int(2), cap_num);
								dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * mp::pow(cpp_int(2), cap_num);
							}
						}

						for (i = 0; i < NUMINPUTS; i++) {
							c1 = getvar(one, i);
							c2 = getvar(two, i);
							for (j = dist; j <= N; j++) {
								if (c1 == '2' || c2 == '2') {
									thread_ksi[i].distx[j - 1] += dc_ksig[j - 1];
								}
								else if (c1 != c2) {
									thread_ksi[i].distx[j - 1] += dist_ksig[j - 1];
								}
							}
						}
					}
				}
			}
			{
				tbb::spin_mutex::scoped_lock lock(sm);
				kso += thread_kso;
				for (i = 0; i < N; i++) {
					ksi[i] += thread_ksi[i];
				}
			}
			FREE(r);
		}
		});
}

/*OpenMP版本計算KSO與KSI*/
void omp_cal_ksignature(pset_family pf, ksig_value& kso, TABLE& ksi, TABLE pt, TABLE kt, int N, pset sup) {
	for (int j = 0; j < NUMINPUTS; j++) {
		if (getvar(sup, j) == '?') {
			for (int i = 0; i < ksi[j].distx.size(); i++) {
				ksi[j].distx[i] = -1;
			}
		}
	}

	omp_lock_t olt;
	omp_init_lock(&olt);

#pragma omp parallel for schedule (dynamic)
	for (int p1 = 0; p1 < pf->count; p1++) {
		cpp_int powdc2;
		ksig_value thread_kso(NUMINPUTS);
		TABLE thread_ksi(NUMINPUTS, ksig_value(N));
		vector<cpp_int> dc_ksig(N), dist_ksig(N);
		int i, j, t;
		pset r = set_new(cube.size);

		int dc_num = dc_count(r, pf->data + pf->wsize * p1);
		/*單一積項*/
		if (dc_num > 0) {
			/*KSO*/
			thread_kso += kt[dc_num];
			/*KSI*/
			powdc2 = mp::pow(mp::cpp_int(2), dc_num - 1);
			for (i = 0; i < dc_num && i < N; i++) {
				dc_ksig[i] = pt[dc_num - 1].distx[i] * powdc2;
			}
			for (i = 1, t = 0; i <= cube.inword && t != dc_num; i++) {
				if (r[i]) {
					for (j = 0; j < varnum; j++) {
						if (is_in_set(r, (i - 1) * BPI + j * 2)) {
							for (int k = 0; k < dc_num && k < N; k++) {
								thread_ksi[(i - 1) * varnum + j].distx[k] += dc_ksig[k];
							}
							t++;
						}
					}
				}
			}
		}
		FREE(r);
		for (int p2 = p1 + 1; p2 < pf->count; p2++) {
			pset one, two;
			int dist, cup_num, cap_num;
			char c1, c2;

			one = pf->data + pf->wsize * p1;
			two = pf->data + pf->wsize * p2;
			dist = cdistN(one, two);
			if (dist == 0) continue;
			cup_num = set_ordc(one, two);
			cap_num = set_anddc(one, two);
			powdc2 = mp::pow(mp::cpp_int(2), cap_num);

			/*kso*/
			for (int i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
				thread_kso.distx[i] += pt[cup_num].distx[j] * powdc2;
			}

			/*ksi*/
			if (dist <= N) {
				for (int i = 0; i < N; i++) {
					if ((i + 1) > dist + cup_num || (i + 1) < dist) {
						dist_ksig[i] = 0;
						dc_ksig[i] = 0;
					}
					else if ((i + 1) - dist == 0) {
						dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * mp::pow(cpp_int(2), cap_num);
						dc_ksig[i] = 0;
					}
					else {
						dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * mp::pow(cpp_int(2), cap_num);
						dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * mp::pow(cpp_int(2), cap_num);
					}
				}
				for (i = 0; i < NUMINPUTS; i++) {
					c1 = getvar(one, i);
					c2 = getvar(two, i);

					if (c1 == '2' || c2 == '2') {
						for (j = dist; j <= N; j++) {
							thread_ksi[i].distx[j - 1] += dc_ksig[j - 1];
						}
					}
					else if (c1 != c2) {
						for (j = dist; j <= N; j++) {
							thread_ksi[i].distx[j - 1] += dist_ksig[j - 1];
						}
					}
				}
			}
		}
		{
			omp_set_lock(&olt);
			kso += thread_kso;
			for (int i = 0; i < N; i++) {
				ksi[i] += thread_ksi[i];
			}
			omp_unset_lock(&olt);
		}
	}
	omp_destroy_lock(&olt);
}

/*C++並行函式庫版本+DCP計算KSO與KSI*/
void chunk_cal_ksignature(pset_family pf, ksig_value& kso, TABLE& ksi, TABLE pt, TABLE kt, int N, pset sup, int tsize) {
	
	for (int j = 0; j < NUMINPUTS; j++) {
		if (getvar(sup, j) == '?') {
			for (int i = 0; i < ksi[j].distx.size(); i++) {
				ksi[j].distx[i] = -1;
			}
		}
	}

	int chunksize = 10; /*資料塊大小*/
	int chunknum = ceil((double)pf->count / chunksize); /*計算on-set的數量為chunksize的幾倍，也為生成執行緒的數目*/
	int cs = ceil((double)pf->count / chunknum); /*計算平均一塊資料塊的積項數量*/

	int TN = chunknum > tsize ? tsize : chunknum;
	TABLE thread_kso(TN, ksig_value(NUMINPUTS));

	atomic<int> Index1 = 0, Index2 = 0;
	int job2_size = (chunknum - 1) * chunknum / 2;
	vector<pair<int, int>> Job2(job2_size, make_pair(-1, -1));

	for (int i = 0, index = 0; i < chunknum; i++) {
		for (int j = i + 1; j < chunknum; j++) {
			Job2[index] = make_pair(i, j);
			index++;
		}
	}

	vector<TABLE> thread_ksi(TN, TABLE(NUMINPUTS, ksig_value(N)));
	vector<thread> thread_vec;
	for (int id = 0; id < TN; id++) {
		thread_vec.emplace_back(std::thread([&](int tid) {
			pset one, two, r = set_new(cube.size);
			int dist, cup_num, cap_num, dc_num;
			int i, j, t, k, e1, e2, se1, se2;
			char c1, c2;

			/*KSI*/
			vector<cpp_int> dc_ksig(N), dist_ksig(N);
			cpp_int powdc2;

			int pid;

			while ((pid = Index1.fetch_add(1)) < chunknum) {
				for (e1 = pid * cs, e2 = 0, one = pf->data + pf->wsize * e1; e1 < pf->count && e2 < cs; e1++, e2++, one += pf->wsize) {
					dc_num = dc_count(r, one);
					if (dc_num > 0) {
						/*KSO*/
						thread_kso[tid] += kt[dc_num];
						powdc2 = mp::pow(mp::cpp_int(2), dc_num - 1);
						for (i = 0; i < dc_num && i < N; i++) {
							dc_ksig[i] = pt[dc_num - 1].distx[i] * powdc2;
						}
						for (i = 1, t = 0; i <= cube.inword && t != dc_num; i++) {
							if (r[i]) {
								for (j = 0; j < varnum; j++) {
									if (is_in_set(r, (i - 1) * BPI + j * 2)) {
										for (k = 0; k < dc_num && k < N; k++) {
											thread_ksi[tid][(i - 1) * varnum + j].distx[k] += dc_ksig[k];
										}
										t++;
									}
								}
							}
						}
					}
					for (se1 = e1 + 1, se2 = e2 + 1, two = pf->data + pf->wsize * se1; se1 < pf->count && se2 < cs; se1++, se2++, two += pf->wsize) {
						dist = cdistN(one, two);
						cup_num = set_ordc(one, two);
						cap_num = set_anddc(one, two);
						powdc2 = mp::pow(mp::cpp_int(2), cap_num);

						if (dist > 0) {
							/*KSO*/
							for (i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
								thread_kso[tid].distx[i] += pt[cup_num].distx[j] * powdc2;
							}
							/*KSI*/
							if (dist <= N) {
								for (i = 0; i < N; i++) {
									if ((i + 1) > dist + cup_num || (i + 1) < dist) {
										dist_ksig[i] = 0;
										dc_ksig[i] = 0;
									}
									else if ((i + 1) - dist == 0) {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = 0;
									}
									else {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * powdc2;
									}
								}
								for (i = 0; i < NUMINPUTS; i++) {
									c1 = getvar(one, i);
									c2 = getvar(two, i);

									if (c1 == '2' || c2 == '2') {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dc_ksig[j - 1];
										}
									}
									else if (c1 != c2) {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dist_ksig[j - 1];
										}
									}
								}
							}
						}
					}
				}
			}
			FREE(r);
			while ((pid = Index2.fetch_add(1)) < job2_size) {
				for (e1 = Job2[pid].first * cs, e2 = 0, one = pf->data + pf->wsize * e1; e1 < pf->count && e2 < cs; e1++, e2++, one += pf->wsize) {
					for (se1 = Job2[pid].second * cs, se2 = 0, two = pf->data + pf->wsize * se1; se1 < pf->count && se2 < cs; se1++, se2++, two += pf->wsize) {
						dist = cdistN(one, two);
						cup_num = set_ordc(one, two);
						cap_num = set_anddc(one, two);
						powdc2 = mp::pow(mp::cpp_int(2), cap_num);

						/*kso*/
						if (dist > 0) {
							/*KSO*/
							for (i = dist - 1, j = 0; j <= cup_num && i >= 0; i++, j++) {
								thread_kso[tid].distx[i] += pt[cup_num].distx[j] * powdc2;
							}
							/*KSI*/
							if (dist <= N) {
								for (i = 0; i < N; i++) {
									if ((i + 1) > dist + cup_num || (i + 1) < dist) {
										dist_ksig[i] = 0;
										dc_ksig[i] = 0;
									}
									else if ((i + 1) - dist == 0) {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = 0;
									}
									else {
										dist_ksig[i] = pt[cup_num].distx[(i + 1) - dist] * powdc2;
										dc_ksig[i] = pt[cup_num - 1].distx[(i + 1) - dist - 1] * powdc2;
									}
								}
								for (i = 0; i < NUMINPUTS; i++) {
									c1 = getvar(one, i);
									c2 = getvar(two, i);

									if (c1 == '2' || c2 == '2') {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dc_ksig[j - 1];
										}
									}
									else if (c1 != c2) {
										for (j = dist; j <= N; j++) {
											thread_ksi[tid][i].distx[j - 1] += dist_ksig[j - 1];
										}
									}
								}
							}
						}
					}
				}
			}
			}, id));
	}
	for (auto& thread : thread_vec) {
		thread.join();
	}
	for (int i = 0; i < TN; i++) {
		kso += thread_kso[i];
		for (int j = 0; j < NUMINPUTS; j++) {
			ksi[j] += thread_ksi[i][j];
		}
	}
}

/*用於比較四種版本計算KSO與KSI的時間*/
void cal_ksignature(pPLA PLA, int Mode) {
	clock_t start_time, single_time = 0, cpp_time = 0, tbb_time = 0, omp_time = 0, other_time = 0;

	int N = NUMINPUTS;

	vector<TABLE> KSO(4, TABLE(NUMOUTPUTS, ksig_value(NUMINPUTS)));
	vector<vector<TABLE>> KSI(4, vector<TABLE>(NUMOUTPUTS, TABLE(NUMINPUTS, ksig_value(N))));

	TABLE pt(NUMINPUTS, ksig_value(NUMINPUTS)), pt1(NUMINPUTS, ksig_value(NUMINPUTS));
	TABLE kt(NUMINPUTS + 1, ksig_value(NUMINPUTS)), kt1(NUMINPUTS + 1, ksig_value(NUMINPUTS));

	parallel_pascal_table(pt, NUMINPUTS);
	parallel_ksig_table(kt, pt, NUMINPUTS + 1);

	int tsize = thread::hardware_concurrency();
	
	pset_family f;
	int i;
	pset p, last;
	pset_family supp;
	if (NUMOUTPUTS == 1) {
		supp = sf_new(1, cube.size);
		sf_addset(supp, set_full(cube.size));
	}
	else {
		supp = support_set(PLA->F, PLA->R);
	}

	for (i = 0, p = supp->data; i < NUMOUTPUTS; i++, p += supp->wsize) {
		/* seperate function */
		{
			start_time = clock();
			f = sf_new(PLA->F->count, PLA->F->sf_size);
			f = sep_sup_output(PLA->F, i + cube.first_part[cube.output], p);
			f = sf_contain(f);
			f = make_disjoint(f);
			other_time += clock() - start_time;
		}

		/*KSIG-Output*/
		/*Single*/
		{
			start_time = clock();
			if (Mode == 0) {
				sigle_kso(f, KSO[0][i], pt, kt);
			}
			else if (Mode == 1) {
				sigle_ksi(f, KSI[0][i], pt, kt, N, p);
			}
			else {
				single_cal_ksignature(f, KSO[0][i], KSI[0][i], pt, kt, N, p);
			}
			single_time += clock() - start_time;
		}
		/*C++*/
		{
			start_time = clock();
			if (Mode == 0) {
				chunk_kso(f, KSO[1][i], pt, kt, tsize);
			}
			else if (Mode == 1) {
				chunk_ksi(f, KSI[1][i], pt, kt, N, p, tsize);
			}
			else {
				chunk_cal_ksignature(f, KSO[1][i], KSI[1][i], pt, kt, N, p, tsize);
			}
			cpp_time += clock() - start_time;
		}
		/*TBB*/
		{
			start_time = clock();
			if (Mode == 0) {
				tbb_kso(f, KSO[2][i], pt, kt);
			}
			else if (Mode == 1) {
				tbb_ksi(f, KSI[2][i], pt, kt, N, p);
			}
			else {
				tbb_cal_ksignature(f, KSO[2][i], KSI[2][i], pt, kt, N, p);
			}
			tbb_time += clock() - start_time;
		}
		/*OMP*/
		{
			start_time = clock();
			if (Mode == 0) {
				omp_kso(f, KSO[3][i], pt, kt);
			}
			else if (Mode == 1) {
				omp_ksi(f, KSI[3][i], pt, kt, N, p);
			}
			else {
				omp_cal_ksignature(f, KSO[3][i], KSI[3][i], pt, kt, N, p);
			}
			omp_time += clock() - start_time;
		}

		sf_free(f);
	}
	if (Mode == 0) {
		printf("KSO\n");
	}
	else if (Mode == 1) {
		printf("KSI\n");
	}
	else {
		printf("KSO & KSI\n");
	}

	printf("Naive time = %.3f \n", single_time / (double)CLOCKS_PER_SEC);
	printf("DCP time = %.3f \n", cpp_time / (double)CLOCKS_PER_SEC);
	printf("TBB time = %.3f \n", tbb_time / (double)CLOCKS_PER_SEC);
	printf("OpenMP time = %.3f \n", omp_time / (double)CLOCKS_PER_SEC);
	
	FREE(supp);
}

/*
	PLA = 需要計算電路的PLA
	mode = 選擇實作版本
	ioA = 選擇計算KSO或KSI, 0->KSO, 1->KSI, 2->KSO與KS
	I
*/
void ksio(pPLA PLA, int mode, int ioA) {
	clock_t start_time, total_time = 0;

	int N = NUMINPUTS;
	TABLE KSO(NUMOUTPUTS, ksig_value(NUMINPUTS));
	vector<TABLE> KSI(NUMOUTPUTS, TABLE(NUMINPUTS, ksig_value(N)));

	TABLE pt(NUMINPUTS, ksig_value(NUMINPUTS));
	TABLE kt(NUMINPUTS + 1, ksig_value(NUMINPUTS));

	pset_family f;
	int i;
	pset p, last;
	pset_family supp;
	if (NUMOUTPUTS == 1) {
		supp = sf_new(1, cube.size);
		sf_addset(supp, set_full(cube.size));
	}
	else {
		supp = support_set(PLA->F, PLA->R);
	}

	if (mode == 0) {
		single_pascal_table(pt, NUMINPUTS);
		single_ksig_table(kt, pt, NUMINPUTS + 1);

		for (i = 0, p = supp->data; i < NUMOUTPUTS; i++, p += supp->wsize) {
			f = sf_new(PLA->F->count, PLA->F->sf_size);
			f = sep_sup_output(PLA->F, i + cube.first_part[cube.output], p);
			f = sf_contain(f);
			f = make_disjoint(f);
			start_time = clock();
			if (ioA == 0) {
				sigle_kso(f, KSO[i], pt, kt);
			}
			else if (ioA == 1) {
				sigle_ksi(f, KSI[i], pt, kt, N, p);
			}
			else {
				single_cal_ksignature(f, KSO[i], KSI[i], pt, kt, N, p);
			}
			total_time += clock() - start_time;
		}
		std::printf("Naive ");
	}
	else if (mode == 1) {
		parallel_pascal_table(pt, NUMINPUTS);
		parallel_ksig_table(kt, pt, NUMINPUTS + 1);

		for (i = 0, p = supp->data; i < NUMOUTPUTS; i++, p += supp->wsize) {
			f = sf_new(PLA->F->count, PLA->F->sf_size);
			f = sep_sup_output(PLA->F, i + cube.first_part[cube.output], p);
			f = sf_contain(f);
			f = make_disjoint(f);
			start_time = clock();
			if (ioA == 0) {
				chunk_kso(f, KSO[i], pt, kt, 24);
			}
			else if (ioA == 1) {
				chunk_ksi(f, KSI[i], pt, kt, N, p, 24);
			}
			else {
				chunk_cal_ksignature(f, KSO[i], KSI[i], pt, kt, N, p, 24);
			}

			total_time += clock() - start_time;
		}
		std::printf("DCP ");
	}
	else if (mode == 2) {
		parallel_pascal_table(pt, NUMINPUTS);
		parallel_ksig_table(kt, pt, NUMINPUTS + 1);

		for (i = 0, p = supp->data; i < NUMOUTPUTS; i++, p += supp->wsize) {
			f = sf_new(PLA->F->count, PLA->F->sf_size);
			f = sep_sup_output(PLA->F, i + cube.first_part[cube.output], p);
			f = sf_contain(f);
			f = make_disjoint(f);
			start_time = clock();
			if (ioA == 0) {
				tbb_kso(f, KSO[i], pt, kt);
			}
			else if (ioA == 1) {
				tbb_ksi(f, KSI[i], pt, kt, N, p);
			}
			else {
				tbb_cal_ksignature(f, KSO[i], KSI[i], pt, kt, N, p);
			}
			total_time += clock() - start_time;
		}
		std::printf("TBB ");
	}
	else {
		parallel_pascal_table(pt, NUMINPUTS);
		parallel_ksig_table(kt, pt, NUMINPUTS + 1);

		for (i = 0, p = supp->data; i < NUMOUTPUTS; i++, p += supp->wsize) {
			f = sf_new(PLA->F->count, PLA->F->sf_size);
			f = sep_sup_output(PLA->F, i + cube.first_part[cube.output], p);
			f = sf_contain(f);
			f = make_disjoint(f);
			start_time = clock();
			if (ioA == 0) {
				omp_kso(f, KSO[i], pt, kt);
			}
			else if (ioA == 1) {
				omp_ksi(f, KSI[i], pt, kt, N, p);
			}
			else {
				omp_cal_ksignature(f, KSO[i], KSI[i], pt, kt, N, p);
			}
			total_time += clock() - start_time;
		}
		std::printf("OpenMP ");
	}
	if (ioA == 0) {
		std::printf("KSO = ");
	}
	else if (ioA == 1) {
		std::printf("KSI = ");
	}
	else {
		std::printf("KSO & KSI = ");
	}

	std::printf("%.3f ", total_time / (double)CLOCKS_PER_SEC);
	FREE(supp);
}