#include <string>
#include "header.h"
#include "main.h"

int     opterr = 1,      /* if error message should be printed */
optind = 1,             /* index into parent argv vector */
optopt,                 /* character checked for validity */
optreset;               /* reset getopt */
string optarg, optarg2;                /* argument associated with option */

#define BADCH   (int)'?'
#define BADARG  (int)':'
#define EMSG    ""

/*Ū���Ѽ�*/
int getopt(int nargc, char* const nargv[], const char* ostr)
{
	static char* place = (char*)EMSG;              /* option letter processing */
	const char* oli;                        /* option letter list index */

	if (optreset || !*place) {              /* update scanning pointer */
		optreset = 0;
		if (optind >= nargc || *(place = nargv[optind]) != '-') {
			place = (char*)EMSG;
			return (-1);
		}
		if (place[1] && *++place == '-') {      /* found "--" */
			++optind;
			place = (char*)EMSG;
			return (-1);
		}
	}                                       /* option letter okay? */
	if ((optopt = (int)*place++) == (int)':' || (optopt == '.') || !(oli = strchr(ostr, optopt))) {
		/*
		* if the user didn't specify '-' as an option,
		* assume it means -1.
		*/
		if (optopt == (int)'-')
			return (-1);
		if (!*place)
			++optind;
		if (opterr && *ostr != ':' && *ostr != '.')
			(void)printf("illegal option -- %c\n", optopt);
		return (BADCH);
	}

	if (*++oli == ':') {
		if (*place)                     /* no white space */
			optarg = place;
		else if (nargc <= ++optind) {   /* no arg */
			place = (char*)EMSG;
			if (*ostr == ':')
				return (BADARG);
			if (opterr)
				(void)printf("option requires an argument -- %c\n", optopt);
			return (BADCH);
		}
		else                            /* white space */
			optarg = nargv[optind];
		place = (char*)EMSG;
		++optind;
	}
	else if (*oli == '.') {
		if (*place)                     /* no white space */
			optarg = place;
		else if (nargc <= (optind++) + 2) {   /* no arg */
			place = (char*)EMSG;
			if (*ostr == '.')
				return (BADARG);
			if (opterr)
				(void)printf("option requires an argument -- %c\n", optopt);
			return (BADCH);
		}
		else {							/* white space */
			optarg = nargv[optind];
			optarg2 = nargv[++optind];
		}
		place = (char*)EMSG;
		optind++;
	}
	else {
		optarg = "";
		if (!*place)
			++optind;
	}

	return (optopt);                        /* dump back option letter */
}

/*Ū��PLA��*/
void get_PLA(char* PLA_address, int pla_type, pPLA* PLA) {
	int input_type = FD_type;
	int out_type = F_type;
	FILE* fp;
	int needs_dcset = TRUE;
	int needs_offset = TRUE;
	int planum = 1;

	if (pla_type == FR_type) {
		input_type = FR_type; //default
		needs_offset = FALSE;
		needs_dcset = FALSE;

	}
	if ((fp = fopen(PLA_address, "r")) == NULL) {
		fprintf(stderr, "�L�k�}�� %s\n", PLA_address);
		exit(1);
	}
	if (read_pla(fp, needs_dcset, needs_offset, input_type, PLA) == EOF) {
		fprintf(stderr, "�L�k���%s����m\n", PLA_address);
		exit(1);
	}

	(*PLA)->filename = _strdup(PLA_address);
	filename = (*PLA)->filename;
}

void general_PLA(pPLA& pla, char* filename, int sym, int input, int output, int onset_size, int offset_size, int hd3) {

	general_setup(input, output);

	srand(time(NULL));

	FILE* output_file = fopen(filename, "w");
	pla = new_PLA();

	pla->F = sf_new(onset_size, input * 2 + 1);
	pla->R = sf_new(offset_size, input * 2 + 1);

	if (sym == SYM) { //�ͦ��H���q��
		GeneralPLA_allrandom(pla->F, pla->R, input, onset_size, offset_size, hd3);
	}
	else { //�ͦ��S����ƹ�٩ʹq��

		while (!GeneralPLA_nonSym(pla->F, pla->R, input, onset_size, offset_size, hd3)) {
			sf_free(pla->F);
			sf_free(pla->R);
			pla->F = sf_new(onset_size, input * 2 + 1);
			pla->R = sf_new(offset_size, input * 2 + 1);
		}
		GeneralPLA_allrandom(pla->F, pla->R, input, onset_size, offset_size, hd3);
		pla->F = sf_rand(pla->F);
		pla->R = sf_rand(pla->R);
	}
	fprint_pla(output_file, pla, FR_type);
	fclose(output_file);
}

void test_bigtesting(pset_family& F, pset_family& R, int c1, int c2, int ordering, int t_) {

	tbb_sym(F, R, c1, c2);
	omp_sym(F, R);
	if (ordering_table[ordering].method == WEIGHT) {
		sym_testing(F, R, c1, c2, ordering, t_);
	}
	else {
		cpp_dcp(F, R, c1, c2, 1, t_);
	}
	sym_naive(F, R);
	cout << endl;
}

int trans_strnum(string num) {
	if (num == "0") {
		return 0;
	}
	for (int i = 0; i < num.size(); i++) {
		if (!isdigit(num[i])) {
			return -1;
		}
	}

	return stoi(num);
}

int main(int argc, char** argv) {
	int i, j;

	int option = 0; /* default -D: ��ٰt���ˬd */
	int pla_type = FD_type; /* default -R: CSF(completely specify function) */

	int method = DCP; /* default -M: ��@��k������*/

	/*����ˬd*/
	int tsize = thread::hardware_concurrency(); /* default -t: �̤j�i�ΰ�����ƶq */
	int chunk_row = 200; /* -c: chunk_row chunk_col , default 200 X 200 */
	int chunk_col = 200;
	int ordering = 1; /* default -o: �H�צV���覡�N��ƶ���J�u�@��C */

	/*K-�S�x�ȭp��*/
	int KSIG = ALL; /* default -K: �p��KSO�PKSI */

	/*�ͦ��q��*/
	int general_pla = FALSE; /* default -G: ���ϥΥͦ����q���i����� */
	int sym = SYM; /* default -S: �ͦ�����٪��q��*/
	int input_size = 1000; /* -i: �ͦ��q������J�ƶq */
	int output_size = 1; /* �ͦ��q������X�ƶq�T�w��1 */
	int onset_size = 10000; /*  -f: �ͦ��q����onset�n���ƶq */
	int offset_size = 10000; /*  -r: �ͦ��q����offset�n���ƶq */
	int hd3 = 1000; /* -h: �ͦ��q����HD�p��3���n���t��ƶq */

	verbose = TRUE; /* default -v: ��X�ثe���窺detail */

	/*-------------------------------���O�P�w-------------------------------*/
	while ((i = getopt(argc, argv, "D:RM:t:c.o:K:GSi:f:r:h:v")) != EOF) {
		switch (i) {
		case 'D':		/* -D ��ܤl�R�O�A�i�ѦҦ�main.h����option_table */
			for (j = 0; j < size(option_table); j++) {
				if (optarg == option_table[j].name) {
					option = j;
					break;
				}
			}
			if (j == size(option_table)) {
				fprintf(stderr, "%s: �L��option���O \"%s\"\n", argv[0], optarg.c_str());
				exit(1);
			}
			break;

		case 'R':		/* -R ��JPLA�榡 */
			pla_type = FR_type; /* ��JPLA��FR_type */
			break;

		case 'M':		/* -M ��@��k */
			for (j = 0; j < size(implements_table); j++) {
				if (optarg == implements_table[j].name) {
					method = j;
					break;
				}
			}
			if (j == size(implements_table)) {
				fprintf(stderr, "%s: �L�Ī���@��k \"%s\"\n", argv[0], optarg.c_str());
				exit(1);
			}
			break;

		case 't':		/* -t �̤j�i�ΰ�����ƶq */
			tsize = trans_strnum(optarg);
			if (tsize <= 0) {
				fprintf(stderr, "%s: ������ƶq�u�ର����� \"%s\"\n", argv[0], optarg.c_str());
				exit(1);
			}
			break;
		case 'c':		/* -c chunksize = chunk_rowXchunk_col */
			chunk_row = trans_strnum(optarg);
			chunk_col = trans_strnum(optarg2);
			if (chunk_row <= 0 || chunk_col <= 0) {
				fprintf(stderr, "%s: chunksize�u�ର����� \"%s %s\"\n", argv[0], optarg.c_str(), optarg2.c_str());
				exit(1);
			}
			break;

		case 'o':		/* -o ��ƶ���J�u�@��C�覡 */
			for (j = 0; j < size(ordering_table); j++) {
				if (optarg == ordering_table[j].name) {
					ordering = j;
					break;
				}
			}
			if (j == size(ordering_table)) {
				fprintf(stderr, "%s: �L�Ī���ƶ���J�u�@��C��k \"%s\"\n", argv[0], optarg.c_str());
				exit(1);
			}
			break;

		case 'K':		/* -K �p��K-�S�x�Ȫ�KSO��KSI�Υ��� */
			for (j = 0; j < size(ksig_table); j++) {
				if (optarg == ksig_table[j].name) {
					KSIG = ksig_table[j].k_mode;
					break;
				}
			}
			if (j == size(ksig_table)) {
				fprintf(stderr, "%s: �L�Ī�K-�S�x�� \"%s\"\n", argv[0], optarg.c_str());
				exit(1);
			}
			break;

		case 'G':		/* -G �ϥΥͦ����q���i����� */
			general_pla = TRUE;
			break;

		case 'S':		/* -S �ͦ��S����ƹ�٩ʪ��q���i����� */
			sym = NSYM;
			break;

		case 'i':		/* -i �ͦ��q������J�ƶq */
			input_size = trans_strnum(optarg);
			if (input_size <= 0) {
				fprintf(stderr, "%s: �ͦ��q������J�ƶq�u�ର����� \"%s\"\n", argv[0], optarg.c_str());
				exit(1);
			}
			break;

		case 'f':		/* -f �ͦ��q����onset�n���ƶq */
			onset_size = trans_strnum(optarg);
			if (onset_size <= 0) {
				fprintf(stderr, "%s: �ͦ��q����onset�n���ƶq�u�ର����� \"%s\"\n", argv[0], optarg.c_str());
				exit(1);
			}
			break;

		case 'r':		/* -r �ͦ��q����offset�n���ƶq */
			offset_size = trans_strnum(optarg);

			if (offset_size <= 0) {
				fprintf(stderr, "%s: �ͦ��q����offset�n���ƶq�u�ର����� \"%s\"\n", argv[0], optarg.c_str());
				exit(1);
			}
			break;

		case 'h':		/* -h �ͦ��q����HD�p��3���n���t��ƶq */
			hd3 = trans_strnum(optarg);
			if (hd3 <= 0) {
				fprintf(stderr, "%s: �ͦ��q����HD�p��3���n���t��ƶq�u�ର����� \"%s\"\n", argv[0], optarg.c_str());
				exit(1);
			}
			break;

		case 'v':		/* -v ������X��detail*/
			verbose = FALSE;
			break;

		default:
			exit(1);
		}
	}

	if ((general_pla || pla_type == FR_type) && option_table[option].name == "ksig") {
		fprintf(stderr, "K-�S�x�ȭp�⤣�i�ϥΥͦ����q��\n");
		exit(1);
	}

	if (optind++ >= argc) {
		fprintf(stderr, "�L��J�q��\n");
		exit(1);
	}

	/*-------------------------------Ū���Υͦ�PLA-------------------------------*/

	pPLA PLA = NIL(PLA_t);
	if (general_pla) {
		if (verbose)
			printf("GENERAL PLA,\nSYM = %s, \nInput = %d, \nOuput = 1, \nOnset size = %d, \nOffset size = %d, \nHD<3 size  = %d, \n", sym ? "Sym" : "Non-Sym", input_size, onset_size, offset_size, hd3);
		general_PLA(PLA, argv[optind - 1], sym, input_size, 1, onset_size, offset_size, hd3);
	}
	else {
		get_PLA(argv[optind - 1], pla_type, &PLA);
	}

	/*-------------------------------operations-------------------------------*/
	switch (option_table[option].key) {

	case KEY_SYMMETRY:

		switch (implements_table[method].key) {

		case ALL:
			if (verbose)
				printf("TESTING ALL SYMMETRY IMPLEMENTS,\nchunk_row = %d,\nchunk_col = %d, \nordering = %s, \nthread_size = %d,\n", chunk_row, chunk_col, ordering_table[ordering].full_name.c_str(), tsize);

			test_bigtesting(PLA->F, PLA->R, chunk_row, chunk_col, ordering, tsize);
			break;

		case NAIVE:
			if (verbose)
				printf("TESTING NAIVE SYMMETRY\n");
			sym_naive(PLA->F, PLA->R);
			break;

		case TBB:
			if (verbose)
				printf("TESTING Intel TBB SYMMETRY, \nchunk_row = %d,\nchunk_col = %d \n", chunk_row, chunk_col);

			tbb_sym(PLA->F, PLA->R, chunk_row, chunk_col);
			break;

		case OPENMP:
			if (verbose)
				printf("TESTING OpenMP SYMMETRY\n");

			omp_sym(PLA->F, PLA->R);
			break;

		case DCP:
			if (verbose)
				printf("TESTING DCP SYMMETRY, \nchunk_row = %d,\nchunk_col = %d, \nordering = %s, \nthread_size = %d,\n", chunk_row, chunk_col, ordering_table[ordering].full_name.c_str(), tsize);
			if (ordering_table[ordering].method == WEIGHT) {
				sym_testing(PLA->F, PLA->R, chunk_row, chunk_col, ordering, tsize);
			}
			else {
				cpp_dcp(PLA->F, PLA->R, chunk_row, chunk_col, ordering, tsize);
			}
			break;
		}
		break;
	case KEY_KSIG:
		if (implements_table[method].key == ALL) {
			if (verbose)
				printf("TESTING ALL IMPLEMENTS %s\n", ksig_table[KSIG].full_name.c_str());
			cal_ksignature(PLA, KSIG);
		}
		else {
			if (verbose)
				printf("TESTING %s IMPLEMENTS %s\n", implements_table[method].full_name.c_str(), ksig_table[KSIG].full_name.c_str());
			ksio(PLA, implements_table[method].id, KSIG);
		}
		break;

	case KEY_TEST_SYM:
		if (verbose)
			printf("�ϥ�Naive�BIntel TBB�BOpenMP�PDCP(unsort & sort)�i�����ˬd����, \nchunk_row = %d,\nchunk_col = %d, \nordering = %s, \nthread_size = %d,\n", chunk_row, chunk_col, ordering_table[ordering].full_name.c_str(), tsize);
		symmetry_benchmarking(PLA->F, PLA->R, chunk_row, chunk_col, ordering, tsize);
		break;

	case KEY_TEST_LARGECIRCUIT:
		if (verbose)
			printf("�ϥ�Naive�BIntel TBB�BOpenMP�PDCP�i�����ˬd����, \nchunk_row = %d,\nchunk_col = %d, \nordering = %s, \nthread_size = %d,\n", chunk_row, chunk_col, ordering_table[ordering].full_name.c_str(), tsize);
		symmetry_LargeCircuit(PLA->F, PLA->R, chunk_row, chunk_col, ordering, tsize);
		break;

	case KEY_TEST_KSIG:
		if (verbose)
			printf("�ϥ�Naive�BIntel TBB�BOpenMP�PDCP�i��K-�S�x�ȹ���, \nchunk_row = %d,\nchunk_col = %d, \nordering = %s, \nthread_size = %d,\n", chunk_row, chunk_col, ordering_table[ordering].full_name.c_str(), tsize);
		cal_ksignature(PLA, KSIG);
		break;

	case KEY_CHUNK:
		if (verbose)
			printf("�ϥΤ��Pchunksize�i�����ˬd����, \nordering = %s, \nthread_size = %d,\n", ordering_table[ordering].full_name.c_str(), tsize);

		testing_chunksize(PLA->F, PLA->R, ordering, tsize);
		break;

	case KEY_THREAD:
		if (verbose)
			printf("�ϥΤ��P������ƶq�i�����ˬd����, \nchunk_row = %d,\nchunk_col = %d, \nordering = %s,\n", chunk_row, chunk_col, ordering_table[ordering].full_name.c_str());

		testing_threadnum(PLA->F, PLA->R, chunk_row, chunk_col, ordering);
		break;

	case KEY_ORDERING:
		if (verbose)
			printf("�ϥΤ��P��ƶ���J�u�@��C����k�i�����ˬd����, \nchunk_row = %d,\nchunk_col = %d, \nthread_size = %d,\n", chunk_row, chunk_col, tsize);

		testing_ordering(PLA->F, PLA->R, chunk_row, chunk_col, tsize);
		break;

	case KEY_SYMRATE:
		if (verbose)
			printf("�ϥγ��������ƶ��i���٩ʤU������ˬd, \nchunk_row = %d,\nchunk_col = %d, \nthread_size = %d,\n", chunk_row, chunk_col, 1);
		testing_symrate(PLA->F, PLA->R, chunk_row, chunk_col, ordering);
		break;
	}

	free_PLA(PLA);

	return 0;
}